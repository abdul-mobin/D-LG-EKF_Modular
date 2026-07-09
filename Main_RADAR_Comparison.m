%% 0. Setup paths
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(genpath(repoRoot));

%% 1. Simulation Configuration
rng(42);
config = struct('num_runs', 50, 'azimuth_stds', 0:2:10, 'dt', 0.01, 'dt_radar', 0.1, 't_end', 200);
config.N = round(config.t_end / config.dt) + 1;
config.meas_steps = round(config.dt_radar / config.dt);

% Noise specs
noise.std_v = 0.03; noise.std_omega = deg2rad(0.2);
noise.Q_sys = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);

results_ekf = zeros(length(config.azimuth_stds), 1);
results_lg = zeros(length(config.azimuth_stds), 1);

%% 2. Monte Carlo Sweep
for az_idx = 1:length(config.azimuth_stds)
    current_az_std = deg2rad(config.azimuth_stds(az_idx));
    sum_err = [0, 0]; % [EKF, LG]
    
    for run = 1:config.num_runs
        disp(['Current SD index is' num2str(az_idx) ':: Run index is' num2str(run)]);
        % Generate Truth
        true_state = generate_ground_truth(config, noise);
        
        % Initialize Filters
        [state_ekf, P_ekf] = init_ekf(true_state(:,1));
        [g_lgkf, P_lgkf] = init_lgkf(true_state(:,1));
        
        err_run = [0, 0];
        
        for k = 2:config.N
            % PREDICTION
            [state_ekf, P_ekf] = predict_ekf(state_ekf, P_ekf, noise.Q_sys, config.dt);
            [g_lgkf, P_lgkf]   = predict_lgkf(g_lgkf, P_lgkf, noise.Q_sys, config.dt);
            
            % UPDATE
            if mod(k, config.meas_steps) == 0
                z = get_radar_measurement(true_state(:,k), current_az_std);
                
                [state_ekf, P_ekf] = update_ekf(state_ekf, P_ekf, z, current_az_std);
                [g_lgkf, P_lgkf]   = update_lgkf(g_lgkf, P_lgkf, z, current_az_std);
            end
            
            % Error Accumulation
            lg_vec = lie2vec_radar(g_lgkf);
            err_run(1) = err_run(1) + (true_state(1,k) - state_ekf(1))^2 + (true_state(2,k) - state_ekf(2))^2;
            err_run(2) = err_run(2) + (true_state(1,k) - lg_vec(1))^2 + (true_state(2,k) - lg_vec(2))^2;
        end
        sum_err = sum_err + (err_run ./ config.N);
    end
    results_ekf(az_idx) = sqrt(sum_err(1) / config.num_runs);
    results_lg(az_idx) = sqrt(sum_err(2) / config.num_runs);
end

plot_results(config.azimuth_stds, results_ekf, results_lg);