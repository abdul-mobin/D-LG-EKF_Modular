%% 0. Setup paths
repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(genpath(repoRoot));

%% 1. Simulation Configuration
rng(42);
config = struct('num_runs', 50, 'azimuth_stds', 0:2:10, 'dt', 0.01, 'dt_radar', 0.1, 't_end', 200);
config.N = round(config.t_end / config.dt) + 1;
config.meas_steps = round(config.dt_radar / config.dt);

% Noise specs and process-noise tuning
noise.std_v = 0.03; noise.std_omega = deg2rad(0.2);

% Separate tuning knobs for the EKF and LGKF process-noise covariance.
% These scale the nominal velocity/turn-rate noise terms independently.
config.q_tuning = struct( ...
    'ekf', struct('v', 1.0, 'omega', 1.0), ...
    'lgkf', struct('v', 1.0, 'omega', 1.0));

% Try to load tuned Q matrices from disk when available; otherwise fall back to
% the default manual scaling above.
matfile = fullfile(repoRoot, 'tuned_Q_matrices_2param.mat');
if exist(matfile, 'file')
    loaded = load(matfile);
    noise.Q_sys_ekf = loaded.Q_ekf_tuned;
    noise.Q_sys_lgkf = loaded.Q_lgkf_tuned;
else
    noise.Q_sys_ekf = diag([0, 0, 0, ...
        (config.q_tuning.ekf.v * noise.std_v)^2, ...
        (config.q_tuning.ekf.omega * noise.std_omega)^2]);

    noise.Q_sys_lgkf = diag([0, 0, 0, ...
        (config.q_tuning.lgkf.v * noise.std_v)^2, ...
        (config.q_tuning.lgkf.omega * noise.std_omega)^2]);
end

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
            [state_ekf, P_ekf] = predict_ekf(state_ekf, P_ekf, noise.Q_sys_ekf, config.dt);
            [g_lgkf, P_lgkf]   = predict_lgkf(g_lgkf, P_lgkf, noise.Q_sys_lgkf, config.dt);
            
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