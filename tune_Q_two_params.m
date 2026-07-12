function tune_Q_two_params()
    % Tune the process-noise covariance Q separately for the EKF and LGKF.
    %
    % This script uses the repository's existing filter implementations,
    % trajectory generation, and measurement model instead of a standalone
    % toy example.

    repoRoot = fileparts(mfilename('fullpath'));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    % Nominal process-noise values used to build the base covariance.
    sigma_v = 0.03;
    sigma_omega = deg2rad(0.2);
    Q_true = diag([0, 0, 0, sigma_v^2, sigma_omega^2]);

    % Simulation settings for tuning.
    config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', 200);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    noise.std_v = sigma_v;
    noise.std_omega = sigma_omega;

    % Use a fixed azimuth noise level for tuning.
    az_std = deg2rad(5);

    % Generate a small set of trajectories for optimization.
    rng(42);
    n_traj = 20;
    trajectories = cell(n_traj, 1);
    for i = 1:n_traj
        trajectories{i} = generate_ground_truth(config, noise);
    end

    opts = optimset('Display', 'iter', 'TolX', 1e-3, 'TolFun', 1e-3);

    % Optimize over log-scale factors for the velocity and turn-rate noise.
    x0 = [0, 0];

    objective_ekf = @(x) rmse_for_two_scales(x, Q_true, config, az_std, trajectories, 'ekf');
    x_ekf = fminsearch(objective_ekf, x0, opts);
    Q_ekf_tuned = scale_Q(Q_true, x_ekf);

    objective_lgkf = @(x) rmse_for_two_scales(x, Q_true, config, az_std, trajectories, 'lgkf');
    x_lgkf = fminsearch(objective_lgkf, x0, opts);
    Q_lgkf_tuned = scale_Q(Q_true, x_lgkf);

    fprintf('EKF  : scale_v = %.4f, scale_omega = %.4f\n', exp(x_ekf(1)), exp(x_ekf(2)));
    fprintf('LGKF : scale_v = %.4f, scale_omega = %.4f\n', exp(x_lgkf(1)), exp(x_lgkf(2)));

    save('tuned_Q_matrices_2param.mat', 'Q_ekf_tuned', 'Q_lgkf_tuned');
end

function Q = scale_Q(Q_base, log_scales)
    % Apply multiplicative scaling to the velocity and turn-rate noise entries.
    Q = Q_base;
    Q(4,4) = Q_base(4,4) * exp(log_scales(1));
    Q(5,5) = Q_base(5,5) * exp(log_scales(2));
end

function rmse = rmse_for_two_scales(log_scales, Q_base, config, az_std, trajectories, filter_type)
    Q = scale_Q(Q_base, log_scales);
    n_traj = numel(trajectories);
    errs = zeros(n_traj, 1);

    for i = 1:n_traj
        traj = trajectories{i};
        if strcmp(filter_type, 'ekf')
            est_xy = run_ekf_tuning(traj, Q, config, az_std);
        else
            est_xy = run_lgkf_tuning(traj, Q, config, az_std);
        end

        % Compare estimated and true position over time.
        pos_err = sqrt(sum((est_xy - traj(1:2, :)').^2, 2));
        errs(i) = sqrt(mean(pos_err.^2));
    end

    rmse = mean(errs);
end

function est_xy = run_ekf_tuning(true_state, Q, config, az_std)
    [x, P] = init_ekf(true_state(:, 1));
    est_xy = zeros(config.N, 2);
    est_xy(1, :) = x(1:2).';

    for k = 2:config.N
        [x, P] = predict_ekf(x, P, Q, config.dt);
        if mod(k, config.meas_steps) == 0
            z = get_radar_measurement(true_state(:, k), az_std);
            [x, P] = update_ekf(x, P, z, az_std);
        end
        est_xy(k, :) = x(1:2).';
    end
end

function est_xy = run_lgkf_tuning(true_state, Q, config, az_std)
    [g, P] = init_lgkf(true_state(:, 1));
    est_xy = zeros(config.N, 2);
    lg_vec = lie2vec_radar(g);
    est_xy(1, :) = lg_vec(1:2).';

    for k = 2:config.N
        [g, P] = predict_lgkf(g, P, Q, config.dt);
        if mod(k, config.meas_steps) == 0
            z = get_radar_measurement(true_state(:, k), az_std);
            [g, P] = update_lgkf(g, P, z, az_std);
        end
        lg_vec = lie2vec_radar(g);
        est_xy(k, :) = lg_vec(1:2).';
    end
end
