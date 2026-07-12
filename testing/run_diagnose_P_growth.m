function run_diagnose_P_growth(filter_type)
    % Generate the state/covariance histories needed by diagnose_P_growth.
    %
    % Usage:
    %   run('testing/run_diagnose_P_growth')           % EKF run
    %   run('testing/run_diagnose_P_growth', 'lgkf')   % LGKF run

    if nargin < 1 || isempty(filter_type)
        filter_type = 'ekf';
    end

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', 200);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    noise.std_v = 0.03;
    noise.std_omega = deg2rad(0.2);
    Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
    az_std = deg2rad(5);

    true_state = generate_ground_truth(config, noise);
    [mu_hist, P_hist, true_states] = run_filter_and_collect_histories(true_state, Q, config, az_std, filter_type);

    diagnose_P_growth(mu_hist, P_hist, true_states, config.dt);
end

function [mu_hist, P_hist, true_states] = run_filter_and_collect_histories(true_state, Q, config, az_std, filter_type)
    true_states = true_state(1:5, :)';

    if strcmpi(filter_type, 'lgkf')
        [g, P] = init_lgkf(true_state(:, 1));
        mu_hist = zeros(config.N, 5);
        P_hist = zeros(5, 5, config.N);

        lg_vec = lie2vec_radar(g);
        mu_hist(1, :) = lg_vec.';
        P_hist(:, :, 1) = P;

        for k = 2:config.N
            [g, P] = predict_lgkf(g, P, Q, config.dt);
            if mod(k, config.meas_steps) == 0
                z = get_radar_measurement(true_state(:, k), az_std);
                [g, P] = update_lgkf(g, P, z, az_std);
            end

            lg_vec = lie2vec_radar(g);
            mu_hist(k, :) = lg_vec.';
            P_hist(:, :, k) = P;
        end
    else
        [x, P] = init_ekf(true_state(:, 1));
        mu_hist = zeros(config.N, 5);
        P_hist = zeros(5, 5, config.N);

        mu_hist(1, :) = x.';
        P_hist(:, :, 1) = P;

        for k = 2:config.N
            [x, P] = predict_ekf(x, P, Q, config.dt);
            if mod(k, config.meas_steps) == 0
                z = get_radar_measurement(true_state(:, k), az_std);
                [x, P] = update_ekf(x, P, z, az_std);
            end

            mu_hist(k, :) = x.';
            P_hist(:, :, k) = P;
        end
    end
end
