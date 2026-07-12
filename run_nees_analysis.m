function run_nees_analysis()
    % Example wrapper that computes NEES for one EKF and one LGKF run.
    %
    % This script uses the repository's existing filter routines and a single
    % simulated trajectory. You can adapt it to your Monte-Carlo sweep as needed.

    repoRoot = fileparts(mfilename('fullpath'));
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

    [ekf_states, ekf_P_history] = run_filter_for_nees(true_state, Q, config, az_std, 'ekf');
    [lgkf_states, lgkf_P_history] = run_filter_for_nees(true_state, Q, config, az_std, 'lgkf');

    true_states = true_state(1:5, :)';

    nees_ekf = compute_nees(ekf_states, true_states, ekf_P_history);
    nees_lgkf = compute_nees(lgkf_states, true_states, lgkf_P_history);

    fprintf('EKF NEES summary\n');
    check_consistency(nees_ekf, 5);
    fprintf('\nLGKF NEES summary\n');
    check_consistency(nees_lgkf, 5);
end

function [states, P_history] = run_filter_for_nees(true_state, Q, config, az_std, filter_type)
    if strcmp(filter_type, 'ekf')
        [x, P] = init_ekf(true_state(:, 1));
    else
        [g, P] = init_lgkf(true_state(:, 1));
    end

    states = zeros(config.N, 5);
    P_history = zeros(5, 5, config.N);

    if strcmp(filter_type, 'ekf')
        states(1, :) = x.';
    else
        lg_vec = lie2vec_radar(g);
        states(1, :) = lg_vec.';
    end
    P_history(:, :, 1) = P;

    for k = 2:config.N
        if strcmp(filter_type, 'ekf')
            [x, P] = predict_ekf(x, P, Q, config.dt);
            if mod(k, config.meas_steps) == 0
                z = get_radar_measurement(true_state(:, k), az_std);
                [x, P] = update_ekf(x, P, z, az_std);
            end
            states(k, :) = x.';
        else
            [g, P] = predict_lgkf(g, P, Q, config.dt);
            if mod(k, config.meas_steps) == 0
                z = get_radar_measurement(true_state(:, k), az_std);
                [g, P] = update_lgkf(g, P, z, az_std);
            end
            lg_vec = lie2vec_radar(g);
            states(k, :) = lg_vec.';
        end
        P_history(:, :, k) = P;
    end
end
