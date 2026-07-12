function run_nees_analysis(n_trials)
    % Compute NEES values for many trajectories using the LGKF (and optionally EKF).
    %
    % Usage:
    %   run('testing/run_nees_analysis')          % defaults to 50 trials
    %   run('testing/run_nees_analysis', 200)     % run 200 trajectories

    if nargin < 1 || isempty(n_trials)
        n_trials = 50;
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

    rng(42);
    nees_lgkf_all = [];

    for trial = 1:n_trials
        true_state = generate_ground_truth(config, noise);

        [~, lgkf_P_history] = run_filter_for_nees(true_state, Q, config, az_std, 'lgkf');
        true_states = true_state(1:5, :)';

        % For LGKF, the state is stored in the Lie-group coordinates and must be
        % converted back to the Euclidean state vector before NEES computation.
        [lgkf_states, ~] = run_filter_for_nees(true_state, Q, config, az_std, 'lgkf');
        nees_lgkf = compute_nees(lgkf_states, true_states, lgkf_P_history);
        nees_lgkf_all = [nees_lgkf_all; nees_lgkf];
    end

    fprintf('LGKF NEES summary over %d trajectories\n', n_trials);
    check_consistency(nees_lgkf_all, 5);
    inspect_nees_distribution(nees_lgkf_all);
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
