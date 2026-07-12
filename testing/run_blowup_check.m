function run_blowup_check()
    % Run a short example trajectory and inspect it with check_for_blowup.
    %
    % Usage:
    %   run('testing/run_blowup_check')

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

    [mu_hist, P_hist, theta_hist] = run_filter_and_collect_history(true_state, Q, config, az_std);
    check_for_blowup(mu_hist, P_hist, theta_hist);
end

function [mu_hist, P_hist, theta_hist] = run_filter_and_collect_history(true_state, Q, config, az_std)
    [x, P] = init_ekf(true_state(:, 1));

    mu_hist = zeros(config.N, 5);
    P_hist = zeros(5, 5, config.N);
    theta_hist = zeros(config.N, 1);

    mu_hist(1, :) = x.';
    P_hist(:, :, 1) = P;
    theta_hist(1) = x(3);

    for k = 2:config.N
        [x, P] = predict_ekf(x, P, Q, config.dt);
        if mod(k, config.meas_steps) == 0
            z = get_radar_measurement(true_state(:, k), az_std);
            [x, P] = update_ekf(x, P, z, az_std);
        end

        mu_hist(k, :) = x.';
        P_hist(:, :, k) = P;
        theta_hist(k) = x(3);
    end
end
