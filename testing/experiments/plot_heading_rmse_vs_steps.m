function plot_heading_rmse_vs_steps(n_trials, t_end, dt, az_std_deg)
    % Plot heading RMSE versus step index for EKF and invariant EKF.
    %
    % In this repository, the invariant EKF-style filter is implemented by the
    % LGKF Lie-group filter.
    %
    % Usage:
    %   run('testing/experiments/plot_heading_rmse_vs_steps')
    %   plot_heading_rmse_vs_steps(50, 20, 0.01, 5)

    if nargin < 1 || isempty(n_trials)
        n_trials = 100;
    end
    if nargin < 2 || isempty(t_end)
        t_end = 20;
    end
    if nargin < 3 || isempty(dt)
        dt = 0.01;
    end
    if nargin < 4 || isempty(az_std_deg)
        az_std_deg = 5;
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    config = struct('dt', dt, 'dt_radar', 0.1, 't_end', t_end);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    noise.std_v = 0.03;
    noise.std_omega = deg2rad(0.2);
    Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
    az_std = deg2rad(az_std_deg);

    ekf_sq_err_sum = zeros(config.N, 1);
    lgkf_sq_err_sum = zeros(config.N, 1);

    rng(42);
    for trial = 1:n_trials
        true_state = generate_ground_truth(config, noise);
        [~, ~, theta_ekf, theta_lgkf] = run_filters_for_track_loss(true_state, Q, config, az_std);

        theta_true = true_state(3, :).';
        ekf_err = wrap_to_pi(theta_ekf - theta_true);
        lgkf_err = wrap_to_pi(theta_lgkf - theta_true);

        ekf_sq_err_sum = ekf_sq_err_sum + ekf_err.^2;
        lgkf_sq_err_sum = lgkf_sq_err_sum + lgkf_err.^2;
    end

    ekf_rmse_by_step = sqrt(ekf_sq_err_sum / n_trials);
    lgkf_rmse_by_step = sqrt(lgkf_sq_err_sum / n_trials);

    steps = 1:config.N;
    time = (steps - 1) * config.dt;

    figure('Color', 'w', 'Position', [200, 200, 1000, 600]);
    plot(steps, ekf_rmse_by_step, 'r-', 'LineWidth', 2, 'DisplayName', 'EKF');
    hold on;
    plot(steps, lgkf_rmse_by_step, 'b-.', 'LineWidth', 2, 'DisplayName', 'Invariant EKF / LGKF');
    grid on;
    xlabel('Step index');
    ylabel('Heading RMSE [rad]');
    title(sprintf('Heading RMSE vs Step over %d trials (t_end = %.1f s, dt = %.3f s)', n_trials, t_end, dt));
    legend('Location', 'best');

    fprintf('Heading RMSE at final step (%d): EKF = %.4f rad, invariant EKF/LGKF = %.4f rad\n', ...
        config.N, ekf_rmse_by_step(end), lgkf_rmse_by_step(end));
    fprintf('Total steps = %d (%.1f s / %.3f s)\n', config.N, t_end, dt);
    fprintf('Equivalent final time axis: %.1f s\n', time(end));
end

function theta_wrapped = wrap_to_pi(theta)
    theta_wrapped = mod(theta + pi, 2 * pi) - pi;
end
