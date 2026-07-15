function plot_single_trajectory_comparison(traj_idx, az_std_deg)
    % Plot one ground-truth trajectory and the corresponding EKF/LGKF estimates.
    %
    % Inputs:
    %   traj_idx   - index of the trajectory to visualize (1-based)
    %   az_std_deg - azimuth measurement noise in degrees
    %
    % Example:
    %   plot_single_trajectory_comparison(1, 5)

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', 200);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    sigma_v = 0.03;
    sigma_omega = deg2rad(0.2);
    noise.std_v = sigma_v;
    noise.std_omega = sigma_omega;

    az_std = deg2rad(az_std_deg);
    Q = diag([0, 0, 0, sigma_v^2, sigma_omega^2]);

    % rng(42);
    trajectories = cell(10, 1);
    for i = 1:numel(trajectories)
        trajectories{i} = generate_ground_truth(config, noise);
    end

    true_state = trajectories{traj_idx};
        
    [ekf_xy, lgkf_xy] = run_filters_for_plot(true_state, Q, config, az_std);

    figure('Color', 'w', 'Position', [200, 200, 900, 600]);
    plot(true_state(1, :), true_state(2, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
    hold on;
    plot(ekf_xy(:, 1), ekf_xy(:, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'EKF estimate');
    plot(lgkf_xy(:, 1), lgkf_xy(:, 2), 'b-.', 'LineWidth', 2, 'DisplayName', 'LGKF estimate');

    grid on;
    axis equal;
    xlabel('x [m]');
    ylabel('y [m]');
    title(sprintf('Trajectory %d and Filter Estimates (azimuth std = %.1f°)', traj_idx, az_std_deg));
    legend('Location', 'best');
end

function [ekf_xy, lgkf_xy] = run_filters_for_plot(true_state, Q, config, az_std)
    [x, P] = init_ekf(true_state(:, 1));
    ekf_xy = zeros(config.N, 2);
    ekf_xy(1, :) = x(1:2).';

    [g, P_lg] = init_lgkf(true_state(:, 1));
    lgkf_xy = zeros(config.N, 2);
    lg_vec = lie2vec_radar(g);
    lgkf_xy(1, :) = lg_vec(1:2).';

    for k = 2:config.N
        [x, P] = predict_ekf(x, P, Q, config.dt);
        [g, P_lg] = predict_lgkf(g, P_lg, Q, config.dt);

        if mod(k, config.meas_steps) == 0
            z = get_radar_measurement(true_state(:, k), az_std);
            [x, P] = update_ekf(x, P, z, az_std);
            [g, P_lg] = update_lgkf(g, P_lg, z, az_std);
        end

        ekf_xy(k, :) = x(1:2).';
        lg_vec = lie2vec_radar(g);
        lgkf_xy(k, :) = lg_vec(1:2).';
    end
end
