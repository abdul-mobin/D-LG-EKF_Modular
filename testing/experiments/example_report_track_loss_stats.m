function example_report_track_loss_stats(n_trials)
    % Example: generate several trajectories, run EKF and LGKF, and summarize
    % track-loss statistics with report_track_loss_stats.
    %
    % Run from MATLAB:
    %   run('testing/experiments/example_report_track_loss_stats')

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    if nargin < 1 || isempty(n_trials)
        n_trials = 10;
    end

    config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', 20);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    noise.std_v = 0.03;
    noise.std_omega = deg2rad(0.2);
    Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
    az_std = deg2rad(5);

    all_est_ekf = cell(n_trials, 1);
    all_est_lgkf = cell(n_trials, 1);
    all_true = cell(n_trials, 1);
    all_theta_ekf = cell(n_trials, 1);
    all_theta_lgkf = cell(n_trials, 1);
    all_theta_true = cell(n_trials, 1);

    for i = 1:n_trials
        true_state = generate_ground_truth(config, noise);
        [ekf_xy, lgkf_xy, theta_ekf, theta_lgkf] = run_filters_for_track_loss(true_state, Q, config, az_std);

        all_est_ekf{i} = ekf_xy;
        all_est_lgkf{i} = lgkf_xy;
        all_true{i} = true_state(1:2, :).';
        all_theta_ekf{i} = theta_ekf;
        all_theta_lgkf{i} = theta_lgkf;
        all_theta_true{i} = true_state(3, :).';
    end

    report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, config.dt, 1000, 20);

    ekf_heading_rmse = compute_heading_rmse(all_theta_ekf, all_theta_true);
    lgkf_heading_rmse = compute_heading_rmse(all_theta_lgkf, all_theta_true);
    fprintf('EKF heading RMSE: %.4f rad (%.2f deg)\n', ekf_heading_rmse, rad2deg(ekf_heading_rmse));
    fprintf('LGKF heading RMSE: %.4f rad (%.2f deg)\n', lgkf_heading_rmse, rad2deg(lgkf_heading_rmse));
end
