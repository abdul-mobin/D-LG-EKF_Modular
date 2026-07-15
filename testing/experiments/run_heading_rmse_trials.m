function run_heading_rmse_trials(n_trials, t_end, az_std_deg)
    % Run multiple trajectories and report heading RMSE for EKF and LGKF.
    %
    % Usage:
    %   run('testing/experiments/run_heading_rmse_trials')
    %   run_heading_rmse_trials(100, 200, 5)

    if nargin < 1 || isempty(n_trials)
        n_trials = 100;
    end
    if nargin < 2 || isempty(t_end)
        t_end = 200;
    end
    if nargin < 3 || isempty(az_std_deg)
        az_std_deg = 5;
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', t_end);
    config.N = round(config.t_end / config.dt) + 1;
    config.meas_steps = round(config.dt_radar / config.dt);

    noise.std_v = 0.03;
    noise.std_omega = deg2rad(0.2);
    Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
    az_std = deg2rad(az_std_deg);

    all_theta_ekf = cell(n_trials, 1);
    all_theta_lgkf = cell(n_trials, 1);
    all_theta_true = cell(n_trials, 1);
    rmse_ekf_each = zeros(n_trials, 1);
    rmse_lgkf_each = zeros(n_trials, 1);

    rng(42);
    for i = 1:n_trials
        true_state = generate_ground_truth(config, noise);
        [~, ~, theta_ekf, theta_lgkf] = run_filters_for_track_loss(true_state, Q, config, az_std);

        theta_true = true_state(3, :).';
        all_theta_ekf{i} = theta_ekf;
        all_theta_lgkf{i} = theta_lgkf;
        all_theta_true{i} = theta_true;

        rmse_ekf_each(i) = compute_heading_rmse(theta_ekf, theta_true);
        rmse_lgkf_each(i) = compute_heading_rmse(theta_lgkf, theta_true);
    end

    rmse_ekf_all = compute_heading_rmse(all_theta_ekf, all_theta_true);
    rmse_lgkf_all = compute_heading_rmse(all_theta_lgkf, all_theta_true);

    fprintf('Heading RMSE over %d trials (t_end = %.1f s, az std = %.1f deg)\n', n_trials, t_end, az_std_deg);
    fprintf('EKF  aggregate RMSE: %.4f rad (%.2f deg)\n', rmse_ekf_all, rad2deg(rmse_ekf_all));
    fprintf('LGKF aggregate RMSE: %.4f rad (%.2f deg)\n', rmse_lgkf_all, rad2deg(rmse_lgkf_all));
    fprintf('EKF  trial RMSE mean +- std: %.4f +- %.4f rad\n', mean(rmse_ekf_each), std(rmse_ekf_each));
    fprintf('LGKF trial RMSE mean +- std: %.4f +- %.4f rad\n', mean(rmse_lgkf_each), std(rmse_lgkf_each));
end
