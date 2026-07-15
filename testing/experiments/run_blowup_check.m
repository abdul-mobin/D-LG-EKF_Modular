function run_blowup_check(filter_type, n_trials)
% Run multiple trajectories and inspect them with check_for_blowup.
%
% Usage:
%   run('testing/experiments/run_blowup_check')                 % defaults to EKF, 1 trial
%   run('testing/experiments/run_blowup_check', 'lgkf', 100)    % runs 100 LGKF trials

if nargin < 1 || isempty(filter_type)
    filter_type = 'ekf';
end
if nargin < 2 || isempty(n_trials)
    n_trials = 1;
end

repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repoRoot);
addpath(genpath(repoRoot));

config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', 20);
config.N = round(config.t_end / config.dt) + 1;
config.meas_steps = round(config.dt_radar / config.dt);

noise.std_v = 0.03;
noise.std_omega = deg2rad(0.2);
Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
az_std = deg2rad(5);

rng(40);
nan_count = 0;
inf_count = 0;
negative_eig_count = 0;
tiny_theta_count = 0;

for trial = 1:n_trials
    true_state = generate_ground_truth(config, noise);
    [mu_hist, P_hist, theta_hist] = run_filter_and_collect_history(true_state, Q, config, az_std, filter_type);

    nan_count = nan_count + any(isnan(mu_hist(:))) + any(isnan(P_hist(:)));
    inf_count = inf_count + any(isinf(mu_hist(:))) + any(isinf(P_hist(:)));

    min_eig = min(eig(squeeze(P_hist(:, :, end))));
    if min_eig < 0
        negative_eig_count = negative_eig_count + 1;
    end
    tiny_theta_count = tiny_theta_count + sum(abs(theta_hist) < 1e-3);

    if n_trials == 1
        check_for_blowup(mu_hist, P_hist, theta_hist);

        figure('Color', 'w', 'Position', [200, 200, 900, 600]);
        set(gca, 'Color', 'w');
        plot(true_state(1, :), true_state(2, :), 'k-', 'LineWidth', 2, 'DisplayName', 'True trajectory');
        hold on;
        plot(mu_hist(:, 1), mu_hist(:, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated trajectory');
        axis equal;
        grid on;
        xlabel('x [m]');
        ylabel('y [m]');
        title(sprintf('%s trajectory comparison', upper(filter_type)));
        legend('Location', 'best');
    end

    theta_err_raw = mu_hist(:,3) - true_state(3,:).';
    theta_err_wrapped = wrap_angle(theta_err_raw);

    fprintf('Raw theta error range: [%.2f, %.2f] rad\n', min(theta_err_raw), max(theta_err_raw));
    fprintf('Wrapped theta error range: [%.2f, %.2f] rad\n', min(theta_err_wrapped), max(theta_err_wrapped));
    fprintf('Max |raw - wrapped| difference: %.2f (multiples of 2pi: %.1f)\n', ...
        max(abs(theta_err_raw - theta_err_wrapped)), max(abs(theta_err_raw - theta_err_wrapped))/(2*pi));
end

if n_trials > 1
    fprintf('Batch summary over %d trials (%s)\n', n_trials, upper(filter_type));
    fprintf('Trials with NaN in state/covariance: %d\n', nan_count);
    fprintf('Trials with Inf in state/covariance: %d\n', inf_count);
    fprintf('Trials with negative final covariance eigenvalue: %d\n', negative_eig_count);
    fprintf('Total steps with |theta| < 1e-3: %d\n', tiny_theta_count);
end
end

function theta_wrapped = wrap_angle(theta)
    theta_wrapped = mod(theta + pi, 2*pi) - pi;
end

function [mu_hist, P_hist, theta_hist] = run_filter_and_collect_history(true_state, Q, config, az_std, filter_type)
if strcmpi(filter_type, 'lgkf')
    [g, P] = init_lgkf(true_state(:, 1));
    mu_hist = zeros(config.N, 5);
    P_hist = zeros(5, 5, config.N);
    theta_hist = zeros(config.N, 1);

    lg_vec = lie2vec_radar(g);
    mu_hist(1, :) = lg_vec.';
    P_hist(:, :, 1) = P;
    theta_hist(1) = lg_vec(3);

    for k = 2:config.N
        [g, P] = predict_lgkf(g, P, Q, config.dt);
        if mod(k, config.meas_steps) == 0
            z = get_radar_measurement(true_state(:, k), az_std);
            [g, P] = update_lgkf(g, P, z, az_std);
        end

        lg_vec = lie2vec_radar(g);
        mu_hist(k, :) = lg_vec.';
        P_hist(:, :, k) = P;
        theta_hist(k) = lg_vec(3);
    end
else
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
end
