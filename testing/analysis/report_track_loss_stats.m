function report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, T, threshold, persist_steps)
    % Report track-loss statistics for EKF and LGKF over a set of trajectories.
    %
    % Inputs:
    %   all_est_ekf   - cell array of EKF position histories (N x 2)
    %   all_est_lgkf  - cell array of LGKF position histories (N x 2)
    %   all_true      - cell array of true position histories (N x 2)
    %   T             - sampling interval [s]
    %   threshold     - position error threshold [m]
    %   persist_steps - required consecutive steps above threshold

    if nargin < 5 || isempty(threshold)
        threshold = 1000;
    end
    if nargin < 6 || isempty(persist_steps)
        persist_steps = 20;
    end

    N = numel(all_true);
    loss_ekf = false(N, 1);
    time_ekf = nan(N, 1);
    loss_lgkf = false(N, 1);
    time_lgkf = nan(N, 1);

    for i = 1:N
        [loss_ekf(i), time_ekf(i)] = detect_track_loss(all_est_ekf{i}, all_true{i}, T, threshold, persist_steps);
        [loss_lgkf(i), time_lgkf(i)] = detect_track_loss(all_est_lgkf{i}, all_true{i}, T, threshold, persist_steps);
    end

    fprintf('EKF track loss rate: %.1f%% (%d/%d trajectories)\n', 100 * mean(loss_ekf), sum(loss_ekf), N);
    fprintf('LGKF track loss rate: %.1f%% (%d/%d trajectories)\n', 100 * mean(loss_lgkf), sum(loss_lgkf), N);
    fprintf('EKF mean time-to-loss (lost trajectories only): %.1f s\n', mean(time_ekf(loss_ekf)));
    fprintf('LGKF mean time-to-loss (lost trajectories only): %.1f s\n', mean(time_lgkf(loss_lgkf)));
end
