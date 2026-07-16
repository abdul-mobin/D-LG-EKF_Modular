classdef RadarDiagnosticsRunner
    % RadarDiagnosticsRunner owns the remaining experiment and diagnostics workflows.
    %
    % The low-level EKF/LGKF kernels remain functional; this class only wraps
    % orchestration, batch runs, and tuning experiments.

    properties (SetAccess = private)
        BaseRunner
    end

    methods
        function obj = RadarDiagnosticsRunner(repoRoot)
            if nargin < 1 || isempty(repoRoot)
                obj.BaseRunner = RadarExperimentRunner();
            else
                obj.BaseRunner = RadarExperimentRunner(repoRoot);
            end
        end

        function runSingleTrajectoryComparison(obj, trajIdx, azStdDeg)
            if nargin < 2 || isempty(trajIdx)
                trajIdx = 1;
            end
            if nargin < 3 || isempty(azStdDeg)
                azStdDeg = 5;
            end

            config = obj.BaseRunner.buildComparisonConfig(10, [], 0.01, 0.1, 200);
            noise = obj.BaseRunner.buildProcessNoise();
            az_std = deg2rad(azStdDeg);

            trajectories = obj.BaseRunner.generateTrajectories(config, noise, 10, 42);

            true_state = trajectories{trajIdx};
            [ekf_xy, lgkf_xy] = obj.BaseRunner.runPositionTrajectory(true_state, noise, config, az_std);

            figure('Color', 'w', 'Position', [200, 200, 900, 600]);
            plot(true_state(1, :), true_state(2, :), 'k-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
            hold on;
            plot(ekf_xy(:, 1), ekf_xy(:, 2), 'r--', 'LineWidth', 2, 'DisplayName', 'EKF estimate');
            plot(lgkf_xy(:, 1), lgkf_xy(:, 2), 'b-.', 'LineWidth', 2, 'DisplayName', 'LGKF estimate');
            grid on;
            axis equal;
            xlabel('x [m]');
            ylabel('y [m]');
            title(sprintf('Trajectory %d and Filter Estimates (azimuth std = %.1f deg)', trajIdx, azStdDeg));
            legend('Location', 'best');
        end

        function runBlowupCheck(obj, filterType, nTrials)
            if nargin < 2 || isempty(filterType)
                filterType = 'ekf';
            end
            if nargin < 3 || isempty(nTrials)
                nTrials = 1;
            end

            config = obj.BaseRunner.buildComparisonConfig(1, [], 0.01, 0.1, 20);
            noise = obj.BaseRunner.buildProcessNoise();
            az_std = deg2rad(5);

            rng(40);
            trajectories = obj.BaseRunner.generateTrajectories(config, noise, nTrials, 40);
            nanCount = 0;
            infCount = 0;
            negativeEigCount = 0;
            tinyThetaCount = 0;

            [all_states, all_P_history, all_theta_history] = ...
                obj.BaseRunner.runStateHistoryBatch(trajectories, noise, config, az_std, filterType);

            for trial = 1:nTrials
                true_state = trajectories{trial};
                mu_hist = all_states{trial};
                P_hist = all_P_history{trial};
                theta_hist = all_theta_history{trial};

                nanCount = nanCount + any(isnan(mu_hist(:))) + any(isnan(P_hist(:)));
                infCount = infCount + any(isinf(mu_hist(:))) + any(isinf(P_hist(:)));

                minEig = min(eig(squeeze(P_hist(:, :, end))));
                if minEig < 0
                    negativeEigCount = negativeEigCount + 1;
                end
                tinyThetaCount = tinyThetaCount + sum(abs(theta_hist) < 1e-3);

                if nTrials == 1
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
                    title(sprintf('%s trajectory comparison', upper(filterType)));
                    legend('Location', 'best');
                end

                thetaErrRaw = mu_hist(:, 3) - true_state(3, :).';
                thetaErrWrapped = obj.BaseRunner.wrapAngle(thetaErrRaw);

                fprintf('Raw theta error range: [%.2f, %.2f] rad\n', min(thetaErrRaw), max(thetaErrRaw));
                fprintf('Wrapped theta error range: [%.2f, %.2f] rad\n', min(thetaErrWrapped), max(thetaErrWrapped));
                fprintf('Max |raw - wrapped| difference: %.2f (multiples of 2pi: %.1f)\n', ...
                    max(abs(thetaErrRaw - thetaErrWrapped)), max(abs(thetaErrRaw - thetaErrWrapped)) / (2 * pi));
            end

            if nTrials > 1
                fprintf('Batch summary over %d trials (%s)\n', nTrials, upper(filterType));
                fprintf('Trials with NaN in state/covariance: %d\n', nanCount);
                fprintf('Trials with Inf in state/covariance: %d\n', infCount);
                fprintf('Trials with negative final covariance eigenvalue: %d\n', negativeEigCount);
                fprintf('Total steps with |theta| < 1e-3: %d\n', tinyThetaCount);
            end
        end

        function runDiagnosePGrowth(obj, filterType)
            if nargin < 2 || isempty(filterType)
                filterType = 'ekf';
            end

            config = obj.BaseRunner.buildComparisonConfig(1, [], 0.01, 0.1, 200);
            noise = obj.BaseRunner.buildProcessNoise();
            az_std = deg2rad(5);

            true_state = generate_ground_truth(config, noise);
            [mu_hist, P_hist, ~] = obj.BaseRunner.runStateHistory(true_state, noise, config, az_std, filterType);
            true_states = true_state(1:5, :)';
            diagnose_P_growth(mu_hist, P_hist, true_states, config.dt);
        end

        function runNeesAnalysis(obj, nTrials, filterType)
            if nargin < 2 || isempty(nTrials)
                nTrials = 50;
            end
            if nargin < 3 || isempty(filterType)
                filterType = 'lgkf';
            end

            config = obj.BaseRunner.buildComparisonConfig(nTrials, [], 0.01, 0.1, 200);
            noise = obj.BaseRunner.buildProcessNoise();
            az_std = deg2rad(5);

            trajectories = obj.BaseRunner.generateTrajectories(config, noise, nTrials, 42);
            [all_states, all_P_history, ~] = obj.BaseRunner.runStateHistoryBatch(trajectories, noise, config, az_std, filterType);
            nees_cells = cell(nTrials, 1);
            for trial = 1:nTrials
                true_states = trajectories{trial}(1:5, :)';
                nees_cells{trial} = compute_nees(all_states{trial}, true_states, all_P_history{trial});
            end
            neesAll = vertcat(nees_cells{:});

            fprintf('%s NEES summary over %d trajectories\n', upper(filterType), nTrials);
            check_consistency(neesAll, 5);
            inspect_nees_distribution(neesAll);
        end

        function sweepTrackLossByTEnd(obj, t_end_values, nTrials, threshold, fullLossGroupSteps)
            if nargin < 2 || isempty(t_end_values)
                t_end_values = [20 50 100 200];
            end
            if nargin < 3 || isempty(nTrials)
                nTrials = 500;
            end
            if nargin < 4 || isempty(threshold)
                threshold = 1000;
            end
            if nargin < 5 || isempty(fullLossGroupSteps)
                fullLossGroupSteps = 20;
            end

            noise = obj.BaseRunner.buildProcessNoise();
            az_std = deg2rad(5);

            fprintf('Sweeping track-loss statistics over t_end values...\n');
            fprintf('Threshold = %.1f m, full-loss group = %d steps\n', threshold, fullLossGroupSteps);

            for idx = 1:numel(t_end_values)
                tEnd = t_end_values(idx);
                config = obj.BaseRunner.buildComparisonConfig(nTrials, [], 0.01, 0.1, tEnd);

                all_est_ekf = cell(nTrials, 1);
                all_est_lgkf = cell(nTrials, 1);
                all_true = cell(nTrials, 1);

                trajectories = obj.BaseRunner.generateTrajectories(config, noise, nTrials, 42);
                [all_est_ekf, all_est_lgkf, ~, ~] = obj.BaseRunner.runTrajectoryBatch(trajectories, noise, config, az_std);
                all_true = cellfun(@(state) state(1:2, :).', trajectories, 'UniformOutput', false);

                fprintf('\n t_end = %.1f s \n', tEnd);
                report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, config.dt, threshold, fullLossGroupSteps);
            end
        end

        function tuneQTwoParams(obj, nTraj)
            if nargin < 2 || isempty(nTraj)
                nTraj = 20;
            end

            noise = obj.BaseRunner.buildProcessNoise();
            sigma_v = noise.std_v;
            sigma_omega = noise.std_omega;
            Q_true = diag([0, 0, 0, sigma_v^2, sigma_omega^2]);

            config = obj.BaseRunner.buildComparisonConfig(nTraj, [], 0.01, 0.1, 200);
            az_std = deg2rad(5);

            trajectories = obj.BaseRunner.generateTrajectories(config, noise, nTraj, 42);

            opts = optimset('Display', 'iter', 'TolX', 1e-3, 'TolFun', 1e-3);
            x0 = [0, 0];

            objective_ekf = @(x) obj.BaseRunner.rmseForTwoScales(x, Q_true, config, az_std, trajectories, 'ekf');
            x_ekf = fminsearch(objective_ekf, x0, opts);
            Q_ekf_tuned = obj.BaseRunner.scaleQ(Q_true, x_ekf);

            objective_lgkf = @(x) obj.BaseRunner.rmseForTwoScales(x, Q_true, config, az_std, trajectories, 'lgkf');
            x_lgkf = fminsearch(objective_lgkf, x0, opts);
            Q_lgkf_tuned = obj.BaseRunner.scaleQ(Q_true, x_lgkf);

            fprintf('EKF  : scale_v = %.4f, scale_omega = %.4f\n', exp(x_ekf(1)), exp(x_ekf(2)));
            fprintf('LGKF : scale_v = %.4f, scale_omega = %.4f\n', exp(x_lgkf(1)), exp(x_lgkf(2)));

            save(fullfile(obj.BaseRunner.RepoRoot, 'tuned_Q_matrices_2param.mat'), 'Q_ekf_tuned', 'Q_lgkf_tuned');
        end
    end
end
