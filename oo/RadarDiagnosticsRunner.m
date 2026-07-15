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

            rng(42);
            trajectories = cell(10, 1);
            for i = 1:numel(trajectories)
                trajectories{i} = generate_ground_truth(config, noise);
            end

            true_state = trajectories{trajIdx};
            [ekf_xy, lgkf_xy] = obj.runPositionTrajectory(true_state, noise, config, az_std);

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
            nanCount = 0;
            infCount = 0;
            negativeEigCount = 0;
            tinyThetaCount = 0;

            for trial = 1:nTrials
                true_state = generate_ground_truth(config, noise);
                [mu_hist, P_hist, theta_hist] = obj.runStateHistory(true_state, noise, config, az_std, filterType);

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
                thetaErrWrapped = obj.wrapAngle(thetaErrRaw);

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
            [mu_hist, P_hist, true_states] = obj.runStateHistory(true_state, noise, config, az_std, filterType);
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

            rng(42);
            neesAll = [];

            for trial = 1:nTrials
                true_state = generate_ground_truth(config, noise);
                [states, P_history] = obj.runStateHistory(true_state, noise, config, az_std, filterType);
                true_states = true_state(1:5, :)';
                neesAll = [neesAll; compute_nees(states, true_states, P_history)]; %#ok<AGROW>
            end

            fprintf('%s NEES summary over %d trajectories\n', upper(filterType), nTrials);
            check_consistency(neesAll, 5);
            inspect_nees_distribution(neesAll);
        end

        function sweepTrackLossByTEnd(obj)
            t_end_values = [20 50 100 200];
            nTrials = 500;
            threshold = 1000;
            fullLossGroupSteps = 20;

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

                rng(42);
                for trial = 1:nTrials
                    true_state = generate_ground_truth(config, noise);
                    [ekf_xy, lgkf_xy] = obj.runPositionTrajectory(true_state, noise, config, az_std);

                    all_est_ekf{trial} = ekf_xy;
                    all_est_lgkf{trial} = lgkf_xy;
                    all_true{trial} = true_state(1:2, :).';
                end

                fprintf('\n t_end = %.1f s \n', tEnd);
                report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, config.dt, threshold, fullLossGroupSteps);
            end
        end

        function tuneQTwoParams(obj)
            noise = obj.BaseRunner.buildProcessNoise();
            sigma_v = noise.std_v;
            sigma_omega = noise.std_omega;
            Q_true = diag([0, 0, 0, sigma_v^2, sigma_omega^2]);

            config = obj.BaseRunner.buildComparisonConfig(20, [], 0.01, 0.1, 200);
            az_std = deg2rad(5);

            rng(42);
            nTraj = 20;
            trajectories = cell(nTraj, 1);
            for i = 1:nTraj
                trajectories{i} = generate_ground_truth(config, noise);
            end

            opts = optimset('Display', 'iter', 'TolX', 1e-3, 'TolFun', 1e-3);
            x0 = [0, 0];

            objective_ekf = @(x) obj.rmseForTwoScales(x, Q_true, config, az_std, trajectories, 'ekf');
            x_ekf = fminsearch(objective_ekf, x0, opts);
            Q_ekf_tuned = obj.scaleQ(Q_true, x_ekf);

            objective_lgkf = @(x) obj.rmseForTwoScales(x, Q_true, config, az_std, trajectories, 'lgkf');
            x_lgkf = fminsearch(objective_lgkf, x0, opts);
            Q_lgkf_tuned = obj.scaleQ(Q_true, x_lgkf);

            fprintf('EKF  : scale_v = %.4f, scale_omega = %.4f\n', exp(x_ekf(1)), exp(x_ekf(2)));
            fprintf('LGKF : scale_v = %.4f, scale_omega = %.4f\n', exp(x_lgkf(1)), exp(x_lgkf(2)));

            save(fullfile(obj.BaseRunner.RepoRoot, 'tuned_Q_matrices_2param.mat'), 'Q_ekf_tuned', 'Q_lgkf_tuned');
        end
    end

    methods (Access = private)
        function [ekf_xy, lgkf_xy] = runPositionTrajectory(obj, true_state, noise, config, az_std)
            [x, P_ekf] = init_ekf(true_state(:, 1));
            [g, P_lgkf] = init_lgkf(true_state(:, 1));

            ekf_xy = zeros(config.N, 2);
            lgkf_xy = zeros(config.N, 2);

            ekf_xy(1, :) = x(1:2).';
            lg_vec = lie2vec_radar(g);
            lgkf_xy(1, :) = lg_vec(1:2).';

            for k = 2:config.N
                [x, P_ekf] = predict_ekf(x, P_ekf, noise.Q_sys_ekf, config.dt);
                [g, P_lgkf] = predict_lgkf(g, P_lgkf, noise.Q_sys_lgkf, config.dt);

                if mod(k, config.meas_steps) == 0
                    z = get_radar_measurement(true_state(:, k), az_std);
                    [x, P_ekf] = update_ekf(x, P_ekf, z, az_std);
                    [g, P_lgkf] = update_lgkf(g, P_lgkf, z, az_std);
                end

                ekf_xy(k, :) = x(1:2).';
                lg_vec = lie2vec_radar(g);
                lgkf_xy(k, :) = lg_vec(1:2).';
            end
        end

        function [states, P_history, theta_history] = runStateHistory(obj, true_state, noise, config, az_std, filterType)
            if strcmpi(filterType, 'ekf')
                [x, P] = init_ekf(true_state(:, 1));
            else
                [g, P] = init_lgkf(true_state(:, 1));
            end

            states = zeros(config.N, 5);
            P_history = zeros(5, 5, config.N);
            theta_history = zeros(config.N, 1);

            if strcmpi(filterType, 'ekf')
                states(1, :) = x.';
                theta_history(1) = x(3);
            else
                lg_vec = lie2vec_radar(g);
                states(1, :) = lg_vec.';
                theta_history(1) = lg_vec(3);
            end
            P_history(:, :, 1) = P;

            for k = 2:config.N
                if strcmpi(filterType, 'ekf')
                    [x, P] = predict_ekf(x, P, noise.Q_sys_ekf, config.dt);
                    if mod(k, config.meas_steps) == 0
                        z = get_radar_measurement(true_state(:, k), az_std);
                        [x, P] = update_ekf(x, P, z, az_std);
                    end
                    states(k, :) = x.';
                    theta_history(k) = x(3);
                else
                    [g, P] = predict_lgkf(g, P, noise.Q_sys_lgkf, config.dt);
                    if mod(k, config.meas_steps) == 0
                        z = get_radar_measurement(true_state(:, k), az_std);
                        [g, P] = update_lgkf(g, P, z, az_std);
                    end
                    lg_vec = lie2vec_radar(g);
                    states(k, :) = lg_vec.';
                    theta_history(k) = lg_vec(3);
                end
                P_history(:, :, k) = P;
            end
        end

        function rmse = rmseForTwoScales(obj, logScales, QBase, config, az_std, trajectories, filterType)
            Q = obj.scaleQ(QBase, logScales);
            nTraj = numel(trajectories);
            errs = zeros(nTraj, 1);

            for i = 1:nTraj
                traj = trajectories{i};
                if strcmp(filterType, 'ekf')
                    est_xy = obj.runEkfTuning(traj, Q, config, az_std);
                else
                    est_xy = obj.runLgkfTuning(traj, Q, config, az_std);
                end

                pos_err = sqrt(sum((est_xy - traj(1:2, :)').^2, 2));
                errs(i) = sqrt(mean(pos_err.^2));
            end

            rmse = mean(errs);
        end

        function Q = scaleQ(~, QBase, logScales)
            Q = QBase;
            Q(4, 4) = QBase(4, 4) * exp(logScales(1));
            Q(5, 5) = QBase(5, 5) * exp(logScales(2));
        end

        function est_xy = runEkfTuning(~, true_state, Q, config, az_std)
            [x, P] = init_ekf(true_state(:, 1));
            est_xy = zeros(config.N, 2);
            est_xy(1, :) = x(1:2).';

            for k = 2:config.N
                [x, P] = predict_ekf(x, P, Q, config.dt);
                if mod(k, config.meas_steps) == 0
                    z = get_radar_measurement(true_state(:, k), az_std);
                    [x, P] = update_ekf(x, P, z, az_std);
                end
                est_xy(k, :) = x(1:2).';
            end
        end

        function est_xy = runLgkfTuning(~, true_state, Q, config, az_std)
            [g, P] = init_lgkf(true_state(:, 1));
            est_xy = zeros(config.N, 2);
            lg_vec = lie2vec_radar(g);
            est_xy(1, :) = lg_vec(1:2).';

            for k = 2:config.N
                [g, P] = predict_lgkf(g, P, Q, config.dt);
                if mod(k, config.meas_steps) == 0
                    z = get_radar_measurement(true_state(:, k), az_std);
                    [g, P] = update_lgkf(g, P, z, az_std);
                end
                lg_vec = lie2vec_radar(g);
                est_xy(k, :) = lg_vec(1:2).';
            end
        end

        function wrapped = wrapAngle(~, theta)
            wrapped = mod(theta + pi, 2 * pi) - pi;
        end
    end
end
