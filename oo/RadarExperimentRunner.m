classdef RadarExperimentRunner
    % RadarExperimentRunner centralizes repository setup and experiment orchestration.
    %
    % The numerical filters remain functional; this class only manages paths,
    % defaults, and experiment loops.

    properties (SetAccess = private)
        RepoRoot
    end

    methods
        function obj = RadarExperimentRunner(repoRoot)
            if nargin < 1 || isempty(repoRoot)
                repoRoot = RadarExperimentRunner.findRepoRoot();
            end

            obj.RepoRoot = repoRoot;
            obj.setupPaths();
        end

        function config = buildComparisonConfig(~, numRuns, azimuthStdsDeg, dt, dtRadar, tEnd)
            if nargin < 2 || isempty(numRuns)
                numRuns = 1000;
            end
            if nargin < 3 || isempty(azimuthStdsDeg)
                azimuthStdsDeg = 0:2:10;
            end
            if nargin < 4 || isempty(dt)
                dt = 0.01;
            end
            if nargin < 5 || isempty(dtRadar)
                dtRadar = 0.1;
            end
            if nargin < 6 || isempty(tEnd)
                tEnd = 20;
            end

            config = struct( ...
                'num_runs', numRuns, ...
                'azimuth_stds', azimuthStdsDeg, ...
                'dt', dt, ...
                'dt_radar', dtRadar, ...
                't_end', tEnd);
            config.N = round(config.t_end / config.dt) + 1;
            config.meas_steps = round(config.dt_radar / config.dt);
        end

        function noise = buildProcessNoise(obj)
            noise.std_v = 0.03;
            noise.std_omega = deg2rad(0.2);

            tunedFile = fullfile(obj.RepoRoot, 'tuned_Q_matrices_2param.mat');
            if exist(tunedFile, 'file')
                loaded = load(tunedFile);
                noise.Q_sys_ekf = loaded.Q_ekf_tuned;
                noise.Q_sys_lgkf = loaded.Q_lgkf_tuned;
            else
                noise.Q_sys_ekf = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
                noise.Q_sys_lgkf = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
            end
        end

        function [results_ekf, results_lgkf] = runComparison(obj, config)
            if nargin < 2 || isempty(config)
                config = obj.buildComparisonConfig();
            end

            noise = obj.buildProcessNoise();

            results_ekf = zeros(numel(config.azimuth_stds), 1);
            results_lgkf = zeros(numel(config.azimuth_stds), 1);

            rng(42);
            for az_idx = 1:numel(config.azimuth_stds)
                current_az_std = deg2rad(config.azimuth_stds(az_idx));
                sum_err = [0, 0];

                for run = 1:config.num_runs
                    true_state = generate_ground_truth(config, noise);
                    z_seq = obj.buildMeasurementSequence(true_state, config, current_az_std);
                    [~, ~, ~, ekf_xy] = obj.propagateEkf(true_state, noise.Q_sys_ekf, config, z_seq, current_az_std);
                    [~, ~, ~, lgkf_xy] = obj.propagateLgkf(true_state, noise.Q_sys_lgkf, config, z_seq, current_az_std);

                    true_xy = true_state(1:2, :).';
                    err_run_ekf = sum((true_xy(2:end, :) - ekf_xy(2:end, :)).^2, 2);
                    err_run_lgkf = sum((true_xy(2:end, :) - lgkf_xy(2:end, :)).^2, 2);
                    err_run = [sum(err_run_ekf), sum(err_run_lgkf)];

                    sum_err = sum_err + (err_run ./ config.N);
                end

                results_ekf(az_idx) = sqrt(sum_err(1) / config.num_runs);
                results_lgkf(az_idx) = sqrt(sum_err(2) / config.num_runs);
            end

            plot_results(config.azimuth_stds, results_ekf, results_lgkf);
        end

        function stats = runHeadingRmseTrials(obj, nTrials, tEnd, dt, azStdDeg)
            if nargin < 2 || isempty(nTrials)
                nTrials = 100;
            end
            if nargin < 3 || isempty(tEnd)
                tEnd = 20;
            end
            if nargin < 4 || isempty(dt)
                dt = 0.01;
            end
            if nargin < 5 || isempty(azStdDeg)
                azStdDeg = 5;
            end

            config = obj.buildComparisonConfig(nTrials, [], dt, 0.1, tEnd);
            noise = obj.buildProcessNoise();
            az_std = deg2rad(azStdDeg);

            trajectories = obj.generateTrajectories(config, noise, nTrials, 42);
            [~, ~, all_theta_ekf, all_theta_lgkf] = obj.runTrajectoryBatch(trajectories, noise, config, az_std);
            all_theta_true = cellfun(@(state) state(3, :).', trajectories, 'UniformOutput', false);
            rmse_ekf_each = cellfun(@(theta_hat, theta_true) compute_heading_rmse(theta_hat, theta_true), all_theta_ekf, all_theta_true);
            rmse_lgkf_each = cellfun(@(theta_hat, theta_true) compute_heading_rmse(theta_hat, theta_true), all_theta_lgkf, all_theta_true);

            stats = struct();
            stats.ekfAggregate = compute_heading_rmse(all_theta_ekf, all_theta_true);
            stats.lgkfAggregate = compute_heading_rmse(all_theta_lgkf, all_theta_true);
            stats.ekfTrialMean = mean(rmse_ekf_each);
            stats.ekfTrialStd = std(rmse_ekf_each);
            stats.lgkfTrialMean = mean(rmse_lgkf_each);
            stats.lgkfTrialStd = std(rmse_lgkf_each);
            stats.nTrials = nTrials;
            stats.tEnd = tEnd;
            stats.dt = dt;
            stats.azStdDeg = azStdDeg;

            fprintf('Heading RMSE over %d trials (t_end = %.1f s, az std = %.1f deg)\n', nTrials, tEnd, azStdDeg);
            fprintf('EKF  aggregate RMSE: %.4f rad (%.2f deg)\n', stats.ekfAggregate, rad2deg(stats.ekfAggregate));
            fprintf('LGKF aggregate RMSE: %.4f rad (%.2f deg)\n', stats.lgkfAggregate, rad2deg(stats.lgkfAggregate));
            fprintf('EKF  trial RMSE mean +- std: %.4f +- %.4f rad\n', stats.ekfTrialMean, stats.ekfTrialStd);
            fprintf('LGKF trial RMSE mean +- std: %.4f +- %.4f rad\n', stats.lgkfTrialMean, stats.lgkfTrialStd);
        end

        function stats = runTrackLossSummary(obj, nTrials, tEnd, dt, azStdDeg, threshold, persistSteps)
            if nargin < 2 || isempty(nTrials)
                nTrials = 10;
            end
            if nargin < 3 || isempty(tEnd)
                tEnd = 20;
            end
            if nargin < 4 || isempty(dt)
                dt = 0.01;
            end
            if nargin < 5 || isempty(azStdDeg)
                azStdDeg = 5;
            end
            if nargin < 6 || isempty(threshold)
                threshold = 1000;
            end
            if nargin < 7 || isempty(persistSteps)
                persistSteps = 20;
            end

            config = obj.buildComparisonConfig(nTrials, [], dt, 0.1, tEnd);
            noise = obj.buildProcessNoise();
            az_std = deg2rad(azStdDeg);

            trajectories = obj.generateTrajectories(config, noise, nTrials, 42);
            [all_est_ekf, all_est_lgkf, all_theta_ekf, all_theta_lgkf] = obj.runTrajectoryBatch(trajectories, noise, config, az_std);
            all_true = cellfun(@(state) state(1:2, :).', trajectories, 'UniformOutput', false);
            all_theta_true = cellfun(@(state) state(3, :).', trajectories, 'UniformOutput', false);

            report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, config.dt, threshold, persistSteps);

            stats = struct();
            stats.ekfHeadingRmse = compute_heading_rmse(all_theta_ekf, all_theta_true);
            stats.lgkfHeadingRmse = compute_heading_rmse(all_theta_lgkf, all_theta_true);
            fprintf('EKF heading RMSE: %.4f rad (%.2f deg)\n', stats.ekfHeadingRmse, rad2deg(stats.ekfHeadingRmse));
            fprintf('LGKF heading RMSE: %.4f rad (%.2f deg)\n', stats.lgkfHeadingRmse, rad2deg(stats.lgkfHeadingRmse));
        end

        function [steps, ekfRmseByStep, lgkfRmseByStep] = plotHeadingRmseVsSteps(obj, nTrials, tEnd, dt, azStdDeg)
            if nargin < 2 || isempty(nTrials)
                nTrials = 50;
            end
            if nargin < 3 || isempty(tEnd)
                tEnd = 20;
            end
            if nargin < 4 || isempty(dt)
                dt = 0.01;
            end
            if nargin < 5 || isempty(azStdDeg)
                azStdDeg = 5;
            end

            config = obj.buildComparisonConfig(nTrials, [], dt, 0.1, tEnd);
            noise = obj.buildProcessNoise();
            az_std = deg2rad(azStdDeg);

            trajectories = obj.generateTrajectories(config, noise, nTrials, 42);
            [~, ~, all_theta_ekf, all_theta_lgkf] = obj.runTrajectoryBatch(trajectories, noise, config, az_std);
            all_theta_true = cellfun(@(state) state(3, :).', trajectories, 'UniformOutput', false);
            [ekfRmseByStep, lgkfRmseByStep] = obj.computeStepwiseHeadingRmse(all_theta_ekf, all_theta_lgkf, all_theta_true);
            steps = 1:config.N;

            figure('Color', 'w', 'Position', [200, 200, 1000, 600]);
            plot(steps, ekfRmseByStep, 'r-', 'LineWidth', 2, 'DisplayName', 'EKF');
            hold on;
            plot(steps, lgkfRmseByStep, 'b-.', 'LineWidth', 2, 'DisplayName', 'Invariant EKF / LGKF');
            grid on;
            xlabel('Step index');
            ylabel('Heading RMSE [rad]');
            title(sprintf('Heading RMSE vs Step over %d trials (t_end = %.1f s, dt = %.3f s)', nTrials, tEnd, dt));
            legend('Location', 'best');

            fprintf('Heading RMSE at final step (%d): EKF = %.4f rad, invariant EKF/LGKF = %.4f rad\n', ...
                config.N, ekfRmseByStep(end), lgkfRmseByStep(end));
            fprintf('Total steps = %d (%.1f s / %.3f s)\n', config.N, tEnd, dt);
        end

        function [ekf_xy, lgkf_xy, ekf_heading, lgkf_heading] = runTrajectory(obj, true_state, noise, config, az_std)
            z_seq = obj.buildMeasurementSequence(true_state, config, az_std);
            [~, ~, ekf_heading, ekf_xy] = obj.propagateEkf(true_state, noise.Q_sys_ekf, config, z_seq, az_std);
            [~, ~, lgkf_heading, lgkf_xy] = obj.propagateLgkf(true_state, noise.Q_sys_lgkf, config, z_seq, az_std);
        end

        function [ekf_xy, lgkf_xy] = runPositionTrajectory(obj, true_state, noise, config, az_std)
            [ekf_xy, lgkf_xy, ~, ~] = obj.runTrajectory(true_state, noise, config, az_std);
        end

        function [all_ekf_xy, all_lgkf_xy, all_ekf_heading, all_lgkf_heading] = runTrajectoryBatch(obj, trajectories, noise, config, az_std)
            nTrials = numel(trajectories);
            all_ekf_xy = cell(nTrials, 1);
            all_lgkf_xy = cell(nTrials, 1);
            all_ekf_heading = cell(nTrials, 1);
            all_lgkf_heading = cell(nTrials, 1);

            for i = 1:nTrials
                [all_ekf_xy{i}, all_lgkf_xy{i}, all_ekf_heading{i}, all_lgkf_heading{i}] = ...
                    obj.runTrajectory(trajectories{i}, noise, config, az_std);
            end
        end

        function [states, P_history, theta_history] = runStateHistory(obj, true_state, noise, config, az_std, filterType)
            z_seq = obj.buildMeasurementSequence(true_state, config, az_std);
            if strcmpi(filterType, 'ekf')
                [states, P_history, theta_history, ~] = obj.propagateEkf(true_state, noise.Q_sys_ekf, config, z_seq, az_std);
            else
                [states, P_history, theta_history, ~] = obj.propagateLgkf(true_state, noise.Q_sys_lgkf, config, z_seq, az_std);
            end
        end

        function [all_states, all_P_history, all_theta_history] = runStateHistoryBatch(obj, trajectories, noise, config, az_std, filterType)
            nTrials = numel(trajectories);
            all_states = cell(nTrials, 1);
            all_P_history = cell(nTrials, 1);
            all_theta_history = cell(nTrials, 1);

            for i = 1:nTrials
                [all_states{i}, all_P_history{i}, all_theta_history{i}] = ...
                    obj.runStateHistory(trajectories{i}, noise, config, az_std, filterType);
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

        function trajectories = generateTrajectories(~, config, noise, nTraj, seed)
            if nargin < 5
                seed = [];
            end
            if ~isempty(seed)
                rng(seed);
            end

            trajectories = cell(nTraj, 1);
            for i = 1:nTraj
                trajectories{i} = generate_ground_truth(config, noise);
            end
        end

        function [ekfRmseByStep, lgkfRmseByStep] = computeStepwiseHeadingRmse(obj, all_theta_ekf, all_theta_lgkf, all_theta_true)
            nTrials = numel(all_theta_true);
            nSteps = numel(all_theta_true{1});
            ekf_sq_err_sum = zeros(nSteps, 1);
            lgkf_sq_err_sum = zeros(nSteps, 1);

            for i = 1:nTrials
                ekf_err = obj.wrapToPiLocal(all_theta_ekf{i} - all_theta_true{i});
                lgkf_err = obj.wrapToPiLocal(all_theta_lgkf{i} - all_theta_true{i});
                ekf_sq_err_sum = ekf_sq_err_sum + ekf_err.^2;
                lgkf_sq_err_sum = lgkf_sq_err_sum + lgkf_err.^2;
            end

            ekfRmseByStep = sqrt(ekf_sq_err_sum / nTrials);
            lgkfRmseByStep = sqrt(lgkf_sq_err_sum / nTrials);
        end
    end

    methods (Static)
        function repoRoot = findRepoRoot()
            thisFile = mfilename('fullpath');
            repoRoot = fileparts(fileparts(thisFile));
        end
    end

    methods (Access = private)
        function setupPaths(obj)
            addpath(obj.RepoRoot);
            addpath(genpath(obj.RepoRoot));
        end

        function z_seq = buildMeasurementSequence(~, true_state, config, az_std)
            z_seq = cell(config.N, 1);
            for k = 2:config.N
                if mod(k, config.meas_steps) == 0
                    z_seq{k} = get_radar_measurement(true_state(:, k), az_std);
                end
            end
        end

        function [states, P_history, theta_history, xy] = propagateEkf(~, true_state, Q_ekf, config, z_seq, az_std)
            [x, P] = init_ekf(true_state(:, 1));

            states = zeros(config.N, 5);
            P_history = zeros(5, 5, config.N);
            theta_history = zeros(config.N, 1);
            xy = zeros(config.N, 2);

            states(1, :) = x.';
            P_history(:, :, 1) = P;
            theta_history(1) = x(3);
            xy(1, :) = x(1:2).';

            for k = 2:config.N
                [x, P] = predict_ekf(x, P, Q_ekf, config.dt);
                if ~isempty(z_seq{k})
                    [x, P] = update_ekf(x, P, z_seq{k}, az_std);
                end

                states(k, :) = x.';
                P_history(:, :, k) = P;
                theta_history(k) = x(3);
                xy(k, :) = x(1:2).';
            end
        end

        function [states, P_history, theta_history, xy] = propagateLgkf(~, true_state, Q_lgkf, config, z_seq, az_std)
            [g, P] = init_lgkf(true_state(:, 1));

            states = zeros(config.N, 5);
            P_history = zeros(5, 5, config.N);
            theta_history = zeros(config.N, 1);
            xy = zeros(config.N, 2);

            lg_vec = lie2vec_radar(g);
            states(1, :) = lg_vec.';
            P_history(:, :, 1) = P;
            theta_history(1) = lg_vec(3);
            xy(1, :) = lg_vec(1:2).';

            for k = 2:config.N
                [g, P] = predict_lgkf(g, P, Q_lgkf, config.dt);
                if ~isempty(z_seq{k})
                    [g, P] = update_lgkf(g, P, z_seq{k}, az_std);
                end

                lg_vec = lie2vec_radar(g);
                states(k, :) = lg_vec.';
                P_history(:, :, k) = P;
                theta_history(k) = lg_vec(3);
                xy(k, :) = lg_vec(1:2).';
            end
        end
    end

    methods (Static, Access = private)
        function thetaWrapped = wrapToPiLocal(theta)
            thetaWrapped = mod(theta + pi, 2 * pi) - pi;
        end
    end
end