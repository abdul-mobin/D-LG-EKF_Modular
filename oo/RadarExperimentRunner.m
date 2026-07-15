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
                    [state_ekf, P_ekf] = init_ekf(true_state(:, 1));
                    [g_lgkf, P_lgkf] = init_lgkf(true_state(:, 1));

                    err_run = [0, 0];
                    for k = 2:config.N
                        [state_ekf, P_ekf] = predict_ekf(state_ekf, P_ekf, noise.Q_sys_ekf, config.dt);
                        [g_lgkf, P_lgkf] = predict_lgkf(g_lgkf, P_lgkf, noise.Q_sys_lgkf, config.dt);

                        if mod(k, config.meas_steps) == 0
                            z = get_radar_measurement(true_state(:, k), current_az_std);
                            [state_ekf, P_ekf] = update_ekf(state_ekf, P_ekf, z, current_az_std);
                            [g_lgkf, P_lgkf] = update_lgkf(g_lgkf, P_lgkf, z, current_az_std);
                        end

                        lg_vec = lie2vec_radar(g_lgkf);
                        err_run(1) = err_run(1) + (true_state(1, k) - state_ekf(1))^2 + (true_state(2, k) - state_ekf(2))^2;
                        err_run(2) = err_run(2) + (true_state(1, k) - lg_vec(1))^2 + (true_state(2, k) - lg_vec(2))^2;
                    end

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

            all_theta_ekf = cell(nTrials, 1);
            all_theta_lgkf = cell(nTrials, 1);
            all_theta_true = cell(nTrials, 1);
            rmse_ekf_each = zeros(nTrials, 1);
            rmse_lgkf_each = zeros(nTrials, 1);

            rng(42);
            for i = 1:nTrials
                true_state = generate_ground_truth(config, noise);
                [~, ~, theta_ekf, theta_lgkf] = obj.runTrajectory(true_state, noise, config, az_std);

                theta_true = true_state(3, :).';
                all_theta_ekf{i} = theta_ekf;
                all_theta_lgkf{i} = theta_lgkf;
                all_theta_true{i} = theta_true;

                rmse_ekf_each(i) = compute_heading_rmse(theta_ekf, theta_true);
                rmse_lgkf_each(i) = compute_heading_rmse(theta_lgkf, theta_true);
            end

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

            all_est_ekf = cell(nTrials, 1);
            all_est_lgkf = cell(nTrials, 1);
            all_true = cell(nTrials, 1);
            all_theta_ekf = cell(nTrials, 1);
            all_theta_lgkf = cell(nTrials, 1);
            all_theta_true = cell(nTrials, 1);

            rng(42);
            for i = 1:nTrials
                true_state = generate_ground_truth(config, noise);
                [ekf_xy, lgkf_xy, theta_ekf, theta_lgkf] = obj.runTrajectory(true_state, noise, config, az_std);

                all_est_ekf{i} = ekf_xy;
                all_est_lgkf{i} = lgkf_xy;
                all_true{i} = true_state(1:2, :).';
                all_theta_ekf{i} = theta_ekf;
                all_theta_lgkf{i} = theta_lgkf;
                all_theta_true{i} = true_state(3, :).';
            end

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

            ekf_sq_err_sum = zeros(config.N, 1);
            lgkf_sq_err_sum = zeros(config.N, 1);

            rng(42);
            for trial = 1:nTrials
                true_state = generate_ground_truth(config, noise);
                [~, ~, theta_ekf, theta_lgkf] = obj.runTrajectory(true_state, noise, config, az_std);

                theta_true = true_state(3, :).';
                ekf_err = obj.wrapToPiLocal(theta_ekf - theta_true);
                lgkf_err = obj.wrapToPiLocal(theta_lgkf - theta_true);

                ekf_sq_err_sum = ekf_sq_err_sum + ekf_err.^2;
                lgkf_sq_err_sum = lgkf_sq_err_sum + lgkf_err.^2;
            end

            ekfRmseByStep = sqrt(ekf_sq_err_sum / nTrials);
            lgkfRmseByStep = sqrt(lgkf_sq_err_sum / nTrials);
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
            [x, P_ekf] = init_ekf(true_state(:, 1));
            [g, P_lgkf] = init_lgkf(true_state(:, 1));

            ekf_xy = zeros(config.N, 2);
            lgkf_xy = zeros(config.N, 2);
            ekf_heading = zeros(config.N, 1);
            lgkf_heading = zeros(config.N, 1);

            ekf_xy(1, :) = x(1:2).';
            lg_vec = lie2vec_radar(g);
            lgkf_xy(1, :) = lg_vec(1:2).';
            ekf_heading(1) = x(3);
            lgkf_heading(1) = lg_vec(3);

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
                ekf_heading(k) = x(3);
                lgkf_heading(k) = lg_vec(3);
            end
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
    end

    methods (Static, Access = private)
        function thetaWrapped = wrapToPiLocal(theta)
            thetaWrapped = mod(theta + pi, 2 * pi) - pi;
        end
    end
end