function sweep_track_loss_by_t_end()
    % Sweep track-loss behavior over several values of t_end for EKF and LGKF.
    %
    % Run from MATLAB:
    %   run('testing/experiments/sweep_track_loss_by_t_end')

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    t_end_values = [20 50 100 200];
    n_trials = 500;
    threshold = 1000;
    full_loss_group_steps = 20;

    noise.std_v = 0.03;
    noise.std_omega = deg2rad(0.2);
    Q = diag([0, 0, 0, noise.std_v^2, noise.std_omega^2]);
    az_std = deg2rad(5);

    fprintf('Sweeping track-loss statistics over t_end values...\n');
    fprintf('Threshold = %.1f m, full-loss group = %d steps\n', threshold, full_loss_group_steps);

    for idx = 1:numel(t_end_values)
        t_end = t_end_values(idx);
        config = struct('dt', 0.01, 'dt_radar', 0.1, 't_end', t_end);
        config.N = round(config.t_end / config.dt) + 1;
        config.meas_steps = round(config.dt_radar / config.dt);

        all_est_ekf = cell(n_trials, 1);
        all_est_lgkf = cell(n_trials, 1);
        all_true = cell(n_trials, 1);

        for trial = 1:n_trials
            true_state = generate_ground_truth(config, noise);
            [ekf_xy, lgkf_xy] = run_filters_for_track_loss(true_state, Q, config, az_std);

            all_est_ekf{trial} = ekf_xy;
            all_est_lgkf{trial} = lgkf_xy;
            all_true{trial} = true_state(1:2, :).';
        end

        fprintf('\n t_end = %.1f s \n', t_end);
        report_track_loss_stats(all_est_ekf, all_est_lgkf, all_true, config.dt, threshold, full_loss_group_steps);
    end
end
