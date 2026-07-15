function [ekf_xy, lgkf_xy, ekf_heading, lgkf_heading] = run_filters_for_track_loss(true_state, Q, config, az_std)
    % Run EKF and LGKF over one trajectory and return position and heading histories.

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
        [x, P_ekf] = predict_ekf(x, P_ekf, Q, config.dt);
        [g, P_lgkf] = predict_lgkf(g, P_lgkf, Q, config.dt);

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