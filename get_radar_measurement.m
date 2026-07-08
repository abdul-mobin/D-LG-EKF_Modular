
function z = get_radar_measurement(true_vec, az_std)
    tx = true_vec(1); ty = true_vec(2);
    R = diag([max(az_std^2, 1e-8), 10.0^2, 1.0^2]);
    z = [atan2(ty, tx); sqrt(tx^2 + ty^2); true_vec(4)*cos(atan2(ty, tx)-true_vec(3))] ...
        + sqrt(R)*randn(3,1);
end