function [x, P] = update_ekf(x, P, z, az_std)
    R = diag([max(az_std^2, 1e-8), 10.0^2, 1.0^2]);

    hx = measurement_model_ekf(x);
    H = numerical_measurement_jacobian_ekf(x);

    S = H * P * H' + R;
    K = P * H' / S;
    innov = z - hx;
    innov(1) = atan2(sin(innov(1)), cos(innov(1)));
    x = x + K * innov;
    P = (eye(5) - K * H) * P;
end

function hx = measurement_model_ekf(x)
    alpha = atan2(x(2), x(1));
    rho = sqrt(x(1)^2 + x(2)^2);
    hx = [alpha; rho; x(4) * cos(alpha - x(3))];
end

function H = numerical_measurement_jacobian_ekf(x)
    eps = 1e-6;
    h0 = measurement_model_ekf(x);
    H = zeros(3, 5);

    for i = 1:5
        xp = x;
        xp(i) = xp(i) + eps;
        hp = measurement_model_ekf(xp);

        delta = atan2(sin(hp(1) - h0(1)), cos(hp(1) - h0(1)));
        H(1, i) = delta / eps;
        H(2:3, i) = (hp(2:3) - h0(2:3)) / eps;
    end
end