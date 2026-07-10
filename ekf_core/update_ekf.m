function [x, P] = update_ekf(x, P, z, az_std)
    R = diag([max(az_std^2, 1e-8), 10.0^2, 1.0^2]);
    
    hx = [atan2(x(2), x(1)); sqrt(x(1)^2 + x(2)^2); x(4)*cos(atan2(x(2), x(1))-x(3))];
    H = zeros(3,5);
    rho_sq = max(x(1)^2 + x(2)^2, 1e-8);
    rho = sqrt(rho_sq);
    alpha = atan2(x(2), x(1));
    delta = alpha - x(3);

    H(1,1) = -x(2) / rho_sq;
    H(1,2) = x(1) / rho_sq;

    H(2,1) = x(1) / rho;
    H(2,2) = x(2) / rho;

    H(3,1) = x(4) * sin(delta) * x(2) / rho_sq;
    H(3,2) = -x(4) * sin(delta) * x(1) / rho_sq;
    H(3,3) = x(4) * sin(delta);
    H(3,4) = cos(delta);

    S = H * P * H' + R;
    K = P * H' / S;
    innov = z - hx; innov(1) = atan2(sin(innov(1)), cos(innov(1)));
    x = x + K * innov;
    P = (eye(5) - K * H) * P;
end