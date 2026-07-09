function [x, P] = update_ekf(x, P, z, az_std)
    R = diag([max(az_std^2, 1e-8), 10.0^2, 1.0^2]);
    hx = [atan2(x(2), x(1)); sqrt(x(1)^2 + x(2)^2); x(4)*cos(atan2(x(2), x(1))-x(3))];
    H = zeros(3,5); eps = 1e-5;
    for i=1:5
        xp = x; xp(i) = xp(i)+eps;
        hp = [atan2(xp(2), xp(1)); sqrt(xp(1)^2 + xp(2)^2); xp(4)*cos(atan2(xp(2), xp(1))-xp(3))];
        H(:,i) = (hp - hx) / eps;
    end
    S = H * P * H' + R;
    K = P * H' / S;
    innov = z - hx; innov(1) = atan2(sin(innov(1)), cos(innov(1)));
    x = x + K * innov;
    P = (eye(5) - K * H) * P;
end