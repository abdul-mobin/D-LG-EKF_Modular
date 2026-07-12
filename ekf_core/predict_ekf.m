function [x, P] = predict_ekf(x, P, Q, dt)
    f = state_transition_ekf(x, dt);
    F = numerical_state_jacobian_ekf(x, dt);

    x = f;
    P = F * P * F' + Q;
end

function f = state_transition_ekf(x, dt)
    v = x(4);
    w = x(5);
    th = x(3);

    if abs(w) > 1e-5
        f = [x(1) + (v/w) * (sin(th + w*dt) - sin(th));
             x(2) + (v/w) * (cos(th) - cos(th + w*dt));
             th + w*dt;
             v;
             w];
    else
        f = [x(1) + v*cos(th)*dt;
             x(2) + v*sin(th)*dt;
             th + w*dt;
             v;
             w];
    end
end

function F = numerical_state_jacobian_ekf(x, dt)
    eps = 1e-6;
    f0 = state_transition_ekf(x, dt);
    F = zeros(5, 5);

    for i = 1:5
        xp = x;
        xp(i) = xp(i) + eps;
        fp = state_transition_ekf(xp, dt);
        F(:, i) = (fp - f0) / eps;
    end
end

