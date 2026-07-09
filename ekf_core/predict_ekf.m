function [x, P] = predict_ekf(x, P, Q, dt)
    v = x(4); w = x(5); th = x(3);
    if abs(w) > 1e-5
        f = [x(1) + (v/w)*(sin(th+w*dt) - sin(th)); x(2) + (v/w)*(cos(th) - cos(th+w*dt)); th+w*dt; v; w];
    else
        f = [x(1) + v*cos(th)*dt; x(2) + v*sin(th)*dt; th+w*dt; v; w];
    end
    % Numerical Jacobian
    F = eye(5); eps = 1e-5;
    for i=1:5
        xp = x; xp(i) = xp(i) + eps;
        vp = xp(4); wp = xp(5); thp = xp(3);
        if abs(wp) > 1e-5
            fn = [xp(1) + (vp/wp)*(sin(thp+wp*dt) - sin(thp)); xp(2) + (vp/wp)*(cos(thp) - cos(thp+wp*dt)); thp+wp*dt; vp; wp];
        else
            fn = [xp(1) + vp*cos(thp)*dt; xp(2) + vp*sin(thp)*dt; thp+wp*dt; vp; wp];
        end
        F(:,i) = (fn - f) / eps;
    end
    x = f; P = F * P * F' + Q;
end

