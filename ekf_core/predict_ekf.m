function [x, P] = predict_ekf(x, P, Q, dt)
    v = x(4); w = x(5); th = x(3);

    if abs(w) > 1e-5
        f = [x(1) + (v/w)*(sin(th+w*dt) - sin(th)); x(2) + (v/w)*(cos(th) - cos(th+w*dt)); th+w*dt; v; w];
    else
        f = [x(1) + v*cos(th)*dt; x(2) + v*sin(th)*dt; th+w*dt; v; w];
    end

    % Analytical Jacobian for the CTRV model.
    F = eye(5);
    if abs(w) > 1e-5
        dx_dth = (v/w) * (cos(th - w*dt) - cos(th));
        dx_dv = (sin( w*dt-th) - sin(th)) / w;
        dx_dw = v * ((w*dt*cos(th + w*dt) - sin(th + w*dt) - sin(th)) / (w^2));

        dy_dth = (v/w) * (sin(th + w*dt) - sin(th));
        dy_dv = (cos(th) - cos(th + w*dt)) / w;
        dy_dw = v * ((w*dt*sin(th + w*dt) - cos(th) + cos(th + w*dt)) / (w^2));
    else
        dx_dth = -v * sin(th) * dt;
        dx_dv = cos(th) * dt;
        dx_dw = -v * sin(th) * (dt^2 / 2);

        dy_dth = v * cos(th) * dt;
        dy_dv = sin(th) * dt;
        dy_dw = v * cos(th) * (dt^2 / 2);
    end

    F(1,3) = dx_dth;
    F(1,4) = dx_dv;
    F(1,5) = dx_dw;
    F(2,3) = dy_dth;
    F(2,4) = dy_dv;
    F(2,5) = dy_dw;
    F(3,5) = dt;

    x = f; P = F * P * F' + Q;
end

