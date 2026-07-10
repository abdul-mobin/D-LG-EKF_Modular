% test_ekf_jacobians.m
% Compare analytical EKF Jacobians against numerical finite-difference Jacobians.

clearvars; clc;
repoRoot = fileparts(mfilename('fullpath'));
addpath(genpath(repoRoot));

states = [
    1.0, 2.0, 0.3, 4.0, 0.1;
    1.0, 2.0, 0.3, 4.0, 1e-6;
    5.0, -3.0, -1.2, 2.5, 0.5;
];
dt = 0.1;

fprintf('EKF Jacobian comparison test\n');
for idx = 1:size(states,1)
    x = states(idx, :)';
    fprintf('\nState %d: [%s]\n', idx, num2str(x', ' %.4f'));

    F_num = numerical_F(x, dt);
    F_ana = analytical_F(x, dt);
    H_num = numerical_H(x);
    H_ana = analytical_H(x);

    fprintf('F diff max abs = %.3e\n', max(abs(F_num(:) - F_ana(:))));
    fprintf('H diff max abs = %.3e\n', max(abs(H_num(:) - H_ana(:))));
    fprintf('F numeric:\n'); disp(F_num);
    fprintf('F analytic:\n'); disp(F_ana);
    fprintf('H numeric:\n'); disp(H_num);
    fprintf('H analytic:\n'); disp(H_ana);
end

function F = analytical_F(x, dt)
    v = x(4); w = x(5); th = x(3);
    F = eye(5);
    if abs(w) > 1e-5
        dx_dth = (v/w) * (cos(th + w*dt) - cos(th));
        dx_dv = (sin(th + w*dt) - sin(th)) / w;
        dx_dw = v * ((w*dt*cos(th + w*dt) - sin(th + w*dt) + sin(th)) / (w^2));

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
end

function H = analytical_H(x)
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
end

function F = numerical_F(x, dt)
    eps = 1e-6;
    f0 = state_transition(x, dt);
    F = zeros(5);
    for i = 1:5
        xp = x;
        xp(i) = xp(i) + eps;
        fn = state_transition(xp, dt);
        F(:, i) = (fn - f0) / eps;
    end
end

function H = numerical_H(x)
    eps = 1e-6;
    h0 = measurement_model(x);
    H = zeros(3,5);
    for i = 1:5
        xp = x;
        xp(i) = xp(i) + eps;
        hp = measurement_model(xp);
        H(:, i) = (hp - h0) / eps;
    end
end

function f = state_transition(x, dt)
    v = x(4); w = x(5); th = x(3);
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

function h = measurement_model(x)
    alpha = atan2(x(2), x(1));
    rho = sqrt(x(1)^2 + x(2)^2);
    h = [alpha; rho; x(4) * cos(alpha - x(3))];
end