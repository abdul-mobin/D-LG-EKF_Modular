function [x, P] = init_ekf(true_initial_state)
    % Initialize with some offset from truth to simulate real-world uncertainty
    % [x; y; theta; v; omega]
    offset = [50; 50; deg2rad(2); 2; 0];
    x = true_initial_state + offset;
    
    % Initial Covariance matrix
    P = diag([1e4, 1e4, 0.1, 100, 0.01]);
end