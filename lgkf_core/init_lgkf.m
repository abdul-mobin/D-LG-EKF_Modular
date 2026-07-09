function [g, P] = init_lgkf(true_initial_state)
    % Initialize with the same offset as EKF for a fair comparison
    offset = [50; 50; deg2rad(2); 2; 0];
    initial_guess = true_initial_state + offset;
    
    % Convert vector to Lie Group matrix (using your utility function)
    g = vec2lie_radar(initial_guess);
    
    % Initial Covariance matrix
    P = diag([1e4, 1e4, 0.1, 100, 0.01]);
end