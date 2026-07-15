function nees = compute_nees(est_states, true_states, P_history)
    % Compute the Normalized Estimation Error Squared (NEES) for each time step.
    %
    % Inputs:
    %   est_states  - N x 5 matrix of estimated states
    %   true_states - N x 5 matrix of true states
    %   P_history   - 5 x 5 x N covariance history
    %
    % Output:
    %   nees        - N x 1 vector of NEES values

    N = size(est_states, 1);
    nees = zeros(N, 1);

    for k = 1:N
        err = (est_states(k, :) - true_states(k, :)).';
        P = squeeze(P_history(:, :, k));
        nees(k) = err' / P * err;
    end
end
