function diagnose_P_growth(mu_history, P_history, true_states, T)
    % Diagnose covariance growth and compare it to observed estimation error.
    %
    % Inputs:
    %   mu_history - N x 5 estimated state history in the same coordinate frame
    %                as the true state vector [x, y, theta, v, omega]
    %   P_history  - 5 x 5 x N covariance history
    %   true_states- N x 5 true state history
    %   T          - sampling interval
    %
    % Example:
    %   diagnose_P_growth(mu_hist, P_hist, true_states, config.dt)

    if nargin < 4 || isempty(T)
        T = 1;
    end

    N = size(P_history, 3);
    if size(mu_history, 1) ~= N || size(true_states, 1) ~= N
        error('mu_history and true_states must contain one row per time step.');
    end

    P_diag = zeros(N, 5);
    for k = 1:N
        P_diag(k, :) = diag(squeeze(P_history(:, :, k)))';
    end

    time = (0:N-1) * T;

    figure('Color', 'w', 'Position', [200, 200, 1000, 800]);
    labels = {'x', 'y', 'theta', 'v', 'omega'};
    for i = 1:5
        subplot(5, 1, i);
        semilogy(time, max(P_diag(:, i), 1e-12));
        title(['P(' num2str(i) ',' num2str(i) ') -- ' labels{i}]);
        xlabel('Time [s]');
        ylabel('Variance');
        grid on;
    end

    err = mu_history - true_states;
    err(:, 3) = wrap_angle(err(:, 3));

    fprintf('\n--- Actual error^2 vs reported P diagonal, at final step ---\n');
    for i = 1:5
        fprintf('%s: err^2 = %.4e, P(%d,%d) = %.4e, ratio = %.2e\n', ...
            labels{i}, err(end, i)^2, i, i, P_diag(end, i), err(end, i)^2 / max(P_diag(end, i), 1e-12));
    end

    fprintf('\n--- P(1,1) growth: first vs last step ---\n');
    fprintf('P(1,1) at k=1: %.4e\n', P_diag(1, 1));
    fprintf('P(1,1) at k=end: %.4e\n', P_diag(end, 1));
    fprintf('Growth factor: %.4e\n', P_diag(end, 1) / max(P_diag(1, 1), 1e-12));
end

function theta_wrapped = wrap_angle(theta)
    theta_wrapped = mod(theta + pi, 2*pi) - pi;
end
