function check_for_blowup(mu_history, P_history, theta_history)
    % Diagnose numerical blow-up or covariance issues during filtering.
    %
    % Inputs:
    %   mu_history   - state history matrix of size N x 5 (or N x 5 x 1-style)
    %   P_history    - covariance history of size 5 x 5 x N
    %   theta_history- scalar history of heading values

    fprintf('Any NaN in mu? %d\n', any(isnan(mu_history(:))));
    fprintf('Any Inf in mu? %d\n', any(isinf(mu_history(:))));
    fprintf('Any NaN in P? %d\n', any(isnan(P_history(:))));
    fprintf('Any Inf in P? %d\n', any(isinf(P_history(:))));

    % Check P for PSD violations at every step.
    N = size(P_history, 3);
    min_eig = zeros(N, 1);
    for k = 1:N
        min_eig(k) = min(eig(squeeze(P_history(:, :, k))));
    end
    fprintf('Most negative eigenvalue seen: %.6e\n', min(min_eig));

    % Check for steps with tiny theta values.
    small_theta_frac = mean(abs(theta_history) < 1e-3);
    fprintf('Fraction of steps with |theta| < 1e-3: %.2f%%\n', small_theta_frac * 100);

    % Inspect the magnitude of P and mu over time.
    figure('Color', 'w', 'Position', [200, 200, 900, 600]);
    subplot(2, 1, 1);
    plot(squeeze(P_history(1, 1, :)));
    title('P(1,1) over time');
    xlabel('Step');
    ylabel('P(1,1)');

    subplot(2, 1, 2);
    if size(mu_history, 2) >= 2
        plot(vecnorm(mu_history(:, 1:2), 2, 2));
    else
        plot(mu_history(:, 1));
    end
    title('||position|| over time');
    xlabel('Step');
    ylabel('Position norm');
end
