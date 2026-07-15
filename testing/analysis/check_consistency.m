function check_consistency(nees, dim, alpha)
    % Check whether the observed NEES values are statistically consistent.
    %
    % Inputs:
    %   nees - vector of NEES values
    %   dim  - state dimension used for the chi-squared reference
    %   alpha - significance level (default 0.05)

    if nargin < 3, alpha = 0.05; end

    lower = chi2inv(alpha/2, dim);
    upper = chi2inv(1 - alpha/2, dim);
    frac_in = mean(nees >= lower & nees <= upper);

    fprintf('NEES in-bound fraction: %.2f%% (target ~%.0f%%)\n', frac_in * 100, (1 - alpha) * 100);
    fprintf('Mean NEES: %.2f (ideal = %d)\n', mean(nees), dim);
end
