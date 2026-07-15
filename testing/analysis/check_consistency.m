function check_consistency(nees, dim, alpha)
    % Check whether the observed NEES values are statistically consistent.
    %
    % Inputs:
    %   nees - vector of NEES values
    %   dim  - state dimension used for the chi-squared reference
    %   alpha - significance level (default 0.05)

    if nargin < 3, alpha = 0.05; end

    lower = chi2inv_fallback(alpha/2, dim);
    upper = chi2inv_fallback(1 - alpha/2, dim);
    frac_in = mean(nees >= lower & nees <= upper);

    fprintf('NEES in-bound fraction: %.2f%% (target ~%.0f%%)\n', frac_in * 100, (1 - alpha) * 100);
    fprintf('Mean NEES: %.2f (ideal = %d)\n', mean(nees), dim);
end

function x = chi2inv_fallback(p, k)
    if exist('chi2inv', 'file') == 2
        x = chi2inv(p, k);
        return;
    end

    p = min(max(p, eps), 1 - eps);
    z = sqrt(2) * erfinv(2 * p - 1);
    x = k * (1 - 2 / (9 * k) + z * sqrt(2 / (9 * k)))^3;
    x = max(x, 0);
end
