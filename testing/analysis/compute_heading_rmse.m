function rmse = compute_heading_rmse(est_theta, true_theta)
    % Compute heading RMSE with angle wrapping to [-pi, pi].
    %
    % Inputs can be numeric vectors/matrices with matching size, or cell arrays
    % where each cell contains one trajectory.

    if iscell(est_theta)
        if ~iscell(true_theta) || numel(est_theta) ~= numel(true_theta)
            error('For cell input, est_theta and true_theta must be cell arrays with equal length.');
        end

        sq_err_all = [];
        for i = 1:numel(est_theta)
            est_i = est_theta{i};
            true_i = true_theta{i};
            if ~isequal(size(est_i), size(true_i))
                error('Cell entry %d has mismatched heading sizes.', i);
            end
            err_i = wrap_to_pi(est_i - true_i);
            sq_err_all = [sq_err_all; err_i(:).^2]; %#ok<AGROW>
        end
        rmse = sqrt(mean(sq_err_all));
        return;
    end

    if ~isequal(size(est_theta), size(true_theta))
        error('est_theta and true_theta must have matching size.');
    end

    err = wrap_to_pi(est_theta - true_theta);
    rmse = sqrt(mean(err(:).^2));
end

function theta_wrapped = wrap_to_pi(theta)
    theta_wrapped = mod(theta + pi, 2 * pi) - pi;
end