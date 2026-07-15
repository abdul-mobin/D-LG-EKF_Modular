function inspect_nees_distribution(nees)
    % Summarize the distribution of NEES values and highlight heavy-tail behavior.
    %
    % Inputs:
    %   nees - vector of NEES values
    %
    % Example:
    %   inspect_nees_distribution(nees_lgkf)

    if isempty(nees)
        error('nees must be a non-empty vector.');
    end

    fprintf('Mean: %.2e\n', mean(nees));
    fprintf('Median: %.2e\n', median(nees));
    fprintf('Max: %.2e\n', max(nees));
    fprintf('95th pct: %.2e\n', local_percentile(nees, 95));
    fprintf('99th pct: %.2e\n', local_percentile(nees, 99));
    fprintf('99.9th pct: %.2e\n', local_percentile(nees, 99.9));

    % Identify how many steps account for most of the total NEES sum.
    sorted_nees = sort(nees, 'descend');
    cumfrac = cumsum(sorted_nees) / sum(sorted_nees);
    n50 = find(cumfrac > 0.5, 1);
    fprintf('Top %d steps (%.4f%% of all steps) account for 50%% of the total NEES sum\n', ...
        n50, 100 * n50 / numel(nees));
end

function value = local_percentile(data, percentile)
    sorted_data = sort(data(:));
    if isempty(sorted_data)
        value = NaN;
        return;
    end

    if numel(sorted_data) == 1
        value = sorted_data(1);
        return;
    end

    position = 1 + (numel(sorted_data) - 1) * percentile / 100;
    lower_idx = floor(position);
    upper_idx = ceil(position);
    if lower_idx == upper_idx
        value = sorted_data(lower_idx);
    else
        weight = position - lower_idx;
        value = (1 - weight) * sorted_data(lower_idx) + weight * sorted_data(upper_idx);
    end
end
