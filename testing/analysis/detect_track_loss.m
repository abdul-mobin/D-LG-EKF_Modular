function [loss_rate, loss_time] = detect_track_loss(est_xy, true_xy, T, threshold, persist_steps)
    % Detect whether a trajectory has effectively lost track.
    %
    % Inputs:
    %   est_xy        - N x 2 estimated position history
    %   true_xy       - N x 2 true position history
    %   T             - sampling interval [s] for est_xy/true_xy rows
    %   threshold     - position error threshold [m]
    %   persist_steps - number of consecutive samples above threshold required
    %                   before declaring loss
    %
    % Outputs:
    %   loss_rate     - 1 if the trajectory lost track, else 0
    %   loss_time     - first time index at which loss was declared, in seconds

    if nargin < 3 || isempty(T)
        T = 1;
    end
    if nargin < 4 || isempty(threshold)
        threshold = 1000;
    end
    if nargin < 5 || isempty(persist_steps)
        persist_steps = 20;
    end

    err = sqrt(sum((est_xy - true_xy).^2, 2));
    above = err > threshold;

    run_len = 0;
    loss_idx = -1;
    for k = 1:numel(above)
        if above(k)
            run_len = run_len + 1;
            if run_len >= persist_steps
                loss_idx = k - persist_steps + 1;
                break;
            end
        else
            run_len = 0;
        end
    end

    if loss_idx > 0
        loss_rate = 1;
        loss_time = loss_idx * T;
    else
        loss_rate = 0;
        loss_time = NaN;
    end
end
