function sweep_track_loss_by_t_end_oop(t_end_values, nTrials, threshold, fullLossGroupSteps)
    % Object-oriented wrapper for the track-loss sweep.

    if nargin < 1 || isempty(t_end_values)
        t_end_values = [];
    end
    if nargin < 2 || isempty(nTrials)
        nTrials = [];
    end
    if nargin < 3 || isempty(threshold)
        threshold = [];
    end
    if nargin < 4 || isempty(fullLossGroupSteps)
        fullLossGroupSteps = [];
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.sweepTrackLossByTEnd(t_end_values, nTrials, threshold, fullLossGroupSteps);
end