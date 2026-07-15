function sweep_track_loss_by_t_end_oop()
    % Object-oriented wrapper for the track-loss sweep.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.sweepTrackLossByTEnd();
end