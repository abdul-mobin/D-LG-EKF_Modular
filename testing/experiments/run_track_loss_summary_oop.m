function run_track_loss_summary_oop(n_trials, t_end, dt, az_std_deg, threshold, persist_steps)
    % Object-oriented wrapper for the track-loss summary experiment.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarExperimentRunner();
    runner.runTrackLossSummary(n_trials, t_end, dt, az_std_deg, threshold, persist_steps);
end