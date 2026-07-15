function run_heading_rmse_trials_oop(n_trials, t_end, dt, az_std_deg)
    % Object-oriented wrapper for the heading RMSE trial experiment.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarExperimentRunner();
    runner.runHeadingRmseTrials(n_trials, t_end, dt, az_std_deg);
end