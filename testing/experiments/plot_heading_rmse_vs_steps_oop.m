function plot_heading_rmse_vs_steps_oop(n_trials, t_end, dt, az_std_deg)
    % Object-oriented wrapper for plotting heading RMSE versus step index.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarExperimentRunner();
    runner.plotHeadingRmseVsSteps(n_trials, t_end, dt, az_std_deg);
end