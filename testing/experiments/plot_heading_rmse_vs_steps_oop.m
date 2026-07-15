function plot_heading_rmse_vs_steps_oop(n_trials, t_end, dt, az_std_deg)
    % Object-oriented wrapper for plotting heading RMSE versus step index.

    if nargin < 1 || isempty(n_trials)
        n_trials = [];
    end
    if nargin < 2 || isempty(t_end)
        t_end = [];
    end
    if nargin < 3 || isempty(dt)
        dt = [];
    end
    if nargin < 4 || isempty(az_std_deg)
        az_std_deg = [];
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarExperimentRunner();
    runner.plotHeadingRmseVsSteps(n_trials, t_end, dt, az_std_deg);
end