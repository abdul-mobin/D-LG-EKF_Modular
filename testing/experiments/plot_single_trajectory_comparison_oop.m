function plot_single_trajectory_comparison_oop(traj_idx, az_std_deg)
    % Object-oriented wrapper for plotting one trajectory comparison.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.runSingleTrajectoryComparison(traj_idx, az_std_deg);
end