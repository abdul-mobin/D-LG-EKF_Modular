function Main_RADAR_Comparison_OOP()
    % Object-oriented entry point for the radar comparison experiment.

    repoRoot = fileparts(mfilename('fullpath'));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarExperimentRunner();
    config = runner.buildComparisonConfig();
    runner.runComparison(config);
end