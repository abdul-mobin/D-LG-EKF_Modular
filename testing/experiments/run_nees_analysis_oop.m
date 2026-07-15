function run_nees_analysis_oop(n_trials, filter_type)
    % Object-oriented wrapper for NEES analysis.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.runNeesAnalysis(n_trials, filter_type);
end