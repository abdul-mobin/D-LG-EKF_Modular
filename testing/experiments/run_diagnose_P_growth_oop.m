function run_diagnose_P_growth_oop(filter_type)
    % Object-oriented wrapper for covariance-growth diagnostics.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.runDiagnosePGrowth(filter_type);
end