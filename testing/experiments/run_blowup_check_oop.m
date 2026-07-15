function run_blowup_check_oop(filter_type, n_trials)
    % Object-oriented wrapper for the blow-up diagnostic.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.runBlowupCheck(filter_type, n_trials);
end