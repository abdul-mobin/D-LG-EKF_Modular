function tune_Q_two_params_oop()
    % Object-oriented wrapper for process-noise tuning.

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.tuneQTwoParams();
end