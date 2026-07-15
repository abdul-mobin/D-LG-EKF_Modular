function tune_Q_two_params_oop(nTraj)
    % Object-oriented wrapper for process-noise tuning.

    if nargin < 1 || isempty(nTraj)
        nTraj = [];
    end

    repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repoRoot);
    addpath(genpath(repoRoot));

    runner = RadarDiagnosticsRunner();
    runner.tuneQTwoParams(nTraj);
end