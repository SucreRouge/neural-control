addpath(fullfile('..', 'utils'))
addpath(fullfile('..', 'visualization'))
addpath(fullfile('..', 'models'))

fnparams = FitzhughNagumoParams;
fnparams.A    = 0.5;
fnparams.B    = 0.25;
fnparams.Tau  = 1;
fnparams.Iext = 1;

fitzhugh_nagumo(TaskManager.FUN_SOL, fnparams)