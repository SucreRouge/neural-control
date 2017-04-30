addpath('/home/g/Aether/not_lab/src/neural_stimulation_models/utils')
addpath('/home/g/Aether/not_lab/src/neural_stimulation_models/visualization')
addpath('/home/g/Aether/not_lab/src/neural_stimulation_models/models')

fnparams = FitzhughNagumoParams;
fnparams.A    = 0.5;
fnparams.B    = 0.25;
fnparams.Tau  = 1;
fnparams.Iext = 1;

fitzhugh_nagumo(TaskManager.FUN_SOL, fnparams)