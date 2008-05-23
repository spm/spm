function cfg = cfg_example_master
% Master file that collects the cfg_exbranches in conceptually similar
% groups.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_master.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Collect Modules to add two Numbers
add        = cfg_repeat; % A repeat collects a variable number of items from its .values field in its .val field
add.name   = 'Add 2 Numbers';
add.tag    = 'add2';
add.values = {cfg_example_add1 cfg_example_add2}; % Reference the config files for both add modules
add.forcestruct = true; % There is a speciality in cfg_repeat harvest behaviour that makes a difference depending on the number of elements in values. forcestruct == true forces the repeat to behave as if there are more than one distinct values, even if there is only one.
add.help   = {'These modules add 2 numbers and return their sum'};

%% Collect Modules that sum over an Vector of arbitrary Length
sumv        = cfg_repeat;
sumv.name   = 'Sum Vectors';
sumv.tag    = 'sumv';
sumv.values = {cfg_example_sum cfg_example_cumsum1 cfg_example_cumsum2};
sumv.forcestruct = true;
sumv.help   = {'Sum over vectors.'};

%% Collect above Collections and the Div Module
cfg        = cfg_repeat;
cfg.name   = 'Toy Example';
cfg.tag    = 'cfg_toy';
cfg.values = {add sumv cfg_example_div}; % Values in a cfg_repeat can be any cfg_item objects
cfg.forcestruct = true;
cfg.help   = {'The full toy example'};