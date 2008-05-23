function sum = cfg_example_sum
% Example script that creates an cfg_exbranch to sum two numbers. The
% inputs are entered as vector, the output is just a single
% number. This function differs from cfg_example_add2 (except from names)
% only in the specification of input1.num.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_sum.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Input Item
input1         = cfg_entry; % This is the generic data entry item
input1.name    = 'Input a'; % The displayed name
input1.tag     = 'a';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input1.strtype = 'e';       % No restriction on what type of data is entered. This could be used to restrict input to real numbers, integers ...
input1.num     = [1 Inf];   % Number of inputs required (2D-array with exactly one row and two column)
input1.help    = {'This is the input vector.','The elements will be added together.'}; % help text displayed

%% Executable Branch
sum      = cfg_exbranch;       % This is the branch that has information about how to run this module
sum.name = 'sum';             % The display name
sum.tag  = 'cfg_example_sum'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
sum.val  = {input1};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
sum.prog = @cfg_example_run_sum;  % A function handle that will be called with the harvested job to run the computation
sum.vout = @cfg_example_vout_sum; % A function handle that will be called with the harvested job to determine virtual outputs
sum.help = {'Add two numbers.'};

%% Local Functions
% The cfg_example_vout_sum function can go here, it is not useful outside
% the batch environment.
function vout = cfg_example_vout_sum(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

vout = cfg_dep;                        % The dependency object
vout.sname      = 'sum(a)';       % Displayed dependency name
vout.src_output = substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
