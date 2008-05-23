function cumsum = cfg_example_cumsum2
% Example script that creates an cfg_exbranch to sum two numbers. The
% inputs are entered as vector, the output is a vector containing the
% cumulative sums. This function differs from cfg_example_sum (except from
% names) only in the specification of the output subscript.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_cumsum2.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Input Item
input1         = cfg_entry; % This is the generic data entry item
input1.name    = 'Input a'; % The displayed name
input1.tag     = 'a';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input1.strtype = 'e';       % No restriction on what type of data is entered. This could be used to restrict input to real numbers, integers ...
input1.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and two column)
input1.help    = {'This is a vector element.'}; % help text displayed

%% Collect Single Numbers to Vector
invec        = cfg_repeat;
invec.name   = 'Vector';
invec.tag    = 'unused'; % According to the harvest rules for cfg_repeat, this tag will never appear in the output, because there is only one cfg_item to repeat and forcestruct is false
invec.values = {input1};
invec.num    = [1 Inf];  % At least one input should be there
invec.help   = {'Enter a vector element by element.'};

%% Executable Branch
cumsum      = cfg_exbranch;       % This is the branch that has information about how to run this module
cumsum.name = 'cumsum2';             % The display name
cumsum.tag  = 'cfg_example_cumsum2'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
cumsum.val  = {invec};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
cumsum.prog = @cfg_example_run_cumsum2;  % A function handle that will be called with the harvested job to run the computation
cumsum.vout = @cfg_example_vout_cumsum2; % A function handle that will be called with the harvested job to determine virtual outputs
cumsum.help = {'Compute the cumulative sum of a vector.'};

%% Local Functions
% The cfg_example_vout_cumsum2 function can go here, it is not useful outside
% the batch environment.
function vout = cfg_example_vout_cumsum2(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.
for k = 1:numel(job.a)
    % Create a subscript reference for each element of the cumsum vector
    vout(k) = cfg_dep;                               % The dependency object
    vout(k).sname      = sprintf('cumsum(a(1:%d))', k); % Displayed dependency name
    vout(k).src_output = substruct('.','cs','()',{k}); % The output subscript reference. The length of the output vector depends on the unknown length of the input vector. Therefore, the generic ':' subscript is needed.
end;