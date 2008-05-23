function div = cfg_example_div
% Example script that creates an cfg_exbranch to compute mod and rem of two
% natural numbers. The inputs are entered as two single numbers, the output
% is a struct with two fields 'mod' and 'rem'.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_div.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Input Items
% Input a
input1         = cfg_entry; % This is the generic data entry item
input1.name    = 'Input a'; % The displayed name
input1.tag     = 'a';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input1.strtype = 'n';       % Natural number required
input1.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
input1.help    = {'This is input a.','This input will be divided by the other input.'}; % help text displayed

% Input b
input2         = cfg_entry; % This is the generic data entry item
input2.name    = 'Input b'; % The displayed name
input2.tag     = 'b';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input2.strtype = 'n';       % Natural number required
input2.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
input2.help    = {'This is input b.','The other input will be divided by this input.'}; % help text displayed

%% Executable Branch
div      = cfg_exbranch;       % This is the branch that has information about how to run this module
div.name = 'div';             % The display name
div.tag  = 'cfg_example_div'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
div.val  = {input1 input2};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
div.prog = @cfg_example_run_div;  % A function handle that will be called with the harvested job to run the computation
div.vout = @cfg_example_vout_div; % A function handle that will be called with the harvested job to determine virtual outputs
div.help = {'Compute mod and rem of two numbers.'};

%% Local Functions
% The cfg_example_vout_div function can go here, it is not useful outside
% the batch environment.
function vout = cfg_example_vout_div(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

vout(1) = cfg_dep;                         % The dependency object
vout(1).sname      = 'a div b: mod';       % Displayed dependency name
vout(1).src_output = substruct('.','mod'); % The output subscript reference. 
vout(2) = cfg_dep;                         % The dependency object
vout(2).sname      = 'a div b: rem';       % Displayed dependency name
vout(2).src_output = substruct('.','rem'); % The output subscript reference. 
