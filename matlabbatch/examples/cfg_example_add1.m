function add1 = cfg_example_add1
% Example script that creates an cfg_exbranch to sum two numbers. The
% inputs are entered as two single numbers, the output is just a single
% number.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_example_add1.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

%% Input Items
% Input a
input1         = cfg_entry; % This is the generic data entry item
input1.name    = 'Input a'; % The displayed name
input1.tag     = 'a';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input1.strtype = 'e';       % No restriction on what type of data is entered. This could be used to restrict input to real numbers, integers ...
input1.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
input1.help    = {'This is input a.','This input will be added to the other input.'}; % help text displayed

% Input b
input2         = cfg_entry; % This is the generic data entry item
input2.name    = 'Input b'; % The displayed name
input2.tag     = 'b';       % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
input2.strtype = 'e';       % No restriction on what type of data is entered. This could be used to restrict input to real numbers, integers ...
input2.num     = [1 1];     % Number of inputs required (2D-array with exactly one row and one column)
input2.help    = {'This is input b.','This input will be added to the other input.'}; % help text displayed

%% Executable Branch
add1      = cfg_exbranch;       % This is the branch that has information about how to run this module
add1.name = 'Add1';             % The display name
add1.tag  = 'cfg_example_add1'; % The name appearing in the harvested job structure. This name must be unique among all items in the val field of the superior node
add1.val  = {input1 input2};    % The items that belong to this branch. All items must be filled before this branch can run or produce virtual outputs
add1.prog = @cfg_example_run_add1;  % A function handle that will be called with the harvested job to run the computation
add1.vout = @cfg_example_vout_add1; % A function handle that will be called with the harvested job to determine virtual outputs
add1.help = {'Add two numbers.'};

%% Local Functions
% The cfg_example_vout_add1 function can go here, it is not useful outside
% the batch environment.
function vout = cfg_example_vout_add1(job)
% Determine what outputs will be present if this job is run. In this case,
% the structure of the inputs is fixed, and the output is always a single
% number. Note that input items may not be numbers, they can also be
% dependencies.

vout = cfg_dep;                        % The dependency object
vout.sname      = 'Add1: a + b';       % Displayed dependency name
vout.src_output = substruct('()',{1}); % The output subscript reference. This could be any reference into the output variable created during computation
