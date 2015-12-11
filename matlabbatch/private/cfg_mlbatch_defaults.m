function cfg_defaults = cfg_mlbatch_defaults

% function cfg_defaults = cfg_mlbatch_defaults
% This file contains defaults that control the behaviour and appearance 
% of matlabbatch.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_mlbatch_defaults.m 6641 2015-12-11 10:17:11Z volkmar $

rev = '$Rev: 6641 $'; %#ok

try
    % Font definition for cfg_ui user interface
    % cfg_defaults.cfg_ui.Xfont is a param/value list of property values as returned by uisetfont
    % lfont: used in lists, value edit dialogues etc.
    cfg_defaults.cfg_ui.lfont = {'FontAngle','normal',...
        'FontName',get(0,'factoryTextFontName'),...
        'FontSize',12,...
        'FontUnits','points',...
        'FontWeight','normal'};
    % bfont: used for buttons
    cfg_defaults.cfg_ui.bfont = {'FontAngle',get(0, 'factoryUicontrolFontAngle'),...
        'FontName',get(0,'factoryUicontrolFontName'),...
        'FontSize',get(0, 'factoryUicontrolFontSize'),...
        'FontUnits',get(0, 'factoryUicontrolFontUnits'),...
        'FontWeight',get(0, 'factoryUicontrolFontWeight')};
    % Toggle ExpertEdit mode. Value can be 'on' or 'off'
    cfg_defaults.cfg_ui.ExpertEdit = 'off';
catch
    cfg_defaults.cfg_ui = false;
end

% cfg_util
% Parallel execution of independent modules
% Currently, this does not run modules in parallel, but it may reorder
% execution order of modules: all modules without dependencies will be run
% before modules with dependencies will be harvested again. If some modules
% have side effects (e.g. "Change Directory") that are not encoded as
% dependency, this may lead to unwanted results. Disabling parallel
% execution incurs an overhead during job execution because the job
% must be harvested more often.
cfg_defaults.cfg_util.runparallel = false;
% cfg_util('genscript',...) hook to add application specific initialisation
% code to generated scripts. This must be a function handle that takes no
% input arguments and returns two output arguments - the code to be
% inserted as a string array, and a flag indicating whether cfg_util should
% add job execution code (flag == true) or not (flag == false).
% The generated code will be inserted after the for loop which collects the
% input. In the generated code, variables 'jobs' and 'inputs' can be
% referenced. These will hold the jobs and corresponding inputs.
cfg_defaults.cfg_util.genscript_run = [];
% A diary of command window output can be kept for jobs. This is useful for
% debugging various problems related to job execution.
cfg_defaults.cfg_util.run_diary = false;

% Message defaults
cfg_defaults.msgdef.identifier  = 'cfg_defaults:defaultmessage';
cfg_defaults.msgdef.level       = 'info'; % one of 'info', 'warning', 'error'
cfg_defaults.msgdef.destination = 'stdout'; % one of 'none', 'stdout',
                                            % 'stderr', 'syslog'. Errors
                                            % will always be logged to
                                            % the command window, and
                                            % additionally to syslog, if specified
cfg_defaults.msgdef.verbose     = 'off';
cfg_defaults.msgdef.backtrace   = 'off';

cfg_defaults.msgcfg(1)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(1).identifier  = 'matlabbatch:run:jobfailederr';
cfg_defaults.msgcfg(1).level       = 'error';
cfg_defaults.msgcfg(2)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(2).identifier  = 'matlabbatch:cfg_util:addapp:done';
cfg_defaults.msgcfg(2).destination = 'none';
cfg_defaults.msgcfg(3)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(3).identifier  = 'matlabbatch:initialise:invalid';
cfg_defaults.msgcfg(3).level       = 'error';
cfg_defaults.msgcfg(4)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(4).identifier  = 'matlabbatch:fopen';
cfg_defaults.msgcfg(4).level       = 'error';
cfg_defaults.msgcfg(5)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(5).identifier  = 'cfg_getfile:notcolumncell';
cfg_defaults.msgcfg(5).level       = 'error';
cfg_defaults.msgcfg(6)             = cfg_defaults.msgdef;
cfg_defaults.msgcfg(6).identifier  = 'matlabbatch:subsref:cfg_dep:multisubs';
cfg_defaults.msgcfg(6).level       = 'warning';
cfg_defaults.msgcfg(6).backtrace   = 'on';

cfg_defaults.msgtpl( 1)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 1).identifier  = '^matlabbatch:subsasgn';
cfg_defaults.msgtpl( 1).level       = 'error';
cfg_defaults.msgtpl( 2)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 2).identifier  = '^matlabbatch:subsref';
cfg_defaults.msgtpl( 2).level       = 'error';
cfg_defaults.msgtpl( 3)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 3).identifier  = '^matlabbatch:constructor';
cfg_defaults.msgtpl( 3).level       = 'error';
cfg_defaults.msgtpl( 4)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 4).identifier  = '^matlabbatch:deprecated';
cfg_defaults.msgtpl( 4).destination = 'none';
cfg_defaults.msgtpl( 5)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 5).identifier  = '^matlabbatch:usage';
cfg_defaults.msgtpl( 5).level       = 'error';
cfg_defaults.msgtpl( 6)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 6).identifier  = '^matlabbatch:setval';
cfg_defaults.msgtpl( 6).destination = 'none';
cfg_defaults.msgtpl( 7)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 7).identifier  = '^matlabbatch:run:nomods';
cfg_defaults.msgtpl( 7).level       = 'info';
cfg_defaults.msgtpl( 8)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 8).identifier  = '^matlabbatch:cfg_struct2cfg';
cfg_defaults.msgtpl( 8).destination = 'none';
cfg_defaults.msgtpl( 9)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl( 9).identifier  = '^MATLAB:inputdlg';
cfg_defaults.msgtpl( 9).level       = 'error';
cfg_defaults.msgtpl(10)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl(10).identifier  = '^MATLAB:listdlg';
cfg_defaults.msgtpl(10).level       = 'error';
cfg_defaults.msgtpl(11)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl(11).identifier  = '^MATLAB:num2str';
cfg_defaults.msgtpl(11).level       = 'error';
cfg_defaults.msgtpl(12)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl(12).identifier  = '^matlabbatch:ok_subsasgn';
cfg_defaults.msgtpl(12).destination = 'none';
cfg_defaults.msgtpl(13)             = cfg_defaults.msgdef;
cfg_defaults.msgtpl(13).identifier  = 'matlabbatch:checkval:numcheck:transposed';
cfg_defaults.msgtpl(13).destination = 'none';

% value check for cfg_branch/choice/repeat items - set to false after
% configuration has been initialised to speed up job
% initialisation/harvest/run - set this to true if you want to debug some
% configuration or the batch system itself
cfg_defaults.cfg_item.checkval = false;

% verbosity of code generation for cfg_dep objects
% true  - one-line code containing only source information
% false - full code, including source and target specifications
cfg_defaults.cfg_dep.gencode_short = true;
