function DCM = spm_dcm_specify(SPM,xY,settings)
% Specify inputs of an fMRI DCM (wrapper around spm_dcm_specify_ui)
% FORMAT DCM = spm_dcm_specify(SPM,xY,settings)
%
% SPM      - SPM structure or its filename
% xY       - (optional) VOI structures to be inserted into the DCM
% settings - (optional) predefined configuration options
%
% DCM      - DCM structure (see spm_dcm_ui)
%
% Example for a task-based experiment:
% -------------------------------------------------------------------------
% n   = 3;    % number of regions
% nu  = 2;    % number of inputs (experimental conditions)
% TR  = 2;    % volume repetition time (seconds)
% TE  = 0.03; % echo time (seconds)
%
% % Experimental conditions to include from the SPM.
% % To see the conditions' names, load your SPM.mat into the workspace and 
% % inspect SPM.Sess(s).U.name, where s is the session (run) number.
% cond = struct();
% cond(1).name    = 'Condition1'; % desired name for the condition
% cond(1).spmname = {'c1','c2'};  % (optional) corresponding name(s) of 
%                                 % conditions in the SPM.mat file, see
%                                 % SPM.Sess(s).U.name. If multiple names 
%                                 % are provided then they will be combined 
%                                 % by binarizing the regressors and performing 
%                                 an 'OR' operation.
%
% cond(2).name    = 'Condition2';
% cond(2).spmname = {'c3','c4'}; 
% 
% % Connectivity matrices
% a  = ones(n,n);
% b  = zeros(n,n,nu);
% c  = ones(n,nu);
% d  = zeros(n,n,0);
%
% s = struct();
% s.name       = 'test';
% s.cond       = cond;
% s.delays     = repmat(TR/2, 1, n);
% s.TE         = TE;
% s.nonlinear  = false;
% s.two_state  = false;
% s.stochastic = false;
% s.centre     = true;
% s.induced    = 0;
% s.a          = a;
% s.b          = b;
% s.c          = c;
% s.d          = d;
% DCM = spm_dcm_specify(SPM,xY,s);
%
%
% Tips: 
% - You can either select which experimental conditions to include by using
% the s.cond structure, as illustrated above, or by specifying a matrix
% s.u(i,j), which sets whether to include regressor j of condition i from  
% the SPM design matrix. If there are no parametric regressors, then j 
% will always equal one.
%
% - xY is a cell array of strings containing the filenames of the VOIs to 
% include.
%
% Example for a resting state experiment:
% -------------------------------------------------------------------------
% n   = 2;    % number of regions
% nu  = 1;    % number of inputs. For DCM for CSD we have one input: null
% TR  = 2;    % volume repetition time (seconds)
% TE  = 0.03; % echo time (seconds)
% 
% % Connectivity matrices
% a  = ones(n,n);
% b  = zeros(n,n,nu);
% c  = zeros(n,nu);
% d  = zeros(n,n,0);
% 
% % Specify DCM
% s = struct();
% s.name       = model_name;
% s.u          = [];
% s.delays     = repmat(TR/2, 1, n);
% s.TE         = TE;
% s.nonlinear  = false;
% s.two_state  = false;
% s.stochastic = false;
% s.centre     = false;
% s.induced    = 1;       % indicates DCM for CSD
% s.a          = a;
% s.b          = b;
% s.c          = c;
% s.d          = d;
% 
% DCM = spm_dcm_specify(SPM,xY,s);
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2002-2022 Wellcome Centre for Human Neuroimaging

if nargin < 3, settings = struct(); end

%-Interactive window
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
spm_input('Specify DCM:...  ',1,'d');

%-Get design and directory
%--------------------------------------------------------------------------
if ~nargin || isempty(SPM)
    [SPM, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, DCM = []; return; end
end
if ischar(SPM)
    swd = spm_file(SPM,'fpath');
    try
        load(fullfile(swd,'SPM.mat'))
    catch
        error('Cannot read %s.',fullfile(swd,'SPM.mat'));
    end
    SPM.swd = swd;
else
    SPM.swd = pwd;
end

%-Name
%--------------------------------------------------------------------------
try
    name = settings.name;
catch
    name = spm_input('name for DCM_???.mat','+1','s');
end

%-Run specify UI to build DCM
%--------------------------------------------------------------------------
if nargin < 2
    xY = [];
end
DCM = spm_dcm_specify_ui(SPM, xY, settings);

%-Save
%--------------------------------------------------------------------------
save(fullfile(SPM.swd,['DCM_' name '.mat']),'DCM', spm_get_defaults('mat.format'));
