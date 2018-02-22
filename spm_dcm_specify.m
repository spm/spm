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
% Example:
%
%  s = struct();
%  s.name       = 'test';
%  s.u          = [1 1]';
%  s.delays     = [1.2 1.2];
%  s.TE         = 0.05;
%  s.nonlinear  = true;
%  s.two_state  = true;
%  s.stochastic = true;
%  s.centre     = true;
%  s.induced    = 0;
%  s.a          = a;
%  s.b          = b;
%  s.c          = c;
%  s.d          = d;
%  DCM = spm_dcm_specify(SPM,xY,s);
%__________________________________________________________________________
% Copyright (C) 2002-2017 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify.m 7264 2018-02-22 14:43:47Z peter $

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
