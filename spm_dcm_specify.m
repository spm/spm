function DCM = spm_dcm_specify(SPM,xY)
% Specify inputs of a DCM (wrapper around spm_dcm_specify_ui)
% FORMAT [DCM] = spm_dcm_specify(SPM,xY)
%
% xY       - (optional) VOI structures to be inserted into the DCM
%
% Returns:
% DCM      - the DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2002-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify.m 6026 2014-05-30 11:09:03Z peter $

%-Interactive window
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
bcolor = get(Finter,'color');
WS     = spm('WinScale');
dx     = 20;

spm_input('Specify DCM:...  ',1,'d');

%==========================================================================
% Get design and directory
%==========================================================================
if ~nargin || isempty(SPM)
    [SPM, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, DCM = []; return; end
end
if ischar(SPM)
    swd = spm_file(SPM,'fpath');
    try
        load(fullfile(swd,'SPM.mat'))
    catch
        error(['Cannot read ' fullfile(swd,'SPM.mat')]);
    end
else
    swd = pwd;
end

%==========================================================================
% Name
%==========================================================================
name = spm_input('name for DCM_???.mat','+1','s');

%==========================================================================
% Run specify UI to build DCM
%==========================================================================
if nargin < 2
    xY = [];
end

DCM = spm_dcm_specify_ui(SPM, xY); 

%==========================================================================
% Save
%==========================================================================
save(fullfile(swd,['DCM_' name '.mat']),'DCM', spm_get_defaults('mat.format'));