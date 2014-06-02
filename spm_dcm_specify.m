function DCM = spm_dcm_specify(SPM,xY)
% Specify inputs of a DCM (wrapper around spm_dcm_specify_ui)
% FORMAT DCM = spm_dcm_specify(SPM,xY)
%
% SPM      - SPM structure or its filename
% xY       - (optional) VOI structures to be inserted into the DCM
%
% DCM      - DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify.m 6031 2014-06-02 12:49:52Z guillaume $

%-Interactive window
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
spm_input('Specify DCM:...  ',1,'d');

% Get design and directory
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
else
    swd = pwd;
end

%-Name
%--------------------------------------------------------------------------
name = spm_input('name for DCM_???.mat','+1','s');

%-Run specify UI to build DCM
%--------------------------------------------------------------------------
if nargin < 2
    xY = [];
end
DCM = spm_dcm_specify_ui(SPM, xY); 

%-Save
%--------------------------------------------------------------------------
save(fullfile(swd,['DCM_' name '.mat']),'DCM', spm_get_defaults('mat.format'));
