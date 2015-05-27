function spm_post_concatenate(P, scans)
% Adjusts an SPM.mat which has concatenated sessions to improve accuracy.
% The high pass filter will be re-specified as if sessions are separate.
%
% P     - filename of the SPM.mat file to adjust
% scans - [1 x n] vector with the original number of scans in each session
%
% The expected workflow is:
%
% 1. Manually specify a GLM with additional unconvolved regressors
%    for each session, except for the last session.
% 2. Run spm_post_concatenate on the saved SPM.mat.
% 3. Estimate the SPM.mat in the normal way.
%
% Tips:
%
% - The session regressors should be binary vectors. E.g. a regressor named
%   'sess1' would include 1s for every scan in session 1 and zeros
%   elsewhere.
%
% - The BOLD-response may overhang from one session to the next. To reduce
%   this, acquire additional volumes at the end of each session and / or
%   add regressors to model the trials at the session borders.
%
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging
%
% Guillaume Flandin & Peter Zeidman
% $Id: spm_post_concatenate.m 6432 2015-05-09 12:58:12Z karl $

% Validate input
if ~nargin || isempty(P)
    [P, sts] = spm_select(1,'^SPM.mat$','select SPM.mat');
end
if iscell(P), P = P{1}; end

% Check SPM is in fMRI mode
[modality,modnum]=spm('CheckModality');
if ~strcmp(modality,'FMRI')
    error('This function only works with SPM in fMRI mode.');
end

% Load
SPM = load(P);
try
    SPM = SPM.SPM;
catch
    error('Selected file is not a valid SPM.mat');
end

% Create backup
[spm_path,name,ext] = fileparts(P);
copyfile(P, fullfile(spm_path,'SPM_backup.mat'));

% -------------------------------------------------------------------------
% High-pass filter
SPM.nscan = scans;
s = cumsum([0 SPM.nscan]);

for i=1:numel(SPM.nscan)
    K(i) = struct('HParam', SPM.xX.K(1).HParam,...
                  'row',    s(i) + (1:SPM.nscan(i)),...
                  'RT',     SPM.xY.RT);
end
SPM.xX.K = spm_filter(K);

% -------------------------------------------------------------------------
% Temporal non-sphericity
SPM.xVi.Vi   = spm_Ce(SPM.nscan,0.2);
SPM.xVi.form = 'AR(0.2)';

% -------------------------------------------------------------------------
% Save
save(P, 'SPM', spm_get_defaults('mat.format'));
disp('SPM.mat adjusted for sessions: Please estimate to apply changes.');
