function elec = spm_eeg_transform_cap(S)
% Transform an electrode cap to match the subject's headshape
% FORMAT shape = spm_eeg_transform_cap(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.standard         - headshape (file) with the standard locations
%   S.custom           - headshape (file) with individually measured locations
%   S.outfile          - file name to save the output                     
%
% Output:
%   sens               - transformed sensors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_transform_cap.m 4020 2010-07-28 12:41:23Z vladimir $

SVNrev = '$Rev: 4020 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Transform cap'); spm('Pointer','Watch');


%-Input
%--------------------------------------------------------------------------
if nargin == 0
    S = [];
end

if ~isfield(S, 'standard')
    [S.standard, sts] = spm_select(1, '.*', 'Select standard locations file');
    if ~sts, return; end
end

if ischar(S.standard)
    S.standard = ft_read_sens(S.standard);
end

if ~isfield(S, 'custom')
    [S.custom, sts] = spm_select(1, '.*', 'Select individual locations file');
    if ~sts, return; end
end

if ischar(S.custom)
    S.custom = ft_read_sens(S.custom);
end

S.standard = ft_convert_units(S.standard, 'mm');
S.custom = ft_convert_units(S.custom, 'mm');

%-Compute
%--------------------------------------------------------------------------
S1 = [];
S1.targetfid.pnt = [];
S1.targetfid.fid = S.custom;

S1.sourcefid.pnt = [];
S1.sourcefid.fid = S.standard;

S1.template = 1;

M1 = spm_eeg_inv_datareg(S1);

elec = ft_transform_sens(M1, S.standard);

%-Change SPM fiducial labels if present to prevent recognition as template.
%--------------------------------------------------------------------------
[sel1, sel2] = spm_match_str(elec.label, {'spmnas', 'spmlpa', 'spmrpa'});
fidlabels = {'nas', 'lpa', 'rpa'};
elec.label(sel1) = fidlabels(sel2);

%-Save
%--------------------------------------------------------------------------
if ~isfield(S, 'outfile')
    [f, p] = uiputfile( ...
        {'*.mat', 'MATLAB File (*.mat)'}, 'Save custom sensors as');
    S.outfile = fullfile(p, f);
end

save(S.outfile, 'elec');

%-Plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics');clf
ind = spm_match_str(elec.label, S.custom.label);
plot3(S.custom.pnt(:,1), S.custom.pnt(:,2), S.custom.pnt(:,3), '.r', 'MarkerSize', 20);
hold on
plot3(elec.pnt(:,1), elec.pnt(:,2), elec.pnt(:,3), '.k');
plot3(elec.pnt(ind, 1), elec.pnt(ind ,2), elec.pnt(ind ,3), '.g', 'MarkerSize', 20);
axis equal off
rotate3d on

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','Transform cap: done'); spm('Pointer','Arrow');
