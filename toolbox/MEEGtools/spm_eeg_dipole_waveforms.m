function sD = spm_eeg_dipole_waveforms(S)
% Function for extracting source data using dipoles.
% FORMAT sD = spm_eeg_dipole_waveforms(S)
%
% S                    - input structure (optional)
% (optional) fields of S:
%   S.D                - MEEG object or filename of M/EEG mat-file
%   S.dipoles          - (optional)
%     Structure describing the dipoles
%     dipoles.pnt      - Nx3 matrix of locations in MNI coordinates
%     dipoles.ori      - Nx3 matrix of orientations in MNI coordinates
%     dipoles.label    - Nx1 cell array of dipole labels
%
% Output:
% sD                   - MEEG object (also written on disk)
% _______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%

% Vladimir Litvak
% $Id: spm_eeg_dipole_waveforms.m 6438 2015-05-18 11:50:42Z vladimir $

[Finter,Fgraph] = spm('FnUIsetup','Dipole waveform extraction', 0);
%%

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

[D, ok] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        errordlg(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        errordlg('The meeg file is corrupt or incomplete');
    end
    return
end


%% ============  Check the input dipole struct. Create a new one if necessary
if ~(isfield(S, 'dipoles') && isfield(S.dipoles, 'pnt') && isfield(S.dipoles, 'ori') && isfield(S.dipoles, 'label') && ...
        numel(S.dipoles.label)>0 && isequal(size(S.dipoles.pnt), [numel(S.dipoles.label), 3]) && ...
        isequal(size(S.dipoles.ori), size(S.dipoles.pnt)))
    S.dipoles = spm_eeg_dipoles_ui;
end

modality = spm_eeg_modality_ui(D, 1, 1);

if strncmp(modality, 'MEG', 3)
    reducerank = 2;
else
    reducerank = 3;
end

%% ============ Find or prepare head model

if ~isfield(D, 'val')
    D.val = 1;
end

if ~isfield(D, 'inv') || ~iscell(D.inv) ||...
        ~(isfield(D.inv{D.val}, 'forward') && isfield(D.inv{D.val}, 'datareg')) ||...
        ~isa(D.inv{D.val}.mesh.tess_ctx, 'char') % detects old version of the struct
    D = spm_eeg_inv_mesh_ui(D, D.val);
    D = spm_eeg_inv_datareg_ui(D, D.val);
    D = spm_eeg_inv_forward_ui(D, D.val);
end

data = spm_eeg_inv_get_vol_sens(D, [], [], [], modality);

vol  = data.(modality).vol;
sens = data.(modality).sens;

if isa(vol, 'char')
    vol = ft_read_vol(vol);
end

chanind = indchantype(D, modality, 'GOOD');

dipoles = [];
dipoles.pos = S.dipoles.pnt;
dipoles.ori = S.dipoles.ori;


dipoles = ft_transform_geometry(inv(data.transforms.toMNI), dipoles);
pnt = dipoles.pos;
label = S.dipoles.label;

ori = dipoles.ori;
for i = 1:numel(label)
    ori(i, :) = ori(i, :)/norm(ori(i, :));
end

[vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', D.chanlabels(chanind));

%% ============ Compute lead fields for the dipoles

nvert = numel(label);

Gxyz = ft_compute_leadfield(pnt, sens, vol, 'reducerank', reducerank, 'dipoleunit', 'nA*m', 'chanunit', D.units(chanind));

G = zeros(size(Gxyz, 1), size(Gxyz, 2)/3);
for i = 1:nvert
    G(:, i) = Gxyz(:, (3*i- 2):(3*i))*ori(i, :)';
end

%% ============  Use the montage functionality to compute source activity.

S = [];
S.D = D;

montage = [];
montage.tra = pinv(G); % This is the key line where the lead field is used to define the transformation
montage.labelorg = D.chanlabels(chanind);
montage.labelnew = label;

S.montage  = montage;

S.keepothers = false;
S.keepsensors = false;

sD = spm_eeg_montage(S);

sD = chantype(sD, ':', 'LFP');
sD = units(sD, ':', 'nA*m');

sD.save;
