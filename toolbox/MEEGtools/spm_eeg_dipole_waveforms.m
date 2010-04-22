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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%

% Vladimir Litvak
% $Id: spm_eeg_dipole_waveforms.m 3833 2010-04-22 14:49:48Z vladimir $

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

[ok, D] = check(D, 'sensfid');

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

pnt = S.dipoles.pnt;
ori = S.dipoles.ori;
label = S.dipoles.label;

modality = spm_eeg_modality_ui(D, 1, 1);

if isequal(modality, 'MEG')
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

for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = ft_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m);
    end
end

sens = datareg.sensors;

M1 = datareg.toMNI;
[U, L, V] = svd(M1(1:3, 1:3));
M1(1:3,1:3) =U*V';

vol = ft_transform_vol(M1, vol);
sens = ft_transform_sens(M1, sens);

chanind = setdiff(meegchannels(D, modality), badchannels(D));

[vol, sens] = ft_prepare_vol_sens(vol, sens, 'channel', D.chanlabels(chanind));

%% ============ Compute lead fields for the dipoles

nvert = numel(label);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, ['Computing ' modality ' leadfields']); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = [1:nvert]; end

Gxyz = zeros(length(chanind), 3*nvert);
for i = 1:nvert

    Gxyz(:, (3*i- 2):(3*i))  = ft_compute_leadfield(pnt(i, :), sens, vol, 'reducerank', reducerank);

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');
spm_progress_bar('Init', nvert, ['Orienting ' modality ' leadfields']); drawnow;

G = zeros(size(Gxyz, 1), size(Gxyz, 2)/3);
for i = 1:nvert

    G(:, i) = Gxyz(:, (3*i- 2):(3*i))*ori(i, :)';

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');drawnow;

%% ============  Use the montage functionality to compute source activity.

S = [];
S.D = D;

montage = [];
montage.tra = pinv(G); % This is the key line where the lead field is used to define the transformation
montage.labelorg = D.chanlabels(chanind);
montage.labelnew = label;

S.montage  = montage;

S.keepothers = 'no';

sD = spm_eeg_montage(S);

sD = chantype(sD, [], 'LFP');

sD.save;
