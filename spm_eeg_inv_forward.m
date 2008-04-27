function D = spm_eeg_inv_forward(varargin)
% FORMAT D = spm_eeg_inv_forward(D,val)
%
%
% D                - input struct
% (optional) fields of S:
% D                - filename of EEG/MEG mat-file
%
% Output:
% D                - EEG/MEG struct with filenames of Gain matrices)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_forward.m 1488 2008-04-27 14:11:48Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});


if strcmp(D.inv{val}.modality, 'MEG')
    error('MEG support is under construction');
end

% Head Geometry (create tesselation file)
%--------------------------------------------------------------------------
[path,nam,ext]             = fileparts(D.inv{val}.mesh.sMRI);
vert = D.inv{val}.mesh.tess_ctx.vert;
face = D.inv{val}.mesh.tess_ctx.face;

% normals
%--------------------------------------------------------------------------
norm = spm_eeg_inv_normals(vert,face);

vol = D.inv{val}.forward.vol;
sens = D.inv{val}.datareg.sensors;

% Forward computation
%--------------------------------------------------------------------------
[vol, sens] = prepare_vol_sens(vol, sens, 'channel', D.inv{val}.forward.channels);
Gxyz  = compute_leadfield(vert, sens, vol) ;
%%
G = [];
for i = 1:size(norm, 1)
    G = [G Gxyz(:, (3*i- 2):(3*i))*norm(i, :)'];
end


% Save
%--------------------------------------------------------------------------
D.inv{val}.forward.gainmat = fullfile(D.path,[nam '_SPMgainmatrix_' num2str(val) '.mat']);
D.inv{val}.forward.gainxyz = fullfile(D.path,[nam '_SPMgainmatxyz_' num2str(val) '.mat']);

save(D.inv{val}.forward.gainmat,'G');
save(D.inv{val}.forward.gainxyz,'Gxyz');

