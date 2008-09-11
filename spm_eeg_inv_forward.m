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
% $Id: spm_eeg_inv_forward.m 2081 2008-09-11 13:04:24Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});


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
[vol, sens] = forwinv_prepare_vol_sens(vol, sens, 'channel', D.inv{val}.forward.channels);

nvert = size(vert, 1);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, 'Computing leadfields'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = [1:nvert]; end

Gxyz = zeros(length(D.inv{val}.forward.channels), 3*nvert);
for i = 1:nvert

    Gxyz(:, (3*i- 2):(3*i))  = forwinv_compute_leadfield(vert(i, :), sens, vol);

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');
%%
spm_progress_bar('Init', nvert, 'Orienting leadfields'); drawnow;

G = zeros(size(Gxyz, 1), size(Gxyz, 2)/3);
for i = 1:nvert
    
    G(:, i) = Gxyz(:, (3*i- 2):(3*i))*norm(i, :)';

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');drawnow;

% Save
%--------------------------------------------------------------------------
D.inv{val}.forward.gainmat = fullfile(D.path,[nam '_SPMgainmatrix_' num2str(val) '.mat']);
%D.inv{val}.forward.gainxyz = fullfile(D.path,[nam '_SPMgainmatxyz_' num2str(val) '.mat']);

save(D.inv{val}.forward.gainmat,'G');
%save(D.inv{val}.forward.gainxyz,'Gxyz');

