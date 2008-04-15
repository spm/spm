function [datareg, mesh] = spm_eeg_inv_template(Msize)
% Fill in mesh fields from saved template files:  This cicumvents the 
% need for a structural MRI and asumes the subject has, roughly the same
% shaped head as the template head.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_template.m 1406 2008-04-15 09:37:59Z vladimir $


% check for mesh size
%--------------------------------------------------------------------------
if nargin == 0
    str   = 'Mesh size (vertices)';
    Msize = spm_input(str,'+1','3000|4000|5000|7200',[1 2 3 4]);
end

mesh.Msize = Msize;

% SPM directory of canonical anatomy
%--------------------------------------------------------------------------
Cdir     = [spm('dir') filesep 'EEGtemplates'];

% fill in fields
%==========================================================================

% head model (sMRI)
%--------------------------------------------------------------------------
mesh.template   = 1;
mesh.sMRI       = fullfile(Cdir,'smri.img');
mesh.msk_iskull = fullfile(Cdir,'smri_iskull.img');
mesh.msk_scalp  = fullfile(Cdir,'smri_scalp.img');
mesh.msk_cortex = fullfile(Cdir,'smri_cortex.img');


% meshes
%--------------------------------------------------------------------------
switch mesh.Msize
    case 1
        Tmesh = load(fullfile(Cdir,'wmeshTemplate_3004d.mat'));
    case 2
        Tmesh = load(fullfile(Cdir,'wmeshTemplate_4004d.mat'));
    case 3
        Tmesh = load(fullfile(Cdir,'wmeshTemplate_5004d.mat'));
    case 4
        Tmesh = load(fullfile(Cdir,'wmeshTemplate_7204d.mat'));
end

% Canonical cortical mesh
%----------------------------------------------------------------------
mesh.tess_mni.vert    = Tmesh.vert;
mesh.tess_mni.face    = uint16(Tmesh.face);

% Cortical mesh from the template
%----------------------------------------------------------------------
mesh.tess_ctx.vert    = Tmesh.vert;
mesh.tess_ctx.face    = uint16(Tmesh.face);

% Scalp mesh from the template
%----------------------------------------------------------------------
Tmesh      = load(fullfile(Cdir,'wmeshTemplate_scalp.mat'));
mesh.tess_scalp.vert  = Tmesh.vert;
mesh.tess_scalp.face  = uint16(Tmesh.face);

% Skull mesh from the template
%----------------------------------------------------------------------
Tmesh      = load(fullfile(Cdir,'wmeshTemplate_skull.mat'));
mesh.tess_iskull.vert = Tmesh.vert;
mesh.tess_iskull.face = uint16(Tmesh.face);


% datareg
%--------------------------------------------------------------------------
datareg.fid_mri    = [ 0   80  -46;
                     -79  -12  -63;
                      79  -12  -63];
datareg.scalpvert  = mesh.tess_scalp.vert;
datareg = datareg;

