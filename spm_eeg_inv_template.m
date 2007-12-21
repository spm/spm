function [D] = spm_eeg_inv_template(varargin);

%==========================================================================
% Fill in mesh fields from saved template files:  This cicumvents the 
% need for a structural MRI and asumes the subject has, roughly the same
% shaped head as the template head.
%
% FORMAT D = spm_eeg_inv_template(D,[val])
% Input:
% D		    - input data struct (optional)
% Output:
% D			- same data struct including the forward solution files and variables
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_template.m 1039 2007-12-21 20:20:38Z karl $

% checks
%--------------------------------------------------------------------------
[D, val] = spm_eeg_inv_check(varargin{:});

% check for mesh size
%--------------------------------------------------------------------------
try
    mesh.Msize = D.inv{val}.mesh.Msize;
    D.inv{val}.mesh.Msize = mesh.Msize;
catch
    str   = 'Mesh size (vertices)';
    D.inv{val}.mesh.Msize = spm_input(str,'+1','3000|4000|5000|7200',[1 2 3 4]);
end


% SPM directory of canonical anatomy
%--------------------------------------------------------------------------
Cdir     = [spm('dir') filesep 'EEGtemplates'];

% fill in fields
%==========================================================================

% head model (sMRI)
%--------------------------------------------------------------------------
D.inv{val}.mesh.template   = 1;
D.inv{val}.mesh.sMRI       = fullfile(Cdir,'smri.img');
D.inv{val}.mesh.msk_iskull = fullfile(Cdir,'smri_iskull.img');
D.inv{val}.mesh.msk_scalp  = fullfile(Cdir,'smri_scalp.img');
D.inv{val}.mesh.msk_cortex = fullfile(Cdir,'smri_cortex.img');


% meshes
%--------------------------------------------------------------------------
switch D.inv{val}.mesh.Msize
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
D.inv{val}.mesh.tess_mni.vert    = Tmesh.vert;
D.inv{val}.mesh.tess_mni.face    = uint16(Tmesh.face);

% Cortical mesh from the template
%----------------------------------------------------------------------
D.inv{val}.mesh.tess_ctx.vert    = Tmesh.vert;
D.inv{val}.mesh.tess_ctx.face    = uint16(Tmesh.face);

% Scalp mesh from the template
%----------------------------------------------------------------------
Tmesh      = load(fullfile(Cdir,'wmeshTemplate_scalp.mat'));
D.inv{val}.mesh.tess_scalp.vert  = Tmesh.vert;
D.inv{val}.mesh.tess_scalp.face  = uint16(Tmesh.face);

% Skull mesh from the template
%----------------------------------------------------------------------
Tmesh      = load(fullfile(Cdir,'wmeshTemplate_skull.mat'));
D.inv{val}.mesh.tess_iskull.vert = Tmesh.vert;
D.inv{val}.mesh.tess_iskull.face = uint16(Tmesh.face);


% datareg
%--------------------------------------------------------------------------
datareg.fid_mri    = [ 0   80  -46;
                     -79  -12  -63;
                      79  -12  -63];
datareg.scalpvert  = D.inv{val}.mesh.tess_scalp.vert;
D.inv{val}.datareg = datareg;

