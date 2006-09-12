function D = spm_eeg_inv_elec_Rsph_ui(S)

% Deals with the solution of the forward problem:
% There are 2 seperate cases:
% 1. individual head model, based on the subject anatomy
% 2. standard head model, based on the template
%
% For an individual head model :
% - project the sensor location, on the scalp surface (EEG), not sure what
% to do for MEG...
% - Estimate the best fitting spheres, and prepare the realistic sphere
% approach
% For a standard head model :
% - load the head model built on the template in MNI space (actually the single subject
% canonical image)
% -select your standard electrode setup
%
% To do all this, I use bits of a routine previously written for the DipFit
% toolbox. I try to render things a bit more coherent...
%
% FORMAT D = spm_eeg_inv_elec_Rsph_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the meshing files and variables
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips
% $Id$

if nargin == 0
    D = spm_eeg_ldata;
elseif nargin == 1
    D = S;
else
    error(sprintf('Trouble reading the data file\n'));
end

% compute forward model
%__________________________________________________________________________
Mesh    = D.inv{D.val}.mesh;
Reg     = load(D.inv{D.val}.datareg.sens_coreg);
el_name = D.channels.name(D.channels.eeg);

if strcmp(questdlg('Use canonical mesh'),'Yes')
    
    % create head structure from previous tesselation  - inner skull
    %----------------------------------------------------------------------
    Vol              = spm_vol(Mesh.msk_iskull);
    Scalp            = load(Mesh.tess_iskull);

    [ctr, Centre]    = spm_eeg_inv_model('CtrBin',Mesh.msk_iskull);

    head.XYZmm       = Scalp.vert';
    head.tri         = Scalp.face';
    head.nr          = [Mesh.Iskull_Nv Mesh.Iskull_Nf];
    head.M           = Vol.mat;
    head.info.str    = 'Tess. inner skull';
    head.info.bin_fn = Vol.fname;
    head.Centre      = Centre;
    model.head       = head;

    % create head structure from previous tesselation  - scalp
    %----------------------------------------------------------------------
    Vol              = spm_vol(Mesh.msk_scalp);
    Scalp            = load(Mesh.tess_scalp);
    [ctr, Centre]    = spm_eeg_inv_model('CtrBin',Mesh.msk_scalp);

    head.XYZmm       = Scalp.vert';
    head.tri         = Scalp.face';
    head.nr          = [Mesh.Scalp_Nv Mesh.Scalp_Nf];
    head.M           = Vol.mat;
    head.info.str    = 'Tess. outer scalp';
    head.info.bin_fn = Vol.fname;
    head.Centre      = Centre;
    model.head(2)    = head;
    
else
    
    Pvol             = strvcat(Mesh.msk_iskull,Mesh.msk_scalp);
    model.head       = spm_eeg_inv_model('GenMesh',Pvol);
    
end

% Place the electrodes in the model
%--------------------------------------------------------------------------
[electrodes,f_el]    = spm_eeg_inv_model('Elec2Scalp',model.head(2), ...
                       Reg.sensreg',el_name);
model.electrodes     = electrodes;
model.flags.fl_elec  = f_el;

% [re]-set the sphere model
%--------------------------------------------------------------------------
try, model           = rmfield(model,'spheres'); end
[spheres,d,L,q,f_Rs] = spm_eeg_inv_Rsph(model,[]);
model.spheres        = spheres;
model.flags.fl_RealS = f_Rs;

% Save in D
D.inv{D.val}.forward = model;





