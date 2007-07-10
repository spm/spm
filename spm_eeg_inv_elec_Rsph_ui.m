function D = spm_eeg_inv_elec_Rsph_ui(varargin)

% Deals with the solution of the forward problem:
% There are 2 seperate cases:
% 1. individual head model, based on the subject anatomy
% 2. standard head model,   based on the template
%
% For an individual head model :
% - project the sensor location, on the scalp surface (EEG)
% - Estimate the best fitting spheres, and prepare the realistic sphere
% approach
% For a standard head model :
% - load the head model built on the template in MNI space (actually the single subject
% canonical image)
% -select your standard electrode setup
% This is not currently implemented but see prgrmamming notes at the end of
% this script
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

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});

% compute forward model
%==========================================================================
Mesh    = D.inv{D.val}.mesh;
Reg     = D.inv{D.val}.datareg.sens_coreg;
el_name = D.channels.name(setdiff(D.channels.eeg, D.channels.Bad));

try 
    canonical = D.inv{D.val}.mesh.canonical;
catch
    D.inv{D.val}.mesh.canonical = 0;
end

% Use canonical mesh or generate new meshes
%------------------------------------------
if D.inv{D.val}.mesh.canonical
    
    % create head structure from previous tesselation  - skull
    %----------------------------------------------------------------------
    Tess             = Mesh.tess_iskull;
    head.XYZmm       = Tess.vert';
    head.tri         = Tess.face';
    head.nr          = [length(Tess.vert) length(Tess.face)];
    head.Centre      = [0 -20 0];
    model.head       = head;

    % create head structure from previous tesselation  - scalp
    %----------------------------------------------------------------------
    Tess             = Mesh.tess_scalp;
    head.XYZmm       = Tess.vert';
    head.tri         = Tess.face';
    head.nr          = [length(Tess.vert) length(Tess.face)];
    model.head(2)    = head;

else

% May involve regenerating an existing individual mesh, but hey!
    Pvol             = strvcat(Mesh.msk_iskull,Mesh.msk_scalp);
    model.head       = spm_eeg_inv_model('GenMesh',Pvol);
    
end

% Place the electrodes in the model
%--------------------------------------------------------------------------
model.electrodes     = spm_eeg_inv_model('Elec2Scalp',model.head(2),Reg',el_name);

% set the sphere model
%--------------------------------------------------------------------------
model.spheres        = spm_eeg_inv_Rsph(model,[]);

% Save in D
%--------------------------------------------------------------------------
D.inv{D.val}.forward = model;
%fprintf('Foward model complete - thank you\n')

return





% This code will constuct a forward model based on a standard head model & 
% electrode placement.
%==========================================================================

% get standard head model
%--------------------------------------------------------------------------
load(fullfile(spm('dir'),'EEGcanonical','standard_head_model.mat'))

% Input standard electrode sets
%--------------------------------------------------------------------------
[set_Nel,set_name]   = spm_eeg_inv_electrset;
el_set               = spm_input('Which set of electrodes','+1','m',set_name);
[el_sphc,el_name]    = spm_eeg_inv_electrset(el_set);
flags_el.q_RealLoc   = 0;
flags_el.q_RealNI    = 0;
[electrodes,f_el]    = spm_eeg_inv_model('Elec2Scalp',model.head(end), ...
                       el_sphc,el_name,flags_el);
model.electrodes     = electrodes;
model.flags.fl_elec  = f_el;


% [re]-set the sphere model
%--------------------------------------------------------------------------
model.spheres        = spm_eeg_inv_Rsph(model,[]);
D.inv{D.val}.forward = model;





