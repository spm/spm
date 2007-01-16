function D = spm_eeg_inv_mesh_ui(varargin)

%=======================================================================
% Cortical Mesh user-interface routine
% Invokes spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(D,val)
% Input:
% D		   - input data struct (optional)
% Output:
% D	       - same data struct including the meshing files and variables
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 716 2007-01-16 21:13:50Z karl $


% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});


% sMRI spatial normalization into MNI T1 template
%--------------------------------------------------------------------------
D.inv{val}.mesh.sMRI = spm_select(1,'image','Select subject''s structural MRI');
D                    = spm_eeg_inv_spatnorm(D);

% get cortical mesh size and compute meshes
%--------------------------------------------------------------------------
D.inv{val}.mesh.Msize = spm_input('Mesh size (vertices)','+1','3000|4000|5000|7200',[1 2 3 4]);
D                     = spm_eeg_inv_meshing(D);

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);
