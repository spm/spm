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
% $Id: spm_eeg_inv_mesh_ui.m 1006 2007-11-20 19:50:53Z karl $


% initialise
%--------------------------------------------------------------------------
[D,val]         = spm_eeg_inv_check(varargin{:});
D.inv{val}.mesh = [];

% get cortical mesh size and compute meshes
%--------------------------------------------------------------------------
D               = spm_eeg_inv_meshing(D);

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);
