function D = spm_eeg_inv_mesh_ui(varargin)
% Cortical Mesh user-interface routine
% Invokes spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(D,val)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the meshing files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 1437 2008-04-17 10:34:39Z christophe $


% initialise
%--------------------------------------------------------------------------
[D,ival]         = spm_eeg_inv_check(varargin{:});
D.inv{ival}.mesh = [];

% get cortical mesh size and compute meshes
%--------------------------------------------------------------------------
D               = spm_eeg_inv_meshing(D);

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);
