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
% $Id: spm_eeg_inv_mesh_ui.m 1488 2008-04-27 14:11:48Z vladimir $


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
