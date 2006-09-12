function D = spm_eeg_inv_mesh_ui(D)

%=======================================================================
% Cortical Mesh user-interface routine
% Invokes spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D			- same data struct including the meshing files and variables
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 621 2006-09-12 17:22:42Z karl $

spm_defaults

if nargin == 0
    try
    D   = spm_select(1, '.mat', 'Select EEG/MEG mat file');
    D   = spm_eeg_ldata(D);
    catch
        error(sprintf('Trouble reading the data file\n'));
    end
end

try
    val = D.val;
catch
    val = length(D.inv);
end

% use mesh sie of 3000
%--------------------------------------------------------------------------
D.inv{val}.mesh.Msize = 1;

% sMRI spatial normalization into MNI T1 template
%------------------------------------------------------------------
D = spm_eeg_inv_spatnorm(D);

% compute meshes
%------------------------------------------------------------------
D = spm_eeg_inv_meshing(D);

% check meshes and display
%------------------------------------------------------------------
D = spm_eeg_inv_checkmeshes(D);

save(D.fname,'D');
