function D = spm_eeg_inv_mesh_ui(varargin)
% Cortical Mesh user-interface routine
% Invokes spatial normalization (if required) and the computation of
% the proper size individual size
%
% FORMAT D = spm_eeg_inv_mesh_ui(D, val, template, Msize)
% Input:
% D        - input data struct (optional)
% Output:
% D        - same data struct including the meshing files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_mesh_ui.m 2452 2008-11-10 18:45:32Z vladimir $


% initialise
%--------------------------------------------------------------------------
[Finter] = spm('FnUIsetup','Define head model',0);

[D,val] = spm_eeg_inv_check(varargin{:});

if val == 0
    val = 1;
end

if ~isfield(D, 'inv')
    D.inv = {struct([])};
end

if ~isfield(D.inv{val}, 'modality')
    modality = spm_eeg_modality_ui(D, 1);
    D.inv{val}(1).modality = modality;
end

if nargin>2
    template = varargin{3};
else
    template = [];
end

if isempty(template)
    template = spm_input('Select head  model', '+1','template|individual', [1 0]);
end

if nargin>3
    Msize = varargin{4};
else
    Msize = spm_input('Mesh size (vertices)', '+1','3000|4000|5000|7200', [1 2 3 4]);
end

if template
    [vol, fid, mesh] = spm_eeg_inv_template(Msize, D.inv{val}.modality);
else
    if strcmp(D.inv{val}(1).modality, 'MEG')
        sMRI =  spm_select(1,'image', 'Select the subject''s structural image');
        [vol, fid, mesh] = spm_eeg_inv_meshing(sMRI, Msize, D.inv{val}.modality);
    else
        warndlg('MRI-based head models are not supported for EEG')
        return;
    end
end

D.inv{val}.mesh = mesh;
D.inv{val}.datareg.fid_mri = fid;
D.inv{val}.forward.vol = vol;

% check meshes and display
%--------------------------------------------------------------------------
spm_eeg_inv_checkmeshes(D);
