function [D val] = spm_eeg_inv_template_ui(varargin)
% Put the template head model in the inv struct
% FORMAT D = spm_eeg_inv_template_ui(D, val, Msize)
%   Msize   - index for precalculated cortical mesh density
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_template_ui.m 1535 2008-05-01 17:08:22Z vladimir $

[D,val] = spm_eeg_inv_check(varargin{:});

if val == 0
    val = 1;
end

if nargin>2
    Msize = varargin{3};
else
    Msize = spm_input('Mesh size (vertices)', '+1','3000|4000|5000|7200', [1 2 3 4]);
end

[eegvol, megvol, fid, mesh] = spm_eeg_inv_template(Msize);

if ~isfield(D, 'inv')
    D.inv = {};
end

D.inv{val}.mesh = mesh;
D.inv{val}.datareg.fid_mri = fid;

if ~isfield(D.inv{val}, 'modality')
    modality = spm_eeg_modality_ui(D, 1);
    D.inv{val}.modality = modality;
end

switch D.inv{val}.modality
    case 'EEG'
        D.inv{val}.forward.vol = eegvol;
    case 'MEG'
        D.inv{val}.forward.vol = megvol;
end
