function D = spm_eeg_inv_datareg_ui(varargin)
% Data registration user-interface routine
% commands the EEG/MEG data co-registration within original sMRI space
%
% FORMAT D = spm_eeg_inv_mesh_ui(D,[val], modality)
% Input:
% Output:
% D         - same data struct including the new required files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_inv_datareg_ui.m 1650 2008-05-15 10:22:31Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    D.inv{val}.mesh.template;
catch
    D.inv{val}.mesh.template = 0;
end

if nargin > 2
    modality = varargin{3};
else
    modality = spm_eeg_modality_ui(D, 1);
end

D.inv{val}.modality = modality;

channels = D.chanlabels;
chanind = strmatch(modality, D.chantype);
chanind = setdiff(chanind, D.badchannels);
D.inv{val}.forward.channels = channels(chanind);

% register
%==========================================================================
S =[];
S.sens = D.sensors(modality);
S.meegfid = D.fiducials;
S.vol = D.inv{val}.forward.vol;
S.mrifid = D.inv{val}.datareg.fid_mri;
S.template = D.inv{val}.mesh.template;
%--------------------------------------------------------------------------

M1 = spm_eeg_inv_datareg(S);

switch D.inv{val}.modality
    case 'EEG'
        D.inv{val}.datareg.sensors = forwinv_transform_sens(M1, S.sens);
        D.inv{val}.datareg.fid_eeg = forwinv_transform_headshape(M1, S.meegfid);
        D.inv{val}.datareg.fid_mri = S.mrifid;
        D.inv{val}.datareg.fromMNI = eye(4);
        D.inv{val}.datareg.toMNI = eye(4);
    case 'MEG'
        D.inv{val}.forward.vol = forwinv_transform_vol(inv(M1), S.vol);
        D.inv{val}.datareg.fid_mri = forwinv_transform_headshape(inv(M1), S.mrifid);
        D.inv{val}.mesh = spm_eeg_inv_transform_mesh(inv(M1), D.inv{val}.mesh);
        D.inv{val}.datareg.sensors = S.sens;
        D.inv{val}.datareg.fid_eeg = S.meegfid; 
        D.inv{val}.datareg.fromMNI = inv(M1);
        D.inv{val}.datareg.toMNI = M1;
end

%%
S =[];

switch D.inv{val}.modality
    case 'EEG'
        S.sens = D.inv{val}.datareg.sensors;
    case 'MEG'    
        cfg = [];
        cfg.style = '3d';
        cfg.rotate = 0;
        cfg.grad = D.inv{val}.datareg.sensors;
   
        lay = ft_prepare_layout(cfg);
                
        S.sens = [];
        S.sens.label = lay.label(:, 1);
        S.sens.pnt = lay.pos;
end
        
S.meegfid = D.inv{val}.datareg.fid_eeg;
S.vol = D.inv{val}.forward.vol;
S.mrifid = D.inv{val}.datareg.fid_mri;
S.mesh = D.inv{val}.mesh;

% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(S);
