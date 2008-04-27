function D = spm_eeg_inv_datareg_ui(varargin)
% Data registration user-interface routine
% commands the EEG/MEG data co-registration within original sMRI space
%
% FORMAT D = spm_eeg_inv_mesh_ui(D,[val])
% Input:
% S         - input data struct (optional)
% Output:
% D         - same data struct including the new required files and variables
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg_ui.m 1488 2008-04-27 14:11:48Z vladimir $

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    D.inv{val}.mesh.template;
catch
    D.inv{val}.mesh.template = 0;
end

iseeg = ~isempty(strmatch('EEG', D.chantype));
ismeg = ~isempty(strmatch('MEG', D.chantype));

if iseeg && ismeg
    modality = spm_input('Which modality?','+1','EEG|MEG');
elseif iseeg
    modality = 'EEG';
elseif ismeg
    modality = 'MEG';
else
    error('No MEEG channels in the data');
end

D.inv{val}.modality = modality;

channels = D.chanlabels;
chanind = strmatch(D.chantype, modality);
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
switch D.inv{val}.modality
    case 'EEG'
        [M1, sens, fid] = spm_eeg_inv_datareg_eeg(S);
    case 'MEG'
        error('MEG support is under construction');
end

D.inv{val}.datareg.sensors = sens;
D.inv{val}.datareg.fid_eeg   = fid;

S.sens = sens;
S.meegfid = fid;
S.mesh = D.inv{val}.mesh;

% check and display registration
%--------------------------------------------------------------------------
spm_eeg_inv_checkdatareg(S);
