function [this] = sensorcoreg(this)
% Coregisters sensor coordinates to standard template
% FORMAT  this = sensorcoreg(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorcoreg.m 1644 2008-05-14 21:14:43Z christophe $

[ok, this] = checkmeeg(struct(this), 'sensfid');

if ~ok
    error('Coregistration cannot be performed due to missing data');
end

this = meeg(this);

if isempty(sensors(this, 'EEG'))
    error('No EEG sensors found');
end


[eegvol, megvol, mrifid, mesh] = spm_eeg_inv_template(1);

S =[];
S.sens = sensors(this, 'EEG');
S.meegfid = fiducials(this);
S.vol = eegvol;
S.mrifid = mrifid;
S.template = 1;
M1 = spm_eeg_inv_datareg(S);

sens = forwinv_transform_sens(M1, S.sens);
fid = forwinv_transform_headshape(M1, S.meegfid);

this = sensors(this, 'EEG', sens);
this = fiducials(this, fid);
% this.M1 = M1;

S.sens = sens;
S.meegfid = fid;
S.mesh = mesh;
S.M1 = M1;

% Work around to pass the coreg parameters...
D = struct(this)
D.other.S = S;
this = meeg(D);

spm_eeg_inv_checkdatareg(S);