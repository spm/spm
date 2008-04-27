function this = sensorcoreg(this)
% Coregisters sensor coordinates to standard template
% FORMAT  this = sensorcoreg(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorcoreg.m 1488 2008-04-27 14:11:48Z vladimir $

[ok, this] = checkmeeg(struct(this), 'sensfid');

if ~ok
    error('Coregistration cannot be performed due to missing data');
end

this = meeg(this);

[vol,mrifid, mesh] = spm_eeg_inv_template(1);

senstypes = {'EEG', 'MEG'};

for i = 1:numel(senstypes)
    if ~isempty(sensors(this, senstypes{i}))
        S =[];
        S.sens = sensors(this, senstypes{i});
        S.meegfid = fiducials(this);
        S.vol = vol;
        S.mrifid = mrifid;
        S.template = 1;
        switch senstypes{i}
            case 'MEG'
                error('MEG coregistration is under construction');
            case 'EEG'
                [M1, sens, fid] = spm_eeg_inv_datareg_eeg(S);
        end
        this = sensors(this, senstypes{i}, sens);
    end
end

S.sens = sens;
S.meegfid = fid;
S.mesh = mesh;

%try
    spm_eeg_inv_checkdatareg(S);
%end

this = fiducials(this, fid);