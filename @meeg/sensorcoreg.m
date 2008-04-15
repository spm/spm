function this = sensorcoreg(this)
% Coregisters sensor coordinates to standard template
% FORMAT  this = sensorcoreg(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorcoreg.m 1406 2008-04-15 09:37:59Z vladimir $

[ok, this] = checkmeeg(struct(this), 'sensfid');

if ~ok
    error('Coregistration cannot be performed due to missing data');
end

this = meeg(this);

[datareg, mesh] = spm_eeg_inv_template(1);

senstypes = {'EEG', 'MEG'};

for i = 1:numel(senstypes)
    if ~isempty(sensors(this, senstypes{i}))
        [M1, sens, fid] = spm_eeg_inv_datareg(sensors(this, senstypes{i}), fiducials(this), datareg, 1);

        this = sensors(this, senstypes{i}, sens);
    end
end

try
    spm_eeg_inv_checkdatareg(mesh, datareg, sens, fid);
end

this = fiducials(this, fid);