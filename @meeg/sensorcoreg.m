function this = sensorcoreg(this)
% Coregisters sensor coordinates to standard template
% FORMAT  this = sensorcoreg(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorcoreg.m 1390 2008-04-14 16:08:09Z vladimir $

if ~isfield(this.sensors, 'pnt') ~isempty(this.sensors.pnt)
    error('No sensor positions are defined');
end

% if ~isempty(strmatch('MEG', this.sensors.type))
%     error('MEG not supported yet');
% end

if isfield(this.other, 'headshape')
    headshape = this.other.headshape;
elseif ~isfield(this.sensors, 'tra')
    headshape = this.sensors.pnt;
else
    headshape = this.fiducials.pnt;
end

if isempty(headshape)
    headshape = sparse(0, 3);
end

[datareg, mesh] = spm_eeg_inv_template(1);

[M1, sensors, fid_eeg, headshape] = ...
    spm_eeg_inv_datareg(this.sensors, this.fiducials.fid.pnt, datareg.fid_mri, headshape, datareg.scalpvert, 1);

datareg.sens_coreg = sensors.pnt;
datareg.fid_coreg = fid_eeg;
datareg.hsp_coreg = headshape;
datareg.label = sensors.label;

spm_eeg_inv_checkdatareg(mesh, datareg);

this.sensors = sensors;
this.fiducials = fid_eeg;