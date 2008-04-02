function this = sensorcoreg(this)
% Coregisters sensor coordinates to standard template
% FORMAT  this = sensorcoreg(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorcoreg.m 1291 2008-04-02 13:58:28Z vladimir $

if ~isfield(this.sensors, 'pnt') ~isempty(this.sensors.pnt)
    error('No sensor positions are defined');
end

if ~isempty(strmatch('MEG', this.sensors.type))
    error('MEG not supported yet');
end

val=1;
this.other(1).inv{val}.mesh.Msize = 1;
this = spm_eeg_inv_template(this, val);

allpoints = [this.sensors.pnt; this.fiducials];

if isfield(this.other, 'headshape')
    headshape = this.other.headshape;
    allpoints = [allpoints; headshape];
else
    eegind = strmatch('EEG', this.sensors.type, 'exact');
    headshape = this.sensors.pnt(eegind, :);
end


% The coregistration function doesn't like coordinates smaller than 1
% So the units are normalized to have values in the order of magnitude of 10.
norm_const = 0.1*mean(std(allpoints));

this.other.inv{val}.datareg.sensors = this.sensors.pnt./norm_const;
this.other.inv{val}.datareg.label = this.sensors.label;
this.other.inv{val}.datareg.fid_eeg = this.fiducials./norm_const;
this.other.inv{val}.datareg.headshape = headshape./norm_const;

this.other.modality = 'EEG'; % for now


[RT,sensors_reg,fid_reg,headshape_reg,orient_reg] = spm_eeg_inv_datareg(this);

this.other.inv{val}.datareg.sens_coreg=sensors_reg;
this.other.inv{val}.datareg.fid_coreg=fid_reg;
this.other.inv{val}.datareg.hsp_coreg=headshape_reg;

spm_eeg_inv_checkdatareg(this);

% 
% % remove effects of inv_datareg
% if exist('invbackup')==1
%     D.inv=invbackup;
% else
%     D=rmfield(D, 'inv');
% end

