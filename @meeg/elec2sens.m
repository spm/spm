function this = elec2sens(this, elec)
% Method for converting configuring the parameters of MEG
% sensors and channels based on the elec struct 
% FORMAT this = elec2sens(this, elec)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: elec2sens.m 1291 2008-04-02 13:58:28Z vladimir $

[sel1, sel2] = spm_match_str(chanlabels(this), elec.label);

% Initialize sensor fields if necessary
if ~isfield(this, 'sensors')
    this.sensors = [];
end

if ~isfield(this.sensors, 'pnt')
    this.sensors.pnt = [];
end

if ~isfield(this.sensors, 'ori')
    this.sensors.ori = [];
end

if ~isfield(this.sensors, 'type')
    this.sensors.type = {};
end

if ~isfield(this.sensors, 'label')
    this.sensors.label = {};
end

% Remove existing EEG sensors
prevEEGind = strmatch(this.sensors.type, 'EEG', 'exact');
this.sensors.pnt(prevEEGind, :) = [];
this.sensors.ori(prevEEGind, :) = [];
this.sensors.type(prevEEGind) = [];
this.sensors.label(prevEEGind) = [];

% Adding EEG to sensors (MEG goes first)
this.sensors.pnt = [this.sensors.pnt; elec.pnt];
this.sensors.ori = [this.sensors.ori; zeros(size(elec.pnt))];
this.sensors.type = [this.sensors.type repmat({'EEG'}, 1, size(elec.pnt, 1))];
this.sensors.label = [this.sensors.label(:)' elec.label(:)'];

if isfield(elec, 'fid');
    this.fiducials = elec.fid;
else
    error('Cannot process sensors without fiducials');
end
