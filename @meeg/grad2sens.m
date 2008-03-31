function this = grad2sens(this, grad)
% Method for converting configuring the parameters of MEG
% sensors and channels based on the grad struct provided by fileio
% FORMAT this = grad2sens(this, grad)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: grad2sens.m 1280 2008-03-31 10:59:39Z vladimir $

[sel1, sel2] = spm_match_str(chanlabels(this), grad.label);

% Check if there are any MEG channels in the data
if isempty(sel1)
    error('The grad struct and meeg object are not compatible');
end

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

% Remove existing MEG sensors
prevMEGind = strmatch(this.sensors.type, 'MEG', 'exact');
this.sensors.pnt(prevMEGind, :) = [];
this.sensors.ori(prevMEGind, :) = [];
this.sensors.type(prevMEGind) = [];
this.sensors.label(prevMEGind) = [];

% Adding MEG to sensors (MEG goes first)
this.sensors.pnt = [grad.pnt; this.sensors.pnt];
this.sensors.ori = [grad.ori; this.sensors.ori];
this.sensors.type = [repmat({'MEG'}, 1, size(this.sensors.pnt, 1)), this.sensors.type];

% For MEG the sensor labels are not really important. At the moment grad
% struct does not contain sensor labels so we'll use arbitrary labels
meglabel={};
for i = 1:size(this.sensors.pnt, 1)
    meglabel = [meglabel {['MEG_sensor_' num2str(i)]}];
end

this.sensors.label = [meglabel this.sensors.label];


% Set the 'tra' vectors for all MEG channels
this = getset(this, 'channels', 'tra', sel1, ...
    mat2cell(grad.tra(sel2, :), ones(1, length(sel2)), size(grad.tra, 2)));


% Set the type of channels found in grad to MEG
this = chantype(this, sel1, 'MEG');

% Check if there any other channels marked as MEG and change theor type.
nonMEGind = setdiff(1:nchannels(this), sel1);
if ~isempty(nonMEGind)
    illegalMEGchan = strmatch('MEG', chantype(this, nonMEGind));
    if ~isempty(illegalMEGchan)
        labels = channlabels(this, nonMEGind(illegalMEGchan));
        if iscell(labels)
            chanstring = '';
            for i = 1:length(labels)
                chanstring = [chanstring labels{i} ' '];
            end
        else
            chanstring = [labels ' '];
        end

        warning(['Channels ' chanstring ' do not appear in grad. Changing their type to ''Other''']);

        this = chantype(this, nonMEGind(illegalMEGchan), 'Other');

        this = getset(this, 'channels', 'tra', nonMEGind(illegalMEGchan), []);
    end
end


