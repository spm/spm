function [result meegstruct]=checkmeeg(meegstruct, option)
% Function for checking the internal struct of meeg objects
% FORMAT [result meegstruct]=checkmeeg(meegstruct, option)
% result - 1 - OK, 0- failed
% meegstruct - the struct to check (is returned modified if necessary)
% option - 'basic' (default) - just check the essential fields
%          'sensfid' - also checks sensor and fiducial definitions
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: checkmeeg.m 2320 2008-10-09 10:01:22Z vladimir $

if nargin==1
    option = 'basic';
end

result=0;

if ~isfield(meegstruct, 'Nsamples')
    disp('checkmeeg: number of samples per trial is missing');
    return;
else
    Nsamples = meegstruct.Nsamples;
end

if ~isfield(meegstruct, 'Fsample') && (Nsamples~=0)
    disp('checkmeeg: sampling rate is missing');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.Fsample = 0;
end

if ~isfield(meegstruct, 'timeOnset')
    meegstruct.timeOnset = 0;
end

if ~isfield(meegstruct, 'trials') && (Nsamples~=0)
    disp('checkmeeg: no trials description');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.trials = struct([]);
else
    Ntrials = length(meegstruct.trials);
    if ~isfield(meegstruct.trials, 'label')
        disp('checkmeeg: no trial label, assigning default');
        [meegstruct.trials.label] = deal('Undefined');
    end
    if ~isfield(meegstruct.trials, 'bad')
        [meegstruct.trials.bad] = deal(0);
    end
    if ~isfield(meegstruct.trials, 'events')
        [meegstruct.trials.events] = deal([]);
    end
    if ~isfield(meegstruct.trials, 'onset')
        [meegstruct.trials.onset] = deal([]);
    end
end

if ~isfield(meegstruct, 'channels') && (Nsamples~=0)
    disp('checkmeeg: no channels description');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.channels = struct([]);
else
    Nchannels = length(meegstruct.channels);
    if ~isfield(meegstruct.channels, 'label')
        disp('checkmeeg: no channel label, assigning default');
        for i = 1:Nchannels
            meegstruct.channels(i).label = ['Ch' num2str(i)];
        end
    end
    if ~isfield(meegstruct.channels, 'bad')
        [meegstruct.channels.bad] = deal(0);
    end
    if ~isfield(meegstruct.channels, 'type')
        disp('checkmeeg: no channel type, assigning default');
        [meegstruct.channels.type] = deal('Other');
    end
    if ~isfield(meegstruct.channels, 'X_plot2D')
        [meegstruct.channels.X_plot2D] = deal([]);
        [meegstruct.channels.Y_plot2D] = deal([]);
    end
    if ~isfield(meegstruct.channels, 'units')
        disp('checkmeeg: no units, assigning default');
        [meegstruct.channels.units] = deal('unknown');
    else
        [meegstruct.channels(find(cellfun('isempty', {meegstruct.channels.units}))).units] = deal('unknown');
    end
end

try
    meegstruct.transform.ID;
catch
    meegstruct.transform.ID = 'time';
    disp('checkmeeg: transform type missing, assigning default');
end

if strncmp(meegstruct.transform.ID, 'TF', 2) % TF and TFphase
    try
        Nfrequencies = length(meegstruct.transform.frequencies);
    catch
        error('Information about frequencies missing');
    end
end

if ~isfield(meegstruct, 'data') && (Nsamples~=0)
    disp('checkmeeg: no data field');
    return;
elseif (Nsamples==0)  % This is to enable creation of empty meeg objects
    meegstruct.data = struct([]);
else
    if ~isfield(meegstruct.data, 'fnamedat')
        disp('checkmeeg: data file name missing');
        return;
    end
    if ~isfield(meegstruct.data, 'datatype')
        disp('checkmeeg: data type missing, assigning default');
        meegstruct.data.datatype = 'float32-le';
    end

    if ~isfield(meegstruct.data, 'scale')
        if strcmp(meegstruct.data.datatype, 'float32-le') || ...
                strcmp(meegstruct.data.datatype, 'float64-le')
            disp('checkmeeg: data scale missing, assigning default');
            meegstruct.data.scale = ones(Nchannels, 1, Ntrials);
        else
            disp('checkmeeg: data scale missing');
            return
        end
    end

    if ~isfield(meegstruct.data, 'y')
        meegstruct.data.y=[];
    end

    if isa(meegstruct.data.y, 'file_array')
        try
            % Try reading data
            meegstruct.data.y(1, 1, 1);
        catch
            meegstruct.data.y = [];
        end
    end

    if ~isa(meegstruct.data.y, 'file_array')
        if isfield(meegstruct, 'path')
            filepath = meegstruct.path;
        else
            filepath = '';
        end
        switch(meegstruct.transform.ID)
            % note: scale no longer used, must insure data is in some float
            % format
            case 'time'
                meegstruct.data.y = file_array(fullfile(filepath, meegstruct.data.fnamedat), ...
                    [Nchannels Nsamples Ntrials], meegstruct.data.datatype);

            case 'TF'
                meegstruct.data.y = file_array(fullfile(filepath, meegstruct.data.fnamedat), ...
                    [Nchannels Nfrequencies Nsamples Ntrials], meegstruct.data.datatype);
                if Ntrials>1
                    expected_size = [Nchannels Nfrequencies Nsamples Ntrials];
                else
                    expected_size = [Nchannels Nfrequencies Nsamples];
                end

            otherwise
                error('Unknown transform type');
        end

    end

    switch(meegstruct.transform.ID)
        case 'time'
            if Ntrials>1
                expected_size = [Nchannels Nsamples Ntrials];
            else
                expected_size = [Nchannels Nsamples];
            end
        case {'TF', 'TFphase'}
            if Ntrials>1
                expected_size = [Nchannels Nfrequencies Nsamples Ntrials];
            else
                expected_size = [Nchannels Nfrequencies Nsamples];
            end

            
        otherwise
            error('Unknown transform type');
    end

    if any(size(meegstruct.data.y) ~= expected_size)
        disp('checkmeeg: data size does not match the header');
        return;
    end
end

if ~isfield(meegstruct, 'type')
    disp('checkmeeg: data type is missing, assigning default');
    % rule of thumb - 10 sec
    if Nsamples == 0
        meegstruct.type = 'continuous';
    elseif Ntrials==1 && (Nsamples/meegstruct.Fsample) > 10
        meegstruct.type = 'continuous';
    elseif numel(unique({meegstruct.trials.label})) == Ntrials
        meegstruct.type = 'evoked';
        if ~isfield(meegstruct.trials, 'repl')
            disp('checkmeeg: missing replications number for evoked data, assuming 1');
            [meegstruct.trials.repl] = deal(1);
        end
    else
        meegstruct.type = 'single';
    end
end

if ~isfield(meegstruct, 'fname')
    meegstruct.fname = '';
end

if ~isfield(meegstruct, 'path')
    meegstruct.path = '';
end

if ~isfield(meegstruct, 'sensors')
    meegstruct.sensors = struct([]);
end

if ~isfield(meegstruct, 'fiducials')
    meegstruct.fiducials = struct([]);
end

if ~isfield(meegstruct, 'artifacts')
    meegstruct.artifacts = struct([]);
end

if ~isfield(meegstruct, 'transform')
    meegstruct.transform = struct('ID', 'time');
end

if ~isfield(meegstruct, 'other')
    meegstruct.other = struct([]);
end

if ~isfield(meegstruct, 'history')
    meegstruct.history = struct([]);
end

if ~isfield(meegstruct, 'cache')
    meegstruct.cache = struct([]);
end

% This makes sure the order of the fields in the struct is right. This
% seems to matter to the class function.

fieldnames_order = {
    'type'
    'Nsamples'
    'Fsample'
    'timeOnset'
    'trials'
    'channels'
    'data'
    'fname'
    'path'
    'sensors'
    'fiducials'
    'artifacts'
    'transform'
    'other'
    'history'
    'cache'};

[sel1, sel2] = spm_match_str(fieldnames_order, fieldnames(meegstruct));
tempcell = struct2cell(meegstruct);
meegstruct = cell2struct(tempcell(sel2), fieldnames_order, 1);

if strcmp(option, 'basic')
    result = 1;
    return;
end

if strcmp(option, 'dcm')
    chantypes = getset(meegstruct, 'channels', 'type');   
    if ~iscell(chantypes)
        chantypes = {chantypes};
    end
    if ismember('LFP', chantypes)
        if ismember('EEG', chantypes) || ismember('MEG', chantypes)
            disp('checkmeeg: DCM does not presently support files with both LFP and scalp channels');
            return;
        else
            result = 1;
            return;
        end
    end
end

if strcmp(option, 'sensfid') || strcmp(option, 'dcm')
    if isempty(meegstruct.sensors)
        disp('checkmeeg: no sensor positions are defined');
        return;
    end

    chantypes = getset(meegstruct, 'channels', 'type');
    eegind = strmatch('EEG', chantypes, 'exact');
    megind = strmatch('MEG', chantypes, 'exact');

    if ~isempty(eegind)
        if ~isfield(meegstruct.sensors, 'eeg') || isempty(meegstruct.sensors.eeg)
            disp('checkmeeg: EEG channel locations are not specified');
            return;
        else
            if ~isempty(setdiff({meegstruct.channels(eegind).label}, meegstruct.sensors.eeg.label))
                disp('checkmeeg: not all EEG channel locations are specified');
                return;
            end
        end
    end

    if ~isempty(megind)
        if ~isfield(meegstruct.sensors, 'meg') || isempty(meegstruct.sensors.meg)
            disp('checkmeeg: MEG channel locations are not specified');
            return;
        else
            if ~isempty(setdiff({meegstruct.channels(megind).label}, meegstruct.sensors.meg.label))
                disp('checkmeeg: not all MEG channel locations are specified');
                return;
            end
        end
    end

    if isempty(meegstruct.fiducials)
        disp('checkmeeg: no fiducials are defined');
        return;
    end

    if ~isfield(meegstruct.fiducials, 'pnt') || isempty(meegstruct.fiducials.pnt)
        if ~isempty(eegind)
            % Copy EEG sensors to fiducials.
            meegstruct.fiducials.pnt = meegstruct.sensors.eeg.pnt;
        else
            meegstruct.fiducials.pnt = sparse(0, 3);
        end
    end

    if ~isfield(meegstruct.fiducials, 'fid') || ...
            ~all(isfield(meegstruct.fiducials.fid, {'pnt', 'label'})) ||...
            (length(meegstruct.fiducials.fid.label) ~= size(meegstruct.fiducials.fid.pnt, 1)) || ...
            length(meegstruct.fiducials.fid.label) < 3
        disp('checkmeeg: at least 3 fiducials with labels are required');
        return
    end

    nzlbl = {'fidnz', 'nz', 'nas'};
    lelbl = {'fidle', 'fidt9', 'lpa', 'lear', 'earl' 'le', 't9'};
    relbl = {'fidre', 'fidt10', 'rpa', 'rear', 'earr', 're', 't10'};

    [sel1, nzind] = spm_match_str(nzlbl, lower(meegstruct.fiducials.fid.label));
    if isempty(nzind)
        disp('checkmeeg: could not find the nasion fiducial');
    else
        nzind = nzind(1);
    end
    
    [sel1, leind] = spm_match_str(lelbl, lower(meegstruct.fiducials.fid.label));
    if isempty(leind)
        disp('checkmeeg: could not find the left fiducial');
    else
        leind = leind(1);
    end

    [sel1, reind] = spm_match_str(relbl, lower(meegstruct.fiducials.fid.label));
    if isempty(reind)
        disp('checkmeeg: could not find the right fiducial');
    else
        reind = reind(1);
    end
    
    restind = setdiff(1:length(meegstruct.fiducials.fid.label), [nzind, leind, reind]);

    meegstruct.fiducials.fid.label = meegstruct.fiducials.fid.label([nzind, leind, reind, restind]);
    meegstruct.fiducials.fid.pnt = meegstruct.fiducials.fid.pnt([nzind, leind, reind, restind], :);

    result = 1;
end



