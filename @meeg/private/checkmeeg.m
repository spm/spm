function [result meegstruct]=checkmeeg(meegstruct, option)
% Function for checking the internal struct of meeg objects
% FORMAT [result meegstruct]=checkmeeg(meegstruct, option)
% result - 1 - OK, 0- failed
% meegstruct - the struct to check (is returned modified if necessary)
% option - 'basic' (default) - just check the essential fields
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: checkmeeg.m 1281 2008-03-31 17:05:49Z vladimir $

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
end

try 
    meegstruct.transform.ID;
catch
    meegstruct.transform.ID = 'time';
    disp('checkmeeg: transform type missing, assigning default');
end

if strcmp(meegstruct.transform.ID, 'TF')
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
        meegstruct.data.datatype = 'float32';
    end

    if ~isfield(meegstruct.data, 'scale')
        if strcmp(meegstruct.data.datatype, 'float32') || ...
                strcmp(meegstruct.data.datatype, 'float64')
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
                    [Nchannels Nsamples Ntrials], spm_type(meegstruct.data.datatype));

            case 'TF'
                meegstruct.data.y = file_array(fullfile(filepath, meegstruct.data.fnamedat), ...
                    [Nchannels Nfrequencies Nsamples Ntrials], spm_type(meegstruct.data.datatype));
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
        case 'TF'
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
    elseif Ntrials==1
        meegstruct.type = 'evoked';
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

if strcmp(option, 'plot2D')
    if (Nsamples == 0)
        result =0;
        return;
    else
        if ~all(isfield(meegstruct.channels, {'X_plot2D', 'Y_plot2D'}))
            [meegstruct.channels(:).X_plot2D] = NaN;
            [meegstruct.channels(:).Y_plot2D] = NaN;
        end
        for i = 1:Nchannels
            X = meegstruct.channels(i).X_plot2D;
            Y = meegstruct.channels(i).Y_plot2D;
            if ~(~isempty(X) && ~isempty(Y) && ...
                    0<=X && X<=1 && 0<=Y && Y<=1)
                meegstruct.channels(i).X_plot2D = NaN;
                meegstruct.channels(i).X_plot2D = NaN;
            end
        end
        result = 1;
        return;
    end
end



