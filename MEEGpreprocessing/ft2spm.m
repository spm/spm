function ft2spm(ftdata, filename, ctf)
% Converter from Fieldtrip (http://www.ru.nl/fcdonders/fieldtrip/)
% data structures to SPM8 file format
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $


disp('Warning: Converting the data to SPM8 format which is in development and');
disp('may not work with some (or most) SPM functions. If you need something that works,');
disp('set your SPM version to SPM5');

% If preprocessing format
if iscell(ftdata.time)

    if length(ftdata.time)>1
        % Initial checks
        if any(diff(cellfun('length', ftdata.time))~=0)
            error('SPM can only handle data with equal trial lengths.');
        else
            times=cell2mat(ftdata.time(:));
            if any(diff(times(:, 1))~=0) || any(diff(times(:, end))~=0)
                error('SPM can only handle data the same trial borders.');
            end
        end
    end
    
    Ntrials=length(ftdata.trial);
    Nchannels  = size(ftdata.trial{1},1);
    Nsamples  = size(ftdata.trial{1},2);

    data = zeros(Nchannels, Nsamples, Ntrials);

    spm_progress_bar('Init', Ntrials, 'converting data to SPM format'); drawnow;
    for n=1:Ntrials
        spm_progress_bar('Set', n); drawnow;
        data(:,:,n) = ftdata.trial{n};
    end
    spm_progress_bar('Clear');
    
    ftdata.time  = ftdata.time{1};
else % timelockanalysis format
    rptind=strmatch('rpt', tokenize(ftdata.dimord, '_'));
    if isempty(rptind)
        rptind=strmatch('subj', tokenize(ftdata.dimord, '_'));
    end

    timeind=strmatch('time', tokenize(ftdata.dimord, '_'));
    chanind=strmatch('chan', tokenize(ftdata.dimord, '_'));

    if ~isempty(rptind)
        if isfield(ftdata, 'trial')
            Ntrials = size(ftdata.trial, rptind);
            data =permute(ftdata.trial, [chanind, timeind, rptind]);
        else
            Ntrials = size(ftdata.individual, rptind);
            data =permute(ftdata.individual, [chanind, timeind, rptind]);
        end
    else
        Ntrials = 1;
        data =permute(ftdata.avg, [chanind, timeind]);
    end
    
    Nchannels  = size(ftdata.trial, chanind);
    Nsamples  = size(ftdata.trial, timeind);
end

% sampling rate in Hz
D.Fsample = ftdata.fsample;

D.timeOnset = ftdata.time(1);

% Number of time bins in peri-stimulus time
D.Nsamples = Nsamples;

% Names of channels in order of the data
D.channels = struct('label', ftdata.label);
[D.channels.modality] = deal('Other');


% name of SPM mat file
[filepath, filename] = fileparts(filename);
D.fname = filename;
D.path = filepath;

megchans=match_str(ftdata.label, ft_channelselection('MEG', ftdata.label));
if isfield(ftdata, 'hdr') &&  isfield(ftdata.hdr, 'grad') && ~isempty(megchans) % MEG
    [D.channels(megchans).type] = deal('MEG');
    [D.channels(megchans).units] = deal('fT');
    % This is a rule of thumb. Look at the distribution of values
    % at one timepoint across channels and trials. Assumes that in MEG
    % datasets most channels will be MEG.
    if (median(abs(reshape(data(:,fix(end/2),:), 1,[])))<1e-10)
        data(megchans, :, :)=1e15*data(megchans, :, :);
        warning('Detected MEG data and converted to the units to fT');
    end
    
    if ~isempty(ftdata.hdr.grad)
        D.sensors.pnt = ftdata.hdr.grad.pnt;
        D.sensors.ori = ftdata.hdr.grad.ori;
        D.sensors.type = cell(size(D.sensors.pnt, 1), 1);
        [D.sensors.type{:}] = deal('MEG');
        [sel1, sel2] = match_str({D.channels.label}, ftdata.hdr.grad.label);        
        tra = num2cell(ftdata.hdr.grad.tra(sel2, :), 2);       
        [D.channels(sel1).tra] = deal(tra{:});
    end
else   %EEG
    [D.channels.type] = deal('EEG');
    [D.channels.units] = deal('\muV');
    if (median(abs(reshape(data(:,fix(end/2),:), 1,[])))<1e-2)
        data=1e6*data;
        warning('Scaled the data for CTF EEG.');
    else
        warning('Assuming that units are uV. Units may not be correct for EEG data.');
    end
end

% In the case of epoched data preprocessed in SPM codes and times are taken
% from the trl table
if issubfield(ftdata, 'cfg.trl') && Ntrials>1
    trl = getsubfield(ftdata, 'cfg.trl');
    D.trials = struct('onset', num2cell(trl(:,1)./D.Fsample));
    if size(trl, 2)<4
        [D.trials.label] = deal('Undefined');
    else
        for i = 1:Ntrials
           D.trials(i).label = num2str(trl(i,4));
        end
    end
    %If there is no trl (data is continuous or not from SPM GUI)
elseif issubfield(ftdata.cfg, 'cfg.event')
    % If data comes from FT spm codes are added to events
    event = getsubfield(ftdata, 'cfg.event');
    if ~isfield(event, 'spmcode')
        event=spm_eeg_recode_events(event);
    end
    for i = 1:Ntrials
        D.trials(i).label = num2str(event(i).spmcode);
        D.trials(i).onset = event(i).sample./D.Fsample + D.timeOnset;
        D.trials(i).events.type = event(i).type;
        D.trials(i).events.value = event(i).value;
        D.trials(i).events.time = event(i).sample./D.Fsample;
    end 
else
    % If there is no information about events the codes are 1s
    warning('No event information found in the FT data');
    [D.trials(1:Ntrials).label]=deal('Undefined');
end

D.data.fnamedat = [spm_str_manip(D.fname, 'r') '.dat'];

D.data.scale = ones(Nchannels, 1, Ntrials);
D.data.datatype  = 'float32';

fpd = fopen(fullfile(D.path, D.data.fnamedat), 'w');

for i = 1:Ntrials
    d = squeeze(data(:, :, i));
    D.data.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.data.datatype);
end

fclose(fpd);

D = meeg(D);

save(fullfile(filepath, filename), 'D');

