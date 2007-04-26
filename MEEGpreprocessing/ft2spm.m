function ft2spm(ftdata, filename, ctf)
% Converter from Fieldtrip (http://www.ru.nl/fcdonders/fieldtrip/)
% data structures to SPM file format
%_______________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging


disp('Converting data to SPM format');


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
    
    D.Nevents=length(ftdata.trial);
    
    nchan  = size(ftdata.trial{1},1);
    ntime  = size(ftdata.trial{1},2);
    
    data = zeros(nchan, ntime, D.Nevents);
    
    spm_progress_bar('Init', D.Nevents, 'converting data to SPM format'); drawnow;
    for n=1:D.Nevents
        spm_progress_bar('Set', n); drawnow;
        data(:,:,n) = ftdata.trial{n};
    end
    spm_progress_bar('Clear');
    
    ftdata.time  = ftdata.time{1};
else % timelockanalysis format
    rptind=strmatch('rpt', tokenize(ftdata.dimord, '_'));
    if isempty(rptind)
        rptind=strmatch('sbj', tokenize(ftdata.dimord, '_'));
    end
    if ~isempty(rptind)
        D.Nevents = size(ftdata.trial, rptind);
    else
        D.Nevents=1;
    end
    
    timeind=strmatch('time', tokenize(ftdata.dimord, '_'));
    chanind=strmatch('chan', tokenize(ftdata.dimord, '_'));
    
    if isfield(ftdata, 'trial')
        data =permute(ftdata.trial, [chanind, timeind, rptind]);
    else
        data =permute(ftdata.avg, [chanind, timeind, rptind]);
    end
end

% sampling rate in Hz
D.Radc = ftdata.fsample;

% Number of time bins in peri-stimulus time
D.Nsamples = length(ftdata.time);

% channel template file
D.channels.ctf = ctf;

% Names of channels in order of the data
D.channels.name = ftdata.label(:)';
D.Nchannels = length(D.channels.name);

% index vector of bad channels
D.channels.Bad = [];

% If the data is continuous un-epoched (rather than epoched or ERP)
% The rule of thumb of total duration longer than 1 sec is used to
% distinguish between continuous and ERP.
if D.Nevents>1 || (ftdata.time(end)-ftdata.time(1))<1
    % nr of time bins before stimulus onset (this excludes the time bin at zero)
    D.events.start = length(find(ftdata.time<-eps));
    % nr of time bins after stimulus onset (this excludes the time bin at zero)
    D.events.stop = length(find(ftdata.time>eps));
end

% name of SPM mat file
D.fname = filename;



if isfield(ftdata, 'hdr') &&  isfield(ftdata.hdr, 'grad') % MEG
    D.modality = 'MEG'; 
    D.units = 'fT';
    % This is a rule of thumb. Look at the distribution of values
    % at one timepoint across channels and trials. Assumes that in MEG
    % datasets most channels will be MEG.
    if (median(reshape(data(:,fix(end/2),:), 1,[]))<1e-10)
        data=1e15*data;
        warning('Detected MEG data and converted to the units to fT');
    end
else   %EEG
    D.modality = 'EEG'; 
    D.units = '\muV';
    warning('Assuming that units are uV. Units may not be correct for EEG data.');
end


% Don't change these
D.channels.reference = 0;
D.channels.ref_name = 'NIL';

% In the case of epoched data preprocessed in SPM codes and times are taken
% from the trl table
if isfield(ftdata, 'cfg') && isfield(ftdata.cfg, 'trl') && D.Nevents>1
    D.events.time = ftdata.cfg.trl(:,1)-ftdata.cfg.trl(:,3);
    if size(ftdata.cfg.trl, 2)<4
        D.events.code = ones(size(D.events.time));
    else
        D.events.code = ftdata.cfg.trl(:,4)';
    end
    %If there is no trl (data is continuous or not from SPM GUI)
elseif isfield(ftdata, 'cfg') && isfield(ftdata.cfg, 'event')
    % If data comes from FT spm codes are added to events
    if ~isfield(ftdata.cfg.event, 'spmcode')
        ftdata.cfg.event=spm_eeg_recode_events(ftdata.cfg.event);
    end
    D.events.time = [ftdata.cfg.event.sample];
    D.events.code = [ftdata.cfg.event.spmcode];
else
    % If there is no information about events the codes are 1s
    warning('No event information found in the FT data');
    D.events.code=ones(1, size(data,3));
end

% FIXME Support for original events for epoched data can be added.

D.events.types = unique(D.events.code);
D.events.Ntypes = length(D.events.types);


% The if is here because spm_eeg_epochs checks whether there is a previous
% reject field and does not replace it. So for continuous data there should
% not be any.
if D.Nevents>1
    D.events.reject = zeros(1, D.Nevents);
end

D.events.repl=[];


if D.Nevents>1 && length(D.events.code)~=D.Nevents
    warning('Number of events in the events field is not consistent with the number of trials');
end


% Name for the data file
D.path = fileparts(D.fname);
D.fname = spm_str_manip(D.fname, 't');
D.fnamedat = [spm_str_manip(D.fname, 'r') '.dat'];


Csetup = load(fullfile(spm('dir'), 'EEGtemplates', D.channels.ctf));

% Map name of channels (EEG only) to channel order specified in channel
% template file.
D.channels.heog = find(strncmpi('heog', D.channels.name, 4));
if isempty(D.channels.heog)
    D.channels.heog = 0;
else
    D.channels.heog = D.channels.heog(1);
end

D.channels.veog = find(strncmpi('veog', D.channels.name, 4));
if isempty(D.channels.veog)
    D.channels.veog = 0;
else
    D.channels.veog = D.channels.veog(1);
end

if Csetup.Nchannels>0
    switch D.modality
        case {'EEG', 'MEG'}
            % FIXME This is not a generic way for all kinds of EEG data.
            D.channels.eeg = setdiff([1:D.Nchannels], [D.channels.veog D.channels.heog]);
            
            for i = D.channels.eeg
                
                index = find(strcmpi(D.channels.name{i}, Csetup.Cnames));
                
                if isempty(index)
                    warning(sprintf('No channel named %s found in channel template file.', D.channels.name{i}));
                else
                    % take only the first found channel descriptor
                    D.channels.order(i) = index(1);
                end
                
            end
    end
end

D.scale = ones(D.Nchannels, 1, D.Nevents);
D.datatype  = 'float32';

fpd = fopen(fullfile(D.path, D.fnamedat), 'w');

for i = 1:D.Nevents
    d = squeeze(data(:, :, i));
    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
end


fclose(fpd);

if str2num(version('-release'))>=14
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end

spm('Pointer', 'Arrow');

