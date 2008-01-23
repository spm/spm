function data = preprocSPM(cfg)

% SPM_EEG_PREPROCESSING
%
% step 1: read header and do channel selection using a GUI
% step 2: read events and do trigger selection and trial definition using a GUI
% step 3: for each trial, read data and filter it

% these options relate to the actual preprocessing, it is neccessary to specify here because of channel selection
if ~isfield(cfg, 'reref'),        cfg.reref = 'no';             end
if ~isfield(cfg, 'refchannel'),   cfg.refchannel = {};          end
if ~isfield(cfg, 'implicitref'),  cfg.implicitref = [];         end
if ~isfield(cfg, 'detrend'),      cfg.detrend = 'no';           end
% set the defaults for the signal processing options
if ~isfield(cfg, 'blc'),          cfg.blc = 'no';               end
if ~isfield(cfg, 'lnfilter'),     cfg.lnfilter = 'no';          end
if ~isfield(cfg, 'dftfilter'),    cfg.dftfilter = 'no';         end
if ~isfield(cfg, 'lpfilter'),     cfg.lpfilter = 'no';          end
if ~isfield(cfg, 'hpfilter'),     cfg.hpfilter = 'no';          end
if ~isfield(cfg, 'bpfilter'),     cfg.bpfilter = 'no';          end
if ~isfield(cfg, 'medianfilter'), cfg.medianfilter  = 'no';     end
if ~isfield(cfg, 'hilbert'),      cfg.hilbert = 'no';           end
if ~isfield(cfg, 'derivative'),   cfg.derivative = 'no';        end
if ~isfield(cfg, 'rectify'),      cfg.rectify = 'no';           end
if ~isfield(cfg, 'boxcar'),       cfg.boxcar = 'no';            end
if ~isfield(cfg, 'absdiff'),      cfg.absdiff = 'no';           end
if ~isfield(cfg, 'timeshift'),    cfg.timeshift = '0';           end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hdr = read_header(cfg.dataset);

[unique_label junk ind]=unique(hdr.label);
if length(unique_label)~=length(hdr.label)
    warning(['Data file contains several channels with ',...
        'the same name. These channels cannot be processed and will be disregarded']);
    % This finds the repeating labels and removes all their occurences
    sortind=sort(ind);
    [junk ind2]=setdiff(hdr.label, unique_label(sortind(find(diff(sortind)==0))));
    hdr.label=hdr.label(ind2);
end

cfg.channel = ft_channelselection(cfg.channel,    hdr.label);

[junk chansel] = match_str(cfg.channel, hdr.label);
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEP 2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % here I would keep open the option for allowing people to use their own trialfuns
% fprintf('determining interesting segments in the data\n');
% if ~isfield(cfg, 'trialfun')
%   [trl, event] = trialfun_general(cfg)
% else
%   [trl, event] = feval(cfg.trialfun, cfg);
% end

% determine the length in samples to which the data should be padded before filtering is applied
if cfg.padding>0
    if strcmp(cfg.lnfilter, 'yes') || ...
            strcmp(cfg.dftfilter, 'yes') || ...
            strcmp(cfg.lpfilter, 'yes') || ...
            strcmp(cfg.hpfilter, 'yes') || ...
            strcmp(cfg.bpfilter, 'yes') || ...
            strcmp(cfg.medianfilter, 'yes')
        padding = round(cfg.padding * hdr.Fs);
    else
        % no filtering will be done, hence no padding is neccessary
        padding = 0;
    end
    % update the configuration (in seconds) for external reference
    cfg.padding = padding / hdr.Fs;
else
    % no padding was requested
    padding = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spm_progress_bar('Init', size(cfg.trl,1), 'reading and preprocessing'); drawnow;
for i=1:size(cfg.trl,1)
    spm_progress_bar('Set', i); drawnow;
    % non-zero padding is used for filtering and line noise removal
    nsamples = cfg.trl(i,2)-cfg.trl(i,1)+1;
    if nsamples>padding
        % the trial is already longer than the total lenght requested
        begsample  = cfg.trl(i,1);
        endsample  = cfg.trl(i,2);
        begpadding = 0;
        endpadding = 0;
    else
        % begpadding+nsamples+endpadding = total length of raw data that will be read
        begpadding = ceil((padding-nsamples)/2);
        endpadding = floor((padding-nsamples)/2);
        begsample  = cfg.trl(i,1) - begpadding;
        endsample  = cfg.trl(i,2) + endpadding;
        if begsample<1
            warning('cannot apply enough padding at begin of file');
            begpadding = begpadding - (1 - begsample);
            begsample  = 1;
        end
        if endsample>(hdr.nSamples*hdr.nTrials)
            warning('cannot apply enough padding at end of file');
            endpadding = endpadding - (endsample - hdr.nSamples*hdr.nTrials);
            endsample  = hdr.nSamples*hdr.nTrials;
        end
    end

    % read the raw data with padding on both sides of the trial
    dat = read_data(cfg.dataset, hdr, begsample, endsample, chansel, cfg.forcecont);

    nsamples = size(dat,2);
    time = (cfg.trl(i,3) - begpadding + (0:(nsamples-1)))/hdr.Fs;


    % FIXME rereferencing
    % do the rereferencing in case of EEG
    %   if ~isempty(cfg.implicitref) && ~any(strmatch(cfg.implicitref,label))
    %     label = {label{:} cfg.implicitref};
    %     dat(end+1,:) = 0;
    %   end
    %   if strcmp(cfg.reref, 'yes'),
    %     cfg.refchannel = channelselection(cfg.refchannel, label);
    %     refindx = match_str(label, cfg.refchannel);
    %     if isempty(refindx)
    %       error('reference channel was not found')
    %     end
    %     dat = avgref(dat, refindx);
    %   end

    % do the filtering
    if strcmp(cfg.medianfilter, 'yes'), dat = medfilt1(dat, cfg.medianfiltord, [], 2); end
    if strcmp(cfg.lpfilter, 'yes'),     dat = ft_lowpassfilter (dat, hdr.Fs, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir); end
    if strcmp(cfg.hpfilter, 'yes'),     dat = ft_highpassfilter(dat, hdr.Fs, cfg.hpfreq, cfg.hpfiltord, cfg.hpfilttype, cfg.hpfiltdir); end
    if strcmp(cfg.bpfilter, 'yes'),     dat = ft_bandpassfilter(dat, hdr.Fs, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir); end
    if strcmp(cfg.lnfilter, 'yes'),     dat = ft_notchfilter   (dat, hdr.Fs, cfg.lnfreq, cfg.lnfiltord); end
    if strcmp(cfg.dftfilter, 'yes'),
        for i=1:length(cfg.dftfreq)
            % filter out the 50Hz noise, optionally also the 100 and 150 Hz harmonics
            dat = ft_dftfilter(dat, hdr.Fs, cfg.dftfreq(i));
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove the filter padding and do the preprocessing on the remaining trial data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if begpadding~=0 || endpadding~=0
        dat = dat(:, (1+begpadding):(end-endpadding));
        time = time((1+begpadding):(end-endpadding));
    end


    time = time - cfg.timeshift;

    % do the rest of the preprocessing
    if strcmp(cfg.detrend, 'yes')
        dat = detrend(dat')';
    end
    if strcmp(cfg.blc, 'yes')
        if isstr(cfg.blcwindow) && strcmp(cfg.blcwindow, 'all')
            blcbegsample = 1;
            blcendsample = nsamples;
            dat          = ft_blc(dat, blcbegsample, blcendsample);
        else
            % determine the begin and endsample of the baseline period and baseline correct for it
            blcbegsample = nearest(time, cfg.blcwindow(1));
            blcendsample = nearest(time, cfg.blcwindow(2));
            dat          = ft_blc(dat, blcbegsample, blcendsample);
        end
    end
    if strcmp(cfg.hilbert, 'yes'),
        dat = abs(hilbert(dat'))';
    end
    if strcmp(cfg.rectify, 'yes'),
        dat = abs(dat);
    end
    if isnumeric(cfg.boxcar)
        numsmp = round(cfg.boxcar*hdr.Fs);
        if ~rem(numsmp,2)
            % the kernel should have an odd number of samples
            numsmp = numsmp+1;
        end
        kernel = ones(1,numsmp) ./ numsmp;
        dat = convn(dat, kernel, 'same');
    end
    if strcmp(cfg.derivative, 'yes'),
        dat = [diff(dat, 1, 2) 0];
    end
    if strcmp(cfg.absdiff, 'yes'),
        % this implements abs(diff(data), which is required for jump detection
        dat = abs([diff(dat, 1, 2) 0]);
    end

    if strcmp(cfg.precision, 'single')
        % convert to single precision, i.e. 4 byte floating point values
        dat = single(dat);
    end

    data.time{i} = time;
    data.trial{i} = dat;

end % for all trials
spm_progress_bar('Clear');

data.label = cfg.channel;
data.hdr = hdr;
data.cfg = cfg;
data.fsample=hdr.Fs;




