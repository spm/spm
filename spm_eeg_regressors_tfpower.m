function res = spm_eeg_regressors_tfpower(S)
% Generate regressors from power in TF dataset
% S                     - input structure
% fields of S:
%    S.D                - M/EEG object
%
%    Additional parameters can be defined specific for each plugin
% Output:
%  res -
%   If no input is provided the plugin returns a cfg branch for itself
%
%   If input is provided the plugin returns
%______________________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak

SVNrev = '$Rev: 6188 $';

if nargin == 0
    %--------------------------------------------------------------------------
    % TF power dataset
    %--------------------------------------------------------------------------
    Dtf        = cfg_files;
    Dtf.tag    = 'Dtf';
    Dtf.name   = 'TF power dataset name';
    Dtf.filter = 'mat';
    Dtf.num    = [1 1];
    Dtf.help   = {'Select the M/EEG mat file containing TF power data'};
    
    %--------------------------------------------------------------------------
    % timewin
    %--------------------------------------------------------------------------
    timewin         = cfg_entry;
    timewin.tag     = 'timewin';
    timewin.name    = 'Time window';
    timewin.help    = {['Start and stop of the time window [ms]. ' ...
        '(Used only for summarised epoched cases)']};
    timewin.strtype = 'r';
    timewin.num     = [1 2];
    timewin.val     = {[-Inf Inf]};
    
    %--------------------------------------------------------------------------
    % freqwin
    %--------------------------------------------------------------------------
    freqwin         = cfg_entry;
    freqwin.tag     = 'freqwin';
    freqwin.name    = 'Frequency window';
    freqwin.help    = {'Start and stop of the frequency window (Hz).'};
    freqwin.strtype = 'r';
    freqwin.num     = [1 2];
    freqwin.val     = {[-Inf Inf]};
    
    %--------------------------------------------------------------------------
    % average
    %--------------------------------------------------------------------------
    average = cfg_menu;
    average.tag = 'average';
    average.name = 'Average over frequency';
    average.labels = {'Yes', 'No'};
    average.val = {1};
    average.values = {1,0};
    average.help = {'Average over the frequency window to produce one regressor',...
        'or output each frequency as a separate regressor'};
    
    %--------------------------------------------------------------------------
    % standardize
    %--------------------------------------------------------------------------
    standardize = cfg_menu;
    standardize.tag = 'standardize';
    standardize.name = 'Standardize';
    standardize.labels = {'Yes', 'No'};
    standardize.val = {0};
    standardize.values = {1,0};
    standardize.help = {'Standardize (zscore) regressor values'};
    
    %--------------------------------------------------------------------------
    % regname
    %--------------------------------------------------------------------------
    regname         = cfg_entry;
    regname.tag     = 'regname';
    regname.name    = 'Regressor name';
    regname.help    = {'Specify the string to be used as regressor name.'};
    regname.strtype = 's';
    regname.num     = [1 Inf];
    regname.val     = {'TFpower'};
    
    tfpower = cfg_branch;
    tfpower.tag = 'tfpower';
    tfpower.name = 'Time-frequency power';
    tfpower.val = {Dtf, spm_cfg_eeg_channel_selector, timewin, freqwin, average, standardize, regname};
    
    res = tfpower;
    
    return
end

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','Time-frequency power regressors');

if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end
if ~isfield(S, 'freqwin'),      S.freqwin  = [-Inf Inf];    end

if iscell(S.Dtf)
    S.Dtf = char(S.Dtf);
end

Dtf  = spm_eeg_load(S.Dtf);
D    = spm_eeg_load(S.D);

if ~isequal(Dtf.transformtype, 'TF');
    error('Time-frequency power dataset is expected as input.')
end

freqind = Dtf.indfrequency(min(S.freqwin)):Dtf.indfrequency(max(S.freqwin)); %BW
if isempty(freqind) || any(isnan(freqind))
    error('Selected frequency window is invalid.');
end

chanind = setdiff(Dtf.selectchannels(spm_cfg_eeg_channel_selector(S.channels)), Dtf.badchannels);

if isempty(chanind)
    error('No channels were selected');
end


if isequal(D.type, 'continuous')
    
    if ~isequal(Dtf.type, 'continuous') || (D.time(1) < Dtf.time(1)) || (D.time(end)>Dtf.time(end))
        error('All times of the input dataset should be within the power dataset.');
    end
    
    data = spm_squeeze(mean(Dtf(chanind, freqind, :), 1), [1 3]);
    
    if S.average
        data = mean(data, 1);
    end
    
    if D.fsample ~= Dtf.fsample
        [data, alpha] = spm_timeseries_resample(data, D.fsample/Dtf.fsample);
    else
        alpha = 1;
    end
    
    start = round(alpha*Dtf.indsample(D.time(1)));
    
    data = data(:, start:(start+D.nsamples-1));
    
else
    if D.ntrials ~= Dtf.ntrials
        error('Trial numbers should be equal between input and power dataset.');
    end
    
    if S.summarise
        timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
        if isempty(timeind) || any(isnan(timeind))
            error('Selected time window is invalid.');
        end
    else
        timeind = 1:D.nsamples;
    end
    
    data = spm_squeeze(mean(Dtf(chanind, freqind, timeind, :), 1), 1);
    
    if S.average
        data = mean(data, 1);
    end
    
    if S.summarise
        data = spm_squeeze(mean(data, 2), 2);
    else
        data = reshape(data, size(data, 1), []);
    end
end


data = data';

if S.standardize
    data = (data - repmat(mean(data, 1), size(data, 1), 1))./repmat(std(data, 1), size(data, 1), 1);
end

res.R     = data;

if S.average
    res.names = {S.regname};
else
    freq = Dtf.frequencies(freqind);
    for i = 1:length(freq)
        res.names{i} = sprintf('%s_%.1fHz', S.regname, freq(i));
    end
end


spm('FigName','Time-frequency power regressors: done');

