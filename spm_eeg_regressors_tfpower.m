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

SVNrev = '$Rev: 6007 $';

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
    timewin.help    = {'Start and stop of the time window [ms]. (Used only for the epoched case)'};
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
    tfpower.val = {Dtf, spm_cfg_eeg_channel_selector, timewin, freqwin, regname};
    
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

freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
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
    
    data = spm_squeeze(mean(mean(Dtf(chanind, freqind, :), 2), 1), 2:4);
    
    if D.fsample ~= Dtf.fsample
        [data, alpha] = spm_resample(data, D.fsample/Dtf.fsample);
    else
        alpha = 1;
    end
    
    start = round(alpha*Dtf.indsample(D.time(1)));
    
    data = data(:, start:(start+D.nsamples-1));
    
else
    if D.ntrials ~= Dtf.ntrials
        error('Trial numbers should be equal between input and power dataset.');
    end
    
    timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
    if isempty(timeind) || any(isnan(timeind))
        error('Selected time window is invalid.');
    end
        
    data = squeeze(mean(mean(mean(Dtf(chanind, freqind, timeind, :), 3), 2), 1));
    
end

res.R     = data(:);

res.names = {S.regname};


spm('FigName','Time-frequency power regressors: done');

