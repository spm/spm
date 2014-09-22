function res = spm_eeg_regressors_chandata(S)
% Generate regressors from channel data
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
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_regressors_chandata.m 6186 2014-09-22 11:31:11Z vladimir $


SVNrev = '$Rev: 6186 $';

if nargin == 0
    %--------------------------------------------------------------------------
    % TF power dataset
    %--------------------------------------------------------------------------
    Dr        = cfg_files;
    Dr.tag    = 'Dr';
    Dr.name   = 'Regressor dataset name';
    Dr.filter = 'mat';
    Dr.num    = [1 1];
    Dr.help   = {'Select the M/EEG mat file containing regressor data'};
    
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
    
    chandata = cfg_branch;
    chandata.tag = 'chandata';
    chandata.name = 'Channel data';
    chandata.val = {Dr, spm_cfg_eeg_channel_selector, timewin};
    
    res = chandata;
    
    return
end

%-Startup
%--------------------------------------------------------------------------
spm('sFnBanner', mfilename, SVNrev);
spm('FigName','Channel data regressors');

if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end

if iscell(S.Dr)
    S.Dr = char(S.Dr);
end

Dr = spm_eeg_load(S.Dr);
D  = spm_eeg_load(S.D);

if isequal(Dr.transformtype, 'TF');
    error('Time domain dataset is expected as input.')
end

chanind = setdiff(Dr.selectchannels(spm_cfg_eeg_channel_selector(S.channels)), Dr.badchannels);

if isempty(chanind)
    error('No channels were selected');
end

res.R = [];

for i = 1:length(chanind)
    if isequal(D.type, 'continuous')
        
        if ~isequal(Dr.type, 'continuous') || (D.time(1) < Dr.time(1)) || (D.time(end)>Dr.time(end))
            error('All times of the input dataset should be within the power dataset.');
        end
        
        data = Dr(chanind(i), :);
        
        if D.fsample ~= Dr.fsample
            [data, alpha] = spm_timeseries_resample(data, D.fsample/Dr.fsample);
        else
            alpha = 1;
        end
        
        start = round(alpha*Dr.indsample(D.time(1)));
        
        data = data(:, start:(start+D.nsamples-1));
        
    else
        if D.ntrials ~= Dr.ntrials
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
        
        data = spm_squeeze(Dr(chanind(i), timeind, :), 1);
        
        if S.summarise
            data = spm_squeeze(mean(data, 1), 1);
        else
            data = reshape(data, 1, []);
        end
        
    end
    
    res.R     = [res.R data(:)];
end

res.names = Dr.chanlabels(chanind);


spm('FigName', 'Channel data regressors: done');
