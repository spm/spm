% demo for creating an SPM M/EEG dataset from arbitrary data using 
% conversion of simple Fieldtrip raw data struct. 
% SPM8 internal format is quite complex but is transparent to the user
% as meeg object's methods take care of maintaining its consistency. 
% The most straightforward way to convert arbitrary data that is available
% as an ASCII file or *.mat file with some variables to SPM8 is to create
% a quite simple Fieldtrip raw data struct and then use SPM's
% spm_eeg_ft2spm to convert this struct to SPM8 file. Missing information
% can then be supplemented using meeg methods and SPM functions.
% Fieldtrip raw struct must contain the following fields:
% .fsample - sampling rate (Hz)
% .trial - cell array of matrices with identical dimensions channels x time
% .time - cell array of time vectors (in sec), the same length as the
%         second dimension of the data. For SPM8 they must be identical.
% .label - cell array of strings, list of channel labels. Same length as
%         the first dimension of the data.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak 
% $Id: spm_eeg_convert_arbitrary_data.m 5404 2013-04-12 15:08:57Z vladimir $


% Initialize SPM
%--------------------------------------------------------------------------
spm('defaults','EEG');

% Some details about the data
%--------------------------------------------------------------------------
Nchannels = 5;
Nsamples  = 500;
Ntrials   = 50;
TimeOnset = -0.1; % in sec
Fsample = 1000;

chlabels = {
            'LFP1'
            'LFP2' 
            'LFP3' 
            'LFP4' 
            'LFP5'
            };

% define the output file name
%--------------------------------------------------------------------------
fname = 'arbitrary_data_example';

% create data array 
%--------------------------------------------------------------------------
data = randn([Nchannels, Nsamples, Ntrials]);

% create the time axis (should be the same for all trials)
%--------------------------------------------------------------------------
timeaxis = [0:(Nsamples-1)]./Fsample + TimeOnset;

% Create the Fieldtrip raw struct
%--------------------------------------------------------------------------

ftdata = [];

for i = 1:Ntrials
   ftdata.trial{i} = squeeze(data(:, :, i));
   ftdata.time{i} = timeaxis;
end


ftdata.fsample = Fsample;
ftdata.label = chlabels;
ftdata.label = ftdata.label(:);

% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, fname);

% Examples of providing additional information in a script
% [] comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
D = type(D, 'single');                        % Sets the dataset type
D = chantype(D, ':', 'LFP');                   % Sets the channel type 
D = conditions(D, 1:Ntrials, 'Condition 1');  % Sets the condition label

% save
%--------------------------------------------------------------------------
save(D);

