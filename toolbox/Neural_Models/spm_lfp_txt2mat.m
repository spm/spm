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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Vladimir Litvak 
% $Id: spm_lfp_txt2mat.m 2255 2008-09-30 15:36:59Z vladimir $


% load data (in this case a .mat file)
%--------------------------------------------------------------------------
load amy_hc_0501_R1

% define the output file name
%--------------------------------------------------------------------------
fname = 'LFP_example';

% define the sampling rate
%--------------------------------------------------------------------------
fsample = 1000;

% define epochs and create data array - this is specific for this data
% replace with your own code
%--------------------------------------------------------------------------
bins  = find(diff(data(:,3)) > 4);
bins  = bins([3,9]);
nbins = 1000*10;

% Define the channels of interest - in this case only 3,4 and 5
%--------------------------------------------------------------------------
Ic = [3 4 5];

% Create the Fieldtrip raw struct
%--------------------------------------------------------------------------
ftdata = [];

for i = 1:length(bins)
   It            = [1:nbins] + bins(i);
   ftdata.trial{i} = data(It, Ic)';
   ftdata.time{i} = [0:(nbins-1)]./fsample;
end


ftdata.fsample = fsample;
ftdata.label = colheaders(Ic);
ftdata.label = ftdata.label(:);

% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, fname);

% Examples of providing additional information in a script
% [] comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
D = type(D, 'single');             % Sets the dataset type
D = chantype(D, [], 'LFP');        % Sets the channel type 
D = conditions(D, [], 'Sound 1');  % Sets the condition label

% save
%--------------------------------------------------------------------------
save(D);