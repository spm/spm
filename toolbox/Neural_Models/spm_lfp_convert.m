% demo for creating SPM data structures: basically an array:
% channels x time x condition x epochs
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_convert.m 1507 2008-04-29 10:44:36Z vladimir $


% load data (in this case a .mat file)
%--------------------------------------------------------------------------
load amy_hc_0501_R1
 
% define epochs and create data array
%--------------------------------------------------------------------------
bins  = find(diff(data(:,3)) > 4);
bins  = bins([3,9]);
nbins = 1000*10;
for i = 1:length(bins)
    It            = [1:nbins] + bins(i);
    D.data(:,:,i) = data(It,:)';
end
 
% channels
%--------------------------------------------------------------------------
channels.Bad   = [];              % indices of bad channels
channels.name  = colheaders;      % a cell of strings
channels.eeg   = [4 5];           % channels that will be used
channels.order = [1:5];
 
    
% events
%--------------------------------------------------------------------------
events.code   = [1 2];                          % a code for each epoch
events.time   = bins + 1;                       % time bins
events.types  = [1 2];                          % a list of codes
events.start  = 1;                              % time bins to PST = 0
events.stop   = size(D.data,2);                 % time bins to end
events.Ntypes = length(events.types);
events.reject = zeros(1,length(events.code));
    
    
% D structure
%--------------------------------------------------------------------------    
D.channels  = channels;
D.events    = events;
 
D.Nchannels = size(D.data,1);
D.Nsamples  = size(D.data,2);
D.Nevents   = size(D.data,3);
D.Radc      = 1000;
D.modality  = 'LFP';
D.units     = '\muV';
D.fname     = 'LFP'
D.path      = pwd;
 
% save
%--------------------------------------------------------------------------
save(D)
