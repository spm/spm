function res = bf_features_cov(BF, S)
% Simple band limited covariance computation
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes  
% $Id: bf_features_cov.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.val = {[0 Inf]};
    foi.help = {'Frequency windows within which to compute covariance over (sec)'};
    
    taper = cfg_menu;
    taper.tag = 'taper';
    taper.name = 'Windowing';
    taper.help = {'Select a window for pre-multiplying the data'};
    taper.labels = {'Hanning', 'None'};
    taper.values = {'hanning', 'none'};
    taper.val = {'hanning'};
    
    cov      = cfg_branch;
    cov.tag  = 'cov';
    cov.name = 'Band limited covariance';
    cov.val  = {foi, taper};
    
    res = cov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;


ntrials = length(S.trials);
nchans  = length(S.channels);
%% now identify frequency bands of interest

nbands = size(S.foi,1);

if length(unique(cellfun(@length, S.samples)))~=1,
    error('all windows must be of equal length');
end;

nwoi            = numel(S.samples);
nsamples        = length(S.samples{1}); %% use length of first window to set up DCT (as all windows fixed at same length)
windowduration  = nsamples/D.fsample;
dctfreq         = (0:nsamples-1)/2/windowduration;           % DCT frequencies (Hz)
dctT            = spm_dctmtx(nsamples,nsamples);

allfreqind=[];

for fband = 1:nbands, %% allows one to break up spectrum and ignore some frequencies
    
    freqrange  = S.foi(fband,:);
    
    j          = find( (dctfreq >= freqrange(1)) & (dctfreq<=freqrange(2)));
  
    allfreqind = sort(unique([allfreqind j]));
    
end; % for fband=1:Nbands

% Hanning operator (if requested)
%----------------------------------------------------------------------

switch lower(S.taper)
    case 'hanning',
        W  = repmat(spm_hanning(nsamples)',nchans,1);
    case 'none'
        W  = ones(nchans,nsamples);
end

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

YY    = 0;
Tband = dctT(:,allfreqind); % filter to this band
for i = 1:ntrials
    for j = 1:nwoi
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));        
        
        Y = detrend(Y', 'constant')';
        
        Y = Y.*W;
        
        dctY = Y*Tband; %% frequency representation
        
        YY = YY+(dctY*dctY');        
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

N = ntrials*nsamples*nwoi;

C = YY./N; 

features.C = C;
features.N = N;

res = features;