function res = bf_features_tdcov(BF, S)
% Simple band limited covariance computation with temporal decomposition
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2020-2023 Wellcome Centre for Human Neuroimaging


if nargin == 0
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.val = {[0 Inf]};
    foi.help = {'Frequency windows within which to compute covariance over (sec)'};
    
    ntmodes = cfg_entry;
    ntmodes.tag = 'ntmodes';
    ntmodes.name = 'Number of temporal modes';
    ntmodes.strtype = 'i';
    ntmodes.num = [1 1];
    ntmodes.val = {[4]};
    ntmodes.help = {'Number of temporal modes, set to 0 for automatic determination'};
    
    taper = cfg_menu;
    taper.tag = 'taper';
    taper.name = 'Windowing';
    taper.help = {'Select a window for pre-multiplying the data'};
    taper.labels = {'Hanning', 'None'};
    taper.values = {'hanning', 'none'};
    taper.val = {'hanning'};
    
    tdcov      = cfg_branch;
    tdcov.tag  = 'tdcov';
    tdcov.name = 'Covariance with temporal decomposition';
    tdcov.val  = {foi, ntmodes, taper};
    
    res = tdcov;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;


ntrials = length(S.trials);
nchans  = length(S.channels);
% now identify frequency bands of interest

nbands = size(S.foi,1);

if length(unique(cellfun(@length, S.samples)))~=1
    error('all windows must be of equal length');
end

nwoi            = numel(S.samples);
nsamples        = length(S.samples{1}); %% use length of first window to set up DCT (as all windows fixed at same length)
windowduration  = nsamples/D.fsample;
dctfreq         = (0:nsamples-1)/2/windowduration;           % DCT frequencies (Hz)
dctT            = spm_dctmtx(nsamples,nsamples);

allfreqind=[];

for fband = 1:nbands %% allows one to break up spectrum and ignore some frequencies
    
    freqrange  = S.foi(fband,:);
    
    j          = find( (dctfreq >= freqrange(1)) & (dctfreq<=freqrange(2)));
    
    allfreqind = sort(unique([allfreqind j]));
    
end % for fband=1:Nbands

% Hanning operator (if requested)
%----------------------------------------------------------------------

switch lower(S.taper)
    case 'hanning'
        W  = repmat(spm_hanning(nsamples)',nchans,1);
    case 'none'
        W  = ones(nchans,nsamples);
end

% Temporal decomposition
%---------------------------------------------------------------------

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Prepating temporal decomposition'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

YY    = 0;
N=0;
Tband = dctT(:,allfreqind); % filter to this band
for i = 1:ntrials                    % pool over trials and windows
    for j = 1:nwoi
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));
        Y = detrend(Y', 'constant')';
        Y = Y.*W;
        dctY = Y*Tband; %% frequency representation
        YY    = YY + dctY'*dctY;
        N     = N + 1;
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end
YY=YY./N;

%

if ~S.ntmodes %% automatically assign appropriate number of temporal modes
    fprintf('Determining number of temporal modes required:\n');
    [U E]  = spm_svd(YY,exp(-8));          % get temporal modes
    if isempty(U) %% fallback
        warning('nothing found using spm svd, using svd');
        [U E]  = svd(YY);          % get temporal modes
    end
    E      = diag(E)/trace(YY);            % normalise variance
    Nr     = min(length(E),16);           % number of temporal modes
    Nr     = max(Nr,1);          %% use at least one mode
else %% use predefined number of modes
    [U E]  = svd(YY);          % get temporal modes
    E      = diag(E)/trace(YY);            % normalise variance
    disp('Fixed number of temporal modes');
    Nr=S.ntmodes;
end

V      = U(:,1:Nr);                     % temporal modes
VE     = sum(E(1:Nr));                  % variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

nsamples = Nr;

% Final Covariance preparation
%---------------------------------------------------------------------

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials, 'Computing covariance'); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end

YY    = 0;
UY    = 0;
Tband = dctT(:,allfreqind); % filter to this band
for i = 1:ntrials
    for j = 1:nwoi
        Y  = squeeze(D(S.channels, S.samples{j}, S.trials(i)));
        
        Y = detrend(Y', 'constant')';
        
        Y = Y.*W;
        
        dctY = (Y*Tband)*V; %% time-frequency representation
        
        YY = YY+(dctY*dctY');
        UY = UY + dctY; % Pooled signal for trial average 'evoked response'
    end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end



spm_progress_bar('Clear');

N = ntrials*nsamples*nwoi;

C = YY./N;

features.C = C;     % Covariance 
features.N = N;     % Nsamples
% features.V = V;     % Temporal projector
% features.W = W;     % Window
% features.T = Tband; % Frequency filter
features.VE = VE;   % Variance explained
features.UY = UY;   % Trial summed 'evoked response' 

res = features;
