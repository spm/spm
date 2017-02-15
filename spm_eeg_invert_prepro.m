function [D] = spm_eeg_invert_prepro(D,val)
%
%% Preprocessing for inversion stage. 
%% includes spatial and temporal dimension reduction
%
%% this version only handles single subject single modality data
%% the removal of many scaling factors makes it easier to compare between forward models
% ReML inversion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert_classic(D)
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG/MEG data
%
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.modality - modality to use in case of multimodal datasets
%
%     inverse.trials - D.events.types to invert
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' minimum norm
%                      'EBB' for empirical bayes beamformer

%    inverse.priors{} - a cell array of anatomical and functional priors to be
%    considered
%
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.lpf    - band-pass filter - low frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter - high frequency cut-off (Hz)
%     inverse.Han    - switch for Hanning window
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%
%     inverse.Nm     - maximum number of channel modes
%     inverse.Nr     - maximum number of temporal modes
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Na     - number of most energetic dipoles
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Qe     - any sensor error components (e.g. empty-room data)
%     inverse.Qe0     - minimum amount of sensor noise power relative to
%                        signal eg 0.1 would correspond to power SNR of 10.0
%     inverse.A       - predefined spatial modes (Nchans*Nmodes) to project
%                       sensor data through
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J{i}   - Conditional expectation (i conditions) J = M*U*Y
%     inverse.L      - Lead field (reduced UL := U*L)
%     inverse.qC     - spatial covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal projector
%     inverse.U(j)   - spatial projector (j modalities) - derived from data
%     inverse.A      - pre-specified spatial projector
%     inverse.Y{i}   - reduced data (i conditions) UY = UL*J + UE
%     inverse.Is     - Indices of active dipoles
%     inverse.It     - Indices of time bins
%     inverse.Ic{j}  - Indices of good channels (j modalities)
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.VE     - variance explained in spatial/temporal subspaces (%)
%     inverse.R2     - variance in subspaces accounted for by model (%)
%     inverse.scale  - scaling of data for each of j modalities
%__________________________________________________________________________
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors.
% No symmetric priors are used in this implementation (just single patches)
% There is an option for a Beamforming prior : inversion type 'EBB'
% also added new beamforming method- using GS rather than ARD- from Juan
% David Martinez Vargas 'EBBgs'
%
% The code was used in
% Lopez, J. D., Penny, W. D., Espinosa, J. J., Barnes, G. R. (2012).
% A general Bayesian treatment for MEG source reconstruction incorporating
% lead field uncertainty.
% Neuroimage 60(2), 1194-1204 doi:10.1016/j.neuroimage.2012.01.077.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
 
% Jose David Lopez, Gareth Barnes, Vladimir Litvak
% $Id: spm_eeg_invert_prepro.m 7017 2017-02-15 12:50:58Z karl $


Nl = length(D);

if Nl>1
    error('function only defined for a single subject');
else
    D=D{1};
end

% D - SPM data structure
%==========================================================================
if nargin > 1
    D.val = val;
elseif ~isfield(D, 'val')
    D.val = 1;
end


val=D.val;

inverse   = D.inv{val}.inverse;

% forward model
%--------------------------------------------------------------------------


% defaults
%--------------------------------------------------------------------------

try, type = inverse.type;   catch, type = 'GS';     end
try, s    = inverse.smooth; catch, s    = 0.6;      end
try, Np   = inverse.Np;     catch, Np   = 256;      end
try, Nr   = inverse.Nr;     catch, Nr   = 16;       end % requested number of temporal modes, could be changed depending on svd
try, xyz  = inverse.xyz;    catch, xyz  = [0 0 0];  end
try, rad  = inverse.rad;    catch, rad  = 128;      end
try, hpf  = inverse.hpf;    catch, hpf  = 48;       end % need to one day put these the correct way round
try, lpf  = inverse.lpf;    catch, lpf  = 0;        end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, Nm   = inverse.Nm;     catch, Nm   = [];       end
try, Nt   = inverse.Nt;     catch, Nt   = [];       end % fixed number of temporal modes
try, Ip   = inverse.Ip;     catch, Ip   = [];       end
try, QE   = inverse.QE;     catch, QE   = 1;        end %  empty room noise measurement
try, Qe0  = inverse.Qe0;    catch, Qe0  = exp(-5);  end % set noise floor at 1/100th signal power i.e. assume amplitude SNR of 10
try, inverse.A;             catch, inverse.A   = [];end % orthogonal channel modes


% get specified modalities to invert (default to all)
%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;    % MEG in this case
Nmax       = 16;                             % max number of temporal modes

% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================
fprintf('Checking leadfields')
[L,D] = spm_eeg_lgainmat(D);                 % Generate/load lead field
Nd    = size(L,2);


% Check gain or lead-field matrices
%------------------------------------------------------------------
if size(modalities,1)>1,
    error('not defined for multiple modalities');
end;

Ic{1}  = setdiff(D.indchantype(modalities), badchannels(D));
Nd    = size(L,2);      % Number of dipoles


%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================

fprintf('Optimising and aligning spatial modes ...\n')

% eliminate low SNR spatial modes
%------------------------------------------------------------------
%disp('FIXING FULL RANK A');
%inverse.A=eye(length(Ic)); Nm=size(inverse.A,1);

if isempty(inverse.A), % no spatial modes prespecified
    if isempty(Nm),    % number of modes not specifiedd
        [U,ss,vv] = spm_svd((L*L'),exp(-16));
        A         = U';                % spatial projector A
        UL        = A*L;
        
    else % number of modes pre-specified
        [U,ss,vv]    = spm_svd((L*L'),0);
        if length(ss)<Nm,
            disp('number available');
            length(ss)
            error('Not this many spatial modes in lead fields');
        end
        ss    = ss(1:Nm);
        disp('using preselected number spatial modes !');
        A     = U(:,1:Nm)';                 % spatial projector A
        UL    = A*L;
    end;
else %% U was specified in input
    disp('Using pre-specified spatial modes');
    if isempty(Nm),
        error('Need to specify number of spatial modes if U is prespecified');
    end
    A  = inverse.A;
    UL = A*L;
end;

Nm     = size(UL,1);         % Number of spatial projectors


clear ss vv

% Report
%----------------------------------------------------------------------
fprintf('Using %d spatial modes',Nm)


%==========================================================================
% Temporal projector
%==========================================================================
AY    = {};                                      % pooled response for MVB
AYYA  = 0;                                       % pooled response for ReML

% Time-window of interest
%----------------------------------------------------------------------

if isempty(woi)
    w  = 1000*[min(D.time) max(D.time)];
else
    w  = woi; %% in milliseconds
end;

It     = (w/1000 - D.timeonset)*D.fsample + 1;
It     = max(1,It(1)):min(It(end), length(D.time));
It     = fix(It);
disp(sprintf('Number of samples %d',length(It)));

% Peristimulus time
%----------------------------------------------------------------------
pst    = 1000*D.time;                   % peristimulus time (ms)
pst    = pst(It);                       % windowed time (ms)
dur    = (pst(end) - pst(1))/1000;      % duration (s)
dct    = (It - It(1))/2/dur;            % DCT frequencies (Hz)
Nb     = length(It);                    % number of time bins

% Serial correlations
%----------------------------------------------------------------------
K      = exp(-(pst - pst(1)).^2/(2*sdv^2)); %% sdv set to 4 by default
K      = toeplitz(K);
qV     = sparse(K*K'); %% Samples* samples covariance matrix- assumes smooth iid

% Confounds and temporal subspace
%----------------------------------------------------------------------
T      = spm_dctmtx(Nb,Nb);         % use plot(T) here
j      = find( (dct >= lpf) & (dct <= hpf) ); %% THis is the wrong way round but leave for nowfor compatibility with spm_eeg_invert
T      = T(:,j);                    % Apply the filter to discrete cosines
dct    = dct(j);                    % Frequencies accepted

% Hanning window
%----------------------------------------------------------------------
if Han
    W  = sparse(1:Nb,1:Nb,spm_hanning(Nb)); %% use hanning unless specified
else
    W  = 1;
end



% get trials or conditions
%----------------------------------------------------------------------
try
    trial = D.inv{D.val}.inverse.trials;
catch
    trial = D.condlist;
end
Ntrialtypes = length(trial);

% get temporal covariance (Y'*Y) to find temporal modes
%======================================================================
YY    = 0;
N     = 0;
badtrialind=D.badtrials;
for j = 1:Ntrialtypes,              % pool over conditions
    c     = D.indtrial(trial{j});   % and trials
    c     = setxor(c,badtrialind);  % ignore bad trials
    Nk    = length(c);
    for k = 1:Nk
        Y  = A*D(Ic{1},It,c(k));
        
        YY = YY + Y'*Y;
        N  = N + 1;
    end
end
YY  = YY./N;


% Apply any Hanning and filtering
%------------------------------------------------------------------
YY  = W'*YY*W;     % Hanning
YTY = T'*YY*T;     % Filter

%======================================================================
if isempty(Nt), %% automatically assign appropriate number of temporal modes    
    [U E]  = spm_svd(YTY,exp(-8));          % get temporal modes
    if isempty(U), %% fallback
        warning('nothing found using spm svd, using svd');
        [U E]  = svd(YTY);          % get temporal modes
    end;
    E      = diag(E)/trace(YTY);            % normalise variance
    Nr     = min(length(E),Nmax);           % number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else %% use predefined number of modes
    [U E]  = svd(YTY);          % get temporal modes
    E      = diag(E)/trace(YTY);            % normalise variance
    disp('Fixed number of temporal modes');
    Nr=Nt;
end;

V      = U(:,1:Nr);                     % temporal modes
VE     = sum(E(1:Nr));                  % variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

% projection and whitening
%----------------------------------------------------------------------
S      = T*V;                           % temporal projector
Vq     = S*pinv(S'*qV*S)*S';            % temporal precision

%disp('FIXIN TEMP PROJECTOR TO FULL RANK');
%S=eye(size(YY,1));

% get spatial covariance (Y*Y') for Gaussian process model
%======================================================================
% loop over Ntrialtypes trial types 
%----------------------------------------------------------------------
UYYU    = 0;
AYYA    = 0;
Nn      = 0;                             % number of samples
AY      = {};
Ntrials = 0;

for j = 1:Ntrialtypes,
    UY{j} = sparse(0);
    c     = D.indtrial(trial{j});
    c     = setxor(c,badtrialind); %% ignore bad trials
    Nk    = length(c);
       
    % loop over epochs
    %------------------------------------------------------------------
    for k = 1:Nk
        
        % stack (scaled aligned data) over modalities
        %--------------------------------------------------------------
        Y = D(Ic{1},It,c(k))*S; %% in temporal subspace
        Y = A*Y; %%  in spatial subspace
  
        
        % accumulate first & second-order responses
        %--------------------------------------------------------------
        Nn      = Nn + Nr;             % number of samples
        
        YY      = Y*Y';                % and covariance
        Ntrials = Ntrials+1;
        
        % accumulate statistics (subject-specific)
        %--------------------------------------------------------------
        UY{j}   = UY{j} + Y;           % condition-specific ERP
        UYYU    = UYYU + YY;           % subject-specific covariance
        
        % and pool for optimisation of spatial priors over subjects
        %--------------------------------------------------------------
        AY{end + 1} = Y;               % pooled response for MVB
        AYYA        = AYYA  + YY;      % pooled response for ReML
        
    end
end
AY  = spm_cat(AY);     %% goes to MVB/GS algorithm
ID  = spm_data_id(AY); %% get a unique ID for these filtered data


%======================================================================

inverse.AY=AY; %% concatenated data in reduced spatial modes 
inverse.AYYA=AYYA;%% sensor level covariance of all trials, conditions and samples (containing Nn temporal modes in total)

inverse.UY     = UY;                    % ERP data (reduced)
inverse.L      = UL;                   % Lead-field (reduced)
inverse.Nn     = Nn; %% number of independent samples (temporal modes* trials* conditions)
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal projector
inverse.U      = {A};                    % spatial projector
inverse.Ic     = Ic; %% good channels
inverse.It     = It; %% time indices
inverse.Nd     = Nd;                   % number of dipoles
inverse.pst    = pst;                  % peristimulus time
inverse.dct    = dct;                  % frequency range
inverse.ID     = ID;                   % data ID
inverse.woi    = w;                    % time-window inverted

inverse.modality = modalities;         % modalities inverted


% save in struct
%----------------------------------------------------------------------
D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';

disp('Done invert prepro');

return
