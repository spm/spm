function [D] = spm_eeg_invert_classic_volumetric(D,val)
%% Volumetric multiple sparse priors
% This version only handles single subject single modality data; the removal
% of many scaling factors makes it easier to compare between forward models
%
% Note also that this funtion performs a bit differently from the
% spm_eeg_invert and spm_eeg_invert_classic functions. Namely:
%
% 1) No temporal filtering is carried out to the data by default
% 2) By default, each lead field element has one associated prior (i.e. no
%    "patches" or graph Laplacians are calculated).
% 3) Loreta-like priors/inversions are not (currently) supported.
%
% Ryan Timms and Gareth Barnes, 2023.
%
% Requires:
%           An SPM object, D
%           An inversion value, val
%
% The usual SPM invert shenanigans applies:
% D{i}.inv{val}.inverse:
%     inverse.modality - modality to use in case of multimodal datasets
%     inverse.trials - D.events.types to invert
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
%                      'IID' minimum norm
%                       'EBB' for empirical bayes beamformer
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.lpf    - band-pass filter - low frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter - high frequency cut-off (Hz)
%     inverse.Han    - switch for Hanning window
%     inverse.Nm     - maximum number of channel modes
%     inverse.Nmax     - maximum number of temporal modes
%     inverse.Nt     - fixed/requested number of temporal modes
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Qe     - any sensor error components (e.g. empty-room data)
%     inverse.Qe0     - minimum amount of sensor noise power relative to
%                        signal eg 0.1 would correspond to power SNR of 10.0
%     inverse.A       - predefined spatial modes (Nchans*Nmodes) to project
%                       sensor data through
%
% Evaluates:
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
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors. No symmetric priors are
% used in this implementation (just single patches) There is an option for
% a Beamforming prior : inversion type 'EBB'.


%==========================================================================
% Check data: only single subject and single modality supported
%==========================================================================
modalities = check_subjects_and_modality(D, val);

%==========================================================================
% Unpack D object: get the inversion index and inversion struct.
%==========================================================================
if nargin > 1
    D.val = val;
elseif ~isfield(D, 'val')
    D.val = 1;
end
val=D.val;
inverse   = D.inv{val}.inverse;

%==========================================================================
% Check function arguments: assign defaults if necessary
%==========================================================================
try     type            = inverse.type;           catch, type            = 'GS';         end
try     hpf             = inverse.hpf;            catch, hpf             = 256;          end
try     lpf             = inverse.lpf;            catch, lpf             = 0;            end
try     sdv             = inverse.sdv;            catch, sdv             = 4;            end
try     Han             = inverse.Han;            catch, Han             = 1;            end
try     woi             = inverse.woi;            catch, woi             = [];           end
try     Nmax            = inverse.Nmax;           catch, Nmax            = 512;          end % max number of temporal modes
try     Nm              = inverse.Nm;             catch, Nm              = [];           end
try     Nt              = inverse.Nt;             catch, Nt              = [];           end %% fixed/requested number of temporal modes
try     Ip              = inverse.Ip;             catch, Ip              = [];           end
try     QE              = inverse.QE;             catch, QE              = 1;            end % empty room noise measurement
try     Qe0             = inverse.Qe0;            catch, Qe0             = exp(-5);      end %% set noise floor at 1/100th signal power i.e. assume amplitude SNR of 10
try     inverse.A       = inverse.A;              catch, inverse.A       = [];           end %% orthogonal channel modes
try     no_temporal_filter = inverse.no_temporal_filter; catch, no_temporal_filter = 1;   end
try     complexind      = inverse.complexind;     catch, complexind      = [];           end

%==========================================================================
% Account for bad channels: get the indices of good channels
%==========================================================================
Ic  = setdiff(D.indchantype(modalities), badchannels(D));

%==========================================================================
% Lead fields: Load and assign values to number of dipoles (Nd). Exclude
% any bad channels at this stage.
%==========================================================================
fprintf('\nLoading or creating lead field matrix')
[L,~] = spm_eeg_lgainmat(D,[],D.chanlabels(Ic)); % Generate/load lead field
Nd=size(L,2); % Number of sources in lead field matrix
Np = Nd; % Number of priors

%==========================================================================
% Calculate the number of priors
%==========================================================================
if ~isempty(Ip)
    Np   = length(Ip);
else
    Ip=ceil([1:Np]*Nd/Np);
end


%==========================================================================
% Spatial projectors: construct a spatial projector matrix, A, and apply
% this to the forward model, i.e. eliminate low SNR spatial modes
%==========================================================================
[A, UL, Is, Ns] = construct_apply_spatial_projector(inverse, Nm, L, Nd);

%==========================================================================
% Time-window of interest (in milliseconds)
%==========================================================================
[w, It] = get_time_window_of_interest(woi, D);

%==========================================================================
% Construct temporal filter
%==========================================================================
if ~no_temporal_filter
    [pst, dct, Nb, qV, T] = construct_temporal_filter(D, It, sdv,lpf, hpf);
else
    T=eye(length(It));
    qV=T;
    pst=0;dct=0;Nb=1;
end 

%==========================================================================
% Construct Hanning window
%==========================================================================
W = construct_Hanning_window(Han, Nb);

%==========================================================================
% Get trials (a.k.a conditions)
%==========================================================================
[trial, Ntrialtypes] = get_trials(D);

%==========================================================================
% Get temporal covariance (Y'*Y)
%==========================================================================
[YY, badtrialind, Ik, Y] = get_temporal_covariance(D, Ntrialtypes, ...
    trial, complexind, Ic, It, A);

%==========================================================================
% Apply any Hanning and filtering
%==========================================================================
YY  = W' * YY * W;     % Hanning
YTY = T' * YY * T;     % Filter

%==========================================================================
% Get the temporal modes from the windowed and filtered data
%==========================================================================
[U, Nr, V, VE] = get_temporal_modes(Nt, YTY, Nmax);

%==========================================================================
% Apply temporal projector
%==========================================================================
S      = T*V;                           % temporal projector
Vq     = S*pinv(S'*qV*S)*S';            % temporal precision

%==========================================================================
% Get Spatial Covariance: Y*Y' for Gaussian process model. n.b. this is
% where the spatial transformation is applied to the sensor level data,
% i.e. Y -> AY
%==========================================================================
[AYYA, Nn, AY,UY, Y] = get_spatial_covariance(Ntrialtypes, D, trial,...
    badtrialind, complexind, Ic, It, Y, S, A, Nr);

%==========================================================================
% Data quality check: ensure that the data are full rank
%==========================================================================
check_data(AYYA, A);

% assuming equal noise over subjects (Qe) and modalities AQ
%--------------------------------------------------------------------------
AQeA   = A*QE*A';           % Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space
Q0     = Qe0*trace(AYYA)*Qe{1}./sum(Nn); %% fixed (min) level of sensor space variance

%==========================================================================
% Create Source Priors
%==========================================================================
if contains(type,'EBBr')
    reglevel=str2num(type(5:end));
    fprintf('\n Using regularizing beamformer prior to keep %d percent variance\n',reglevel)
    type='EBBr';
end

switch(type)
    case {'MSP','GS','ARD'}
        [Qp, LQpL] = build_MSP_source_priors(Np, Nd, UL);
    case {'EBBr'}
        [Qp, LQpL] = build_EBBr_source_priors(Ns, AYYA, reglevel, UL);
    case {'EBB'}
        [Qp, LQpL] = build_EBB_source_priors(Ns, AYYA, UL);
    case {'IID','MMN'}
        [Qp, LQpL] = build_IID_source_priors(Ns,UL);
end
fprintf('Using %d spatial source priors provided\n',length(Qp));

%==========================================================================
% Step 1: Run the inference for the first time.
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};
[QP, LQP, LQPL] = run_first_inference(type, Qp, Ns, AY, UL, Qe, QP, LQP,...
    LQPL, AYYA, LQpL, Nn, Q0);

%==========================================================================
% Step 2: Run the inference for the second time, having collapsed all
% priors into one summed prior.
%==========================================================================
[Np, Ne, Cy, h, F] = run_second_inference(LQPL, Qe, Q0, AYYA, Nn);

%==========================================================================
% Construct MAP weights matrix which maps from sensor to source data
%==========================================================================
[Cp, LCp, M] = build_MAP_inversion(h, Ne, Np, QP, LQP, Cy);

% conditional variance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%----------------------------------------------------------------------
Cq    = Cp - sum(LCp.*M')';

% evaluate conditional expectation
%----------------------------------------------------------------------
% evaluate conditional expectation (of the sum over trials)
%----------------------------------------------------------------------
SSR   = 0;
SST   = 0;
J     = {};

for j = 1:Ntrialtypes
    % trial-type specific source reconstruction
    %------------------------------------------------------------------
    J{j} = M*UY{j};
    % sum of squares
    %------------------------------------------------------------------
    SSR  = SSR + sum(var((UY{j} - UL*J{j}))); % changed variance calculation
    SST  = SST + sum(var( UY{j}));
end


% accuracy; signal to noise (over sources)
%======================================================================
R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE));

%==========================================================================
% Get a unique ID for the filtered data which the inference was run on
%==========================================================================
ID    = spm_data_id(AY);

%==========================================================================
% Save Results
%======================================================================
inverse.type     = type;                 % Inverse model
inverse.M        = M;                    % MAP projector (reduced)
inverse.J        = J;                    % Conditional expectation
inverse.Y        = Y;                    % ERP data (reduced)
inverse.L        = UL;                   % Lead field (reduced)
inverse.qC       = Cq;                   % Spatial covariance
inverse.tempU    = U;                    % Temporal SVD
inverse.E        = V;                    % Temporal modes
inverse.qV       = Vq;                   % Temporal correlations
inverse.T        = S;                    % Temporal projector
inverse.U        = {A};                  % Spatial projector
inverse.Is       = Is;                   % Indices of active dipoles
inverse.It       = It;                   % Indices of time bins
inverse.Ik       = Ik;                   % Indices of trials used
try
    inverse.Ic{1} = Ic;                  % Indices of good channels
catch
    inverse.Ic    = Ic;                 
end
inverse.Nd       = Nd;                   % Number of dipoles
inverse.pst      = pst;                  % Peristimulus time
inverse.dct      = dct;                  % Frequency range
inverse.F        = F;                    % Log-evidence
inverse.ID       = ID;                   % Data ID
inverse.R2       = R2;                   % Variance explained (reduced)
inverse.VE       = R2*VE;                % Variance explained
inverse.woi      = w;                    % Time-window inverted
inverse.Ip       = Ip;                   % Patch locations
inverse.modality = modalities;           % Modalities inverted


% save in struct
%----------------------------------------------------------------------
D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';


return

function [A, UL, Is, Ns] = construct_apply_spatial_projector(inverse,...
    Nm, L, Nd)
% Carries out a PCA on the LL' (aka the Gram matrix) and reduces the lead
% field matrix.
%
% Input arguments:
%   inverse - Inverse model parameters
%   Nm      - Number of spatial modes (if not specified, determined automatically)
%   L       - Lead field matrix
%   Nd      - Number of dipoles
%
% Returns:
%   A       - Spatial projector matrix
%   UL      - Reduced lead field matrix
%   Is      - Indices of the sources to be modeled
%   Ns      - Number of sources to be modeled

if isempty(inverse.A) % no spatial modes pre-specified
    if isempty(Nm) %% number of modes not specified
        [U, ~, ~] = spm_svd((L * L'), exp(-16));
        A         = U';                % spatial projector A
        UL        = A * L;
    else % number of modes pre-specified
        [U, ss, ~] = spm_svd((L * L'), 0);
        if length(ss) < Nm
            disp('number available');
            length(ss)
            error('Not this many spatial modes in lead fields');
        end
        disp('using preselected number spatial modes !');
        A         = U(:, 1:Nm)';        % spatial projector A
        UL        = A * L;
    end
else %% U was specified in input
    disp('Using pre-specified spatial modes');
    if isempty(Nm)
        error('Need to specify number of spatial modes if U is pre-specified');
    end
    
    A   = inverse.A;
    UL  = A * L;
end

Nm          = size(UL, 1);            % Number of spatial projectors
clear ss;

Is          = 1:Nd;                   % Indices of active dipoles - all of them.
Ns          = length(Is);             % Number of sources, Ns
fprintf('Using %d spatial modes', Nm)


function [w, It] = get_time_window_of_interest(woi, D)
% Get the time window of interest in samples and indices.
%
% This function calculates the time window of interest in samples and its
% corresponding indices based on the specified window of interest (woi) and
% the input data struct (D).
%
% Input arguments:
%   woi - Window of interest in milliseconds (1x2 vector) [start_time,
%   end_time] 
%   D   - Input data struct.
%
% Returns:
%   w   - Time window of interest in milliseconds [start_time, end_time] 
%   It  - Indices corresponding to the time window of interest
%
if isempty(woi)
    w      = 1000*[min(D.time) max(D.time)];
else
    w=woi;
end
It     = (w/1000 - D.timeonset)*D.fsample + 1;
It     = max(1,It(1)):min(It(end), length(D.time));
It     = fix(It);
fprintf('\nNumber of samples %d',length(It))

function [trial, Ntrialtypes] = get_trials(D)
% Get trial information and the number of trial types from the input data.
%
% This function extracts trial information and calculates the number of
% trial types from the provided data struct (D).
%
% Input arguments:
%   D - Input data struct containing trial-related information
%
% Returns:
%   trial         - Trial information (cell array of trial types)
%   Ntrialtypes   - Number of trial types
%
try
    trial = D.inv{D.val}.inverse.trials;
catch
    trial = D.condlist;
end
Ntrialtypes=length(trial);

function [YY, badtrialind, Ik, Y] = get_temporal_covariance(D,...
    Ntrialtypes, trial, complexind, Ic, It, A)
% Calculate the temporal covariance matrix and related values.
%
% This function calculates the temporal covariance matrix (YY) based on the
% provided input data struct (D) and various parameters.
%
% Input arguments:
%   D            - Input data struct
%   Ntrialtypes  - Number of trial types
%   trial        - Trial information (cell array of trial types)
%   complexind   - Complex channel indices (optional)
%   Ic           - Indices of good channels
%   It           - Indices of time bins
%   A            - Spatial projector matrix
%
% Returns:
%   YY           - Temporal covariance matrix
%   badtrialind  - Indices of bad trials
%   Ik           - Indices of trials used
%   Y            - Projected data
%

YY=0; % instantiate value of temporal covariance
N=0; % number of trials used in covariance calculation
i=sqrt(-1); % imaginary number, for use with complex data
badtrialind=D.badtrials;
Ik=[]; %% keep a record of trials used
for j = 1:Ntrialtypes                          % pool over conditions
    c     = D.indtrial(trial{j});     % and trials
    [~,ib]=intersect(c,badtrialind); % remove bad trials with the indices ib if there are any
    c=c(setxor(1:length(c),ib));
    Ik=[Ik c];
    Nk    = length(c);
    
    for k = 1:Nk
        if isempty(complexind)
            data=D(Ic,It,c(k));
        else
            data=squeeze(D(Ic,complexind(1,:),c(k))+i.*D(Ic,complexind(2,:),c(k)));
        end
        Y     = A*data;
        YY    = YY + Y'*Y;
        N     = N + 1;
    end
end
YY=YY./N;

function [AYYA, Nn, AY,UY, Y] = get_spatial_covariance(Ntrialtypes,...
    D, trial, badtrialind, complexind, Ic, It, Y, S, A, Nr)
% Calculate spatial covariance matrices and related values.
%
% This function calculates spatial covariance matrices (AYYA, UYYU) and
% various related values based on the provided input data struct (D) and
% parameters.
%
% Input arguments:
%   Ntrialtypes  - Number of trial types
%   D            - Input data struct
%   trial        - Trial information (cell array of trial types)
%   badtrialind  - Indices of bad trials
%   complexind   - Complex channel indices (optional)
%   Ic           - Indices of good channels
%   It           - Indices of time bins
%   Y            - Projected data
%   S            - Temporal projector matrix
%   A            - Spatial projector matrix
%   Nr           - Number of repetitions (samples per trial)
%
% Returns:
%   AYYA         - Pooled response for ReML (spatial covariance)
%   Nn           - Total number of samples
%   AY           - Pooled response for MVB (matrix)
%   UY           - Condition-specific ERP (cell array)
%   Y            - Projected data
%

UYYU = 0;  % Initialize second-order response
AYYA = 0;  % Initialize pooled response for ReML
Nn   = 0;  % Initialize total number of samples
AY = {};   % Initialize pooled response for MVB
UY = {};   % Initialize condition-specific ERP

i = sqrt(-1);  % Imaginary unit

for j = 1:Ntrialtypes % loop over Ntrialtypes trial types
    UY{j} = sparse(0);
    c       = D.indtrial(trial{j});
    [~,ib]=intersect(c,badtrialind); %% remove bad trials ib if there are any
    c=c(setxor(1:length(c),ib));
    Nk    = length(c);
    % loop over epochs
    %------------------------------------------------------------------
    for k = 1:Nk
        % stack (scaled aligned data) over modalities
        %--------------------------------------------------------------
        if isempty(complexind)
            data=D(Ic,It,c(k));
        else
            data=D(Ic,complexind(1,:),c(k))+i.*D(Ic,complexind(2,:),c(k));
        end
        Y       = data*S; %% in temporal subspace
        
        % Now we apply the spatial transformation to the sensor level data,
        % Y
        Y=A*Y; %%  in spatial subspace

        % accumulate first & second-order responses
        %--------------------------------------------------------------
        Nn       = Nn + Nr;         % number of samples
        YY          = Y*Y';                  % and covariance
        % accumulate statistics
        %--------------------------------------------------------------
        UY{j}     = UY{j} + Y;           % condition-specific ERP
        UYYU     = UYYU + YY;
        % and pool for optimisation of spatial priors
        %--------------------------------------------------------------
        AY{end + 1} = Y;                     % pooled response for MVB
        AYYA        = AYYA    + YY;          % pooled response for ReML
    end
end

AY=spm_cat(AY); %% goes to MVB/GS algorithm

function [pst, dct, Nb, qV, T] = construct_temporal_filter(D,...
    It, sdv, lpf, hpf)
% Construct temporal filter parameters.
%
% Input arguments:
%   D     - Input data struct
%   It    - Indices of time bins
%   sdv   - Standard deviation parameter for serial correlations
%   lpf   - Low-pass frequency for DCT filter
%   hpf   - High-pass frequency for DCT filter
%
% Returns:
%   pst   - Peristimulus time (ms)
%   dct   - DCT frequencies (Hz)
%   Nb    - Number of time bins
%   qV    - Samples* samples covariance matrix
%   T     - Temporal projector matrix

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

T      = spm_dctmtx(Nb,Nb);
j      = find( (dct >= lpf) & (dct <= hpf) ); % This is the wrong way round but leave for now for compatibility with spm_eeg_invert
T      = T(:,j);                    % Apply the filter to discrete cosines
dct    = dct(j);                    % Frequencies accepted

function modalities = check_subjects_and_modality(D, val)
% Check to see how many subjects are being fed to the function. This only
% works for a single subject.
Nl = length(D);
if Nl>1
    error('function only defined for a single subject');
end

% Check modalities - this function only works for single modality
%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;
if size(modalities,1)>1
    error('not defined for multiple modalities');
end

function W = construct_Hanning_window(Han, Nb)
% Construct Hanning window if requested
%
% Input arguments:
%   Han   - Flag indicating whether to use Hanning window (1) or not (0)
%   Nb    - Number of time bins
%
% Returns:
%   W     - Hanning window
%
if Han
    W  = sparse(1:Nb,1:Nb,spm_hanning(Nb)); % Use hanning unless specified
else
    W=1;
end

function [U, Nr, V, VE] = get_temporal_modes(Nt, YTY, Nmax)
if isempty(Nt) %% automatically assign appropriate number of temporal modes
    [U, E]  = spm_svd(YTY,exp(-8));          % get temporal modes
    if isempty(U) %% fallback
        warning('nothing found using spm svd, using svd');
        [U, E]  = svd(YTY);          % get temporal modes
    end
    E      = diag(E)/trace(YTY);            % normalise variance
    Nr     = min(length(E),Nmax);           % number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else %% use predefined number of modes
    [U, E]  = svd(YTY);          % get temporal modes
    E      = diag(E)/trace(YTY);            % normalise variance
    disp('Fixed number of temporal modes');
    Nr=Nt;
end

V      = U(:,1:Nr);                     % temporal modes
VE     = sum(E(1:Nr));                  % variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

function [Qp, LQpL] = build_MSP_source_priors(Np, Nd, UL)
% Create MSP spatial basis set in source space
%
% Input arguments:
%   Np    - Number of basis vectors to create
%   Nd    - Dimensionality of the source space
%   UL    - Reduced lead field matrix for linear transformation
%
% Output arguments:
%   Qp    - Cell array containing the created basis vectors
%   LQpL  - Cell array containing linear transformed basis vectors
%

%------------------------------------------------------------------
Qp    = {};
LQpL  = {};
for i = 1:Np
    q=sparse(Nd,1);
    q(i)=1;
    Qp{end + 1}.q   = q;
    LQpL{end + 1}.q = UL*q;
end

function [Qp, LQpL] = build_EBBr_source_priors(Ns, AYYA, reglevel, UL)
% create SMOOTH beamforming prior.
disp('NB  regularizing EBB algorithm, no smoothing !');
%------------------------------------------------------------------
QG=speye(Ns,Ns); %% was surface smoothing but not used in volumetric approach
[u1,s1,~]=svd(AYYA);
eigsum=cumsum(diag(s1))./sum(diag(s1));

usecomp=max(find(eigsum<=reglevel/100));
plot(1:length(s1),eigsum,usecomp,eigsum(usecomp),'*')
fprintf('\nEBB, Keeping %d components\n',usecomp);
ayya=u1(:,1:usecomp)'*AYYA*u1(:,1:usecomp);

InvCov = spm_inv(ayya);
allsource = sparse(Ns,1);
Sourcepower = sparse(Ns,1);
for bk = 1:Ns
    q               = QG(:,bk);
    
    smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
    smthlead=u1(:,1:usecomp)'*smthlead;
    normpower = 1/(smthlead'*smthlead);
    Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
    allsource(bk) = Sourcepower(bk)./normpower;
end
allsource = allsource/max(allsource);   % Normalise
subplot(2,1,1);
plot(Sourcepower);
subplot(2,1,2)
plot(allsource);
Qp{1} = diag(allsource);
LQpL{1} = UL*diag(allsource)*UL';

function [Qp, LQpL] = build_IID_source_priors(Ns, UL)
% create minimum norm prior
%------------------------------------------------------------------
Qp{1}   = speye(Ns,Ns);
LQpL{1} = UL*UL';

function [Qp, LQpL] = build_EBB_source_priors(Ns, AYYA, UL)
% create SMOOTH beamforming prior.
disp('NB  EBB algorithm, no smoothing !');
%------------------------------------------------------------------
QG=speye(Ns,Ns); %% was surface smoothing but not used in volumetric approach


InvCov = spm_inv(AYYA);
allsource = sparse(Ns,1);
Sourcepower = sparse(Ns,1);
for bk = 1:Ns
    q               = QG(:,bk);
    smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
    normpower = 1/(smthlead'*smthlead);
    Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
    allsource(bk) = Sourcepower(bk)./normpower;
end
allsource = allsource/max(allsource);   % Normalise
Qp{1} = diag(allsource);
LQpL{1} = UL*diag(allsource)*UL';

function [QP, LQP, LQPL] = accumulate_inference_results(QP, qp, LQP, UL, LQPL)
QP{end + 1}   = diag(qp);
LQP{end + 1}  = UL*qp;
LQPL{end + 1} = LQP{end}*UL';

function [QP, LQP, LQPL] = run_first_inference(type, Qp, Ns, AY, UL, Qe, QP, LQP, LQPL, AYYA, LQpL, Nn, Q0)
switch(type)
    case {'MSP','GS','EBBgs'}
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = zeros(Ns,Np); %% NB SETTING UP A NEW Q HERE
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        MVB   = spm_mvb(AY,UL,[],Q,Qe,16); %% Qe is identity with unit trace
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        % MVB.cp provides the final weights of the hyperparameters
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    case {'MSP','ARD'}
        % ReML / ARD inversion
        %------------------------------------------------------------------
        [~,h,~,~] = spm_sp_reml(AYYA,[],[Qe LQpL],Nn);
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            if hp(i) > max(hp)/128
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        [QP, LQP, LQPL] = accumulate_inference_results(QP, qp, LQP, UL, LQPL);
end

switch(type)
    case {'IID','MMN','EBB','EBBr'}
        % or ReML - ARD (This is where the inference on h is carried out)
        %------------------------------------------------------------------
        [~,h,~,~] = spm_reml_sc(AYYA,[],[Qe LQpL],Nn,-4,16,Q0);
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        [QP, LQP, LQPL] = accumulate_inference_results(QP, qp, LQP, UL, LQPL);
end

function check_data(AYYA, A)
if rank(AYYA)~=size(A,1)
    warning('AYYA IS RANK DEFICIENT');
end

function [Np, Ne, Cy, h, F] = run_second_inference(LQPL, Qe, Q0, AYYA, Nn)
% re-do ReML (with informative hyperpriors)
%----------------------------------------------------------------------
Np    = length(LQPL);       % Final number of priors
Ne    = length(Qe);         % Sensor noise prior
Q     = [{Q0} LQPL]; %% sensor corvariance prior:  Qe is identity with unit trace, LQPL is in the units of data
[Cy,h,~,F]= spm_reml_sc(AYYA,[],Q,Nn,-4,16,Q0);

function [Cp, LCp, M] = build_MAP_inversion(h, Ne, Np, QP, LQP, Cy)
fprintf('Building reconstruction weights for subject 1\n')
% Recalculate F here
Cp    = sparse(0);
LCp   = sparse(0);
hp    = h(Ne + (1:Np));
for j = 1:Np
    Cp  =  Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

% MAP estimates of instantaneous sources
%======================================================================
% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
M     = LCp'/Cy;