function [D] = spm_eeg_invert(D, val)
% ReML inversion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert(D)
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
%                      'IID' LORETA and minimum norm
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.lpf    - band-pass filter - low frequency cutoff (Hz)
%     inverse.hpf    - band-pass filter - high frequency cutoff (Hz)
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
%     inverse.pQ     - any source priors (eg from fMRI); vector or matrix
%     inverse.Qe     - any sensor error components (eg empty-room data)
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J      - Conditional expectation
%     inverse.L      - Lead field (reduced)
%     inverse.R      - Re-referencing matrix
%     inverse.qC     - spatial covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal subspace
%     inverse.U      - spatial subspace
%     inverse.R      - re-referencing matrix
%     inverse.Y      - reduced data
%     inverse.Is     - Indices of active dipoles
%     inverse.It     - Indices of time bins
%     inverse.Ic     - Indices of good channels
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.VE     - variance explained in spatial/temporal subspaces (%)
%     inverse.R2     - variance in subspaces accounted for by model (%)
%     inverse.Lscale - scaling of leadfields for each modality
%     inverse.Yscale - scaling of data for each modality

%     inverse.
%__________________________________________________________________________
% 
% 1. This routine implements "group-based" inversion, corresponding to
% ill-posed linear models of the following form:
%
% [AY{1}...AY{n}] = L(1} * [J{1}...J{n}]   +  [e{1}...e{n}]
%
% where AY{i} are the spatially normalized or adjusted data from subject i
% that would have been seen if the lead-field L{i} = L{1}. The ensuing
% Gaussian process priors on sources are then used to estimate subject-specific
% MAP estimates of J{i} using
%
% Y{i} = L(i} * J{i}  +  e{i}
%
% Here, A{i}  = L{1}*pinv(L{i}) => 
%       AY{i} = A(i}*L(i}*J{i}
%             = L(1}*J{i}
%
% See: Electromagnetic source reconstruction for group studies.
% Litvak V, Friston K.
% NeuroImage. 2008 Oct 1;42(4):1490-8.
%
%__________________________________________________________________________
% 
% 2. It also implements "fusion" of different types of MEG and EEG data,
% corresponding to ill-posed linear models of the following form:
%
%             y{1,...,t}  = L(1} * J{1,...,t}   +  e{{1,...,t}}
%             y{2,...,t}  = L(2}                   e{{2,...,t}}
%                  .
%                  .
%                  .
%             y{n,...,t}  = L(n}                   e{{n,...,t}}
%
% Under empirical priors on J{1,...,t} for n data-sets with t trial types.
%
% See: MEG and EEG data fusion: Simultaneous localisation of face-evoked 
% responses.
% Henson R, Mouchlianitis E & Friston K.
% Neuroimage. 2009. 47:581-9.
%__________________________________________________________________________
% 
% 3. It also allows incorporation of spatial source priors, eg, from fMRI
% (see spm_eeg_inv_fmripriors.m)
%
% See: A Parametric Empirical Bayesian framework for fMRI-constrained 
% MEG/EEG source reconstruction.
% Henson R, Flandin G, Friston K & Mattout J.
% Human Brain Mapping (in press).
%__________________________________________________________________________
%
% The routine essentially consists of two steps:
%
%   1. Optimisation of spatial source priors over subjects
%   2. Re-inversion of each subject, fusing across all modalities
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_invert.m 3544 2009-11-09 13:45:56Z vladimir $
 
% check whether this is a group inversion
%--------------------------------------------------------------------------
if ~iscell(D), D = {D}; end
Nl = length(D);                                  % number of subjects

% D - SPM data structure
%==========================================================================
if nargin == 1, val = 1; end
inverse    = D{1}.inv{D{1}.val}.inverse;
 
% defaults
%--------------------------------------------------------------------------
try, type  = inverse.type;   catch, type  = 'GS';              end
try, s     = inverse.smooth; catch, s     = 0.6;               end
try, Np    = inverse.Np;     catch, Np    = 256;               end
try, Nm    = inverse.Nm;     catch, Nm    = 128;               end
try, Nr    = inverse.Nr;     catch, Nr    = 16;                end
try, xyz   = inverse.xyz;    catch, xyz   = [0 0 0];           end
try, rad   = inverse.rad;    catch, rad   = 128;               end
try, lpf   = inverse.lpf;    catch, lpf   = 0;                 end
try, hpf   = inverse.hpf;    catch, hpf   = 48;                end
try, sdv   = inverse.sdv;    catch, sdv   = 4;                 end
try, Han   = inverse.Han;    catch, Han   = 1;                 end
try, Na    = inverse.Na;     catch, Na    = 1024;              end
try, woi   = inverse.woi;    catch, woi   = [];                end
try, pQ    = inverse.pQ;     catch, pQ    = [];                end
 
try
    modalities        = inverse.modality;  
catch
    [tmp, modalities] = modality(D{1}, 1, 1);
end

if ischar(modalities)
    modalities = {modalities};
end

Nmod        = numel(modalities);             % number of modalities
Nmax        = Nm;                            % max number of spatial modes
Nrmax       = Nr;                            % max number of temporal modes

fprintf('Checking all leadfields:\n')
for i = 1:Nl    
    for m = 1:Nmod
    
        fprintf('Subject %d: %s\n',i, modalities{m});
        
        % Check gain or lead-field matrices
        %--------------------------------------------------------------------------
        [L D{i}] = spm_eeg_lgainmat(D{i});
        Ic{i}{m} = setdiff(meegchannels(D{i}, modalities{m}), badchannels(D{i}));
        Nd(i)    = size(L,2);                          % number of dipoles
        Nc(i,m)  = length(Ic{i}{m});                   % number of channels
        
        if isempty(Ic{i}{m})
            error(['The specified modality is missing from file ' D{i}.fname]);
        end
        
        if any(diff(Nd))
            warndlg('please ensure subjects have the same number of dipoles')
            return
        end
    
        % Check for null space over sensors and remove it
        %--------------------------------------------------------------------------
        try
            SX      = D{i}.sconfounds{m};
            R{i}{m} = speye(Nc(i,m),Nc(i,m)) - SX*pinv(SX);
        catch
            R{i}{m} = speye(Nc(i,m),Nc(i,m));
        end
    end
end
fprintf(' - done\n')
    
% check restriction; assume radii are the same for all VOI
%--------------------------------------------------------------------------
Nd    = Nd(1);                                     % number of dipoles
Nv    = size(xyz,1);                               % number of VOI
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end
    
    
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================

fprintf('Computing Green function from graph Laplacian:')
%--------------------------------------------------------------------------
vert  = D{1}.inv{D{1}.val}.mesh.tess_mni.vert;
face  = D{1}.inv{D{1}.val}.mesh.tess_mni.face;
A     = spm_eeg_inv_meshdist(vert,face,0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
QG    = QG.*(QG > exp(-8));
QG    = QG*QG;
clear Qi A GL
fprintf(' - done\n')


% Restrict source space by eliminating dipoles
%--------------------------------------------------------------------------
Is    = sparse(Nd,1);
for i = 1:Nv
    Iv = sum([vert(:,1) - xyz(i,1), ...
        vert(:,2) - xyz(i,2), ...
        vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
    Is = Is | Iv;
end
Is    = find(Is);
vert  = vert(Is,:);
QG    = QG(Is,Is);
Ns    = length(Is);

% check for (eg, empty-room) sensor components (in Qe{1})
%==========================================================================
for i = 1:Nl    
    for m = 1:Nmod
        try
            pQe      = D{i}.inv{D{i}.val}.inverse.Qe{m};
            QE{i}{m} = pQe/sqrt(pQe'*pQe);
            if length(QE{i}{m}) ~= Nc(i,m)
                error('Specified error component for modality %s of subject %d does not match number of channels (%d)\n',modalities{m},i,Nc(i,m))
            end
        catch
            pQe = [];
            QE{i}{m} = 1;
        end  
    end
    
    if length(pQe)
        fprintf('Using sensor error component provided for modality %s...\n',modalities{m});
    end
end
  
%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================

for m = 1:Nmod
    
    % Project to channel modes (U); Checking for null space over sensors
    %--------------------------------------------------------------------------
    L       = R{1}{m}*spm_eeg_lgainmat(D{1},Is,D{1}.chanlabels(Ic{1}{m}));
    U{1}{m} = spm_svd(L*L',exp(-16));
    Nm(1,m) = min(size(U{1}{m},2),Nmax);
    U{1}{m} = U{1}{m}(:,1:Nm(1,m));
    A{1,m}  = U{1}{m}';
    UL{m}   = U{1}{m}'*L;
    
    % Scale projected leadfield (in case fusing multiple modalities)
    %--------------------------------------------------------------------------
    Lscale(m) = 1/sqrt(trace(UL{m}*UL{m}')/Nm(1,m));
    UL{m}     = UL{m}*Lscale(m);
    
    fprintf('Using %d spatial modes for modality %s subject %d\n', Nm(1,m), modalities{m}, 1);
    
    % Spatial projectors for subsequent subjects
    %--------------------------------------------------------------------------
    for i = 2:Nl
        L       = R{i}{m}*spm_eeg_lgainmat(D{i},Is,D{i}.chanlabels(Ic{i}{m}));
        U{i}{m} = spm_svd(L*L',exp(-16));
        Nm(i,m) = min(size(U{i}{m},2),Nmax);
        U{i}{m} = U{i}{m}(:,1:Nm(i,m));
        A{i,m}  = UL{m}*pinv(full(L));
        fprintf('Using %d spatial modes for modality %s subject %d\n',Nm(i,m),modalities{m}, i);
    end       
end

    
%==========================================================================
% Temporal projector
%==========================================================================

Nn    = sparse(1,Nl);                       % total number of samples
AY    = {};                                 % aligned data for MVB
AYYA  = sparse(0);                          % aligned covariance for ReML

for m = 1:Nmod                              % error covariance per modality
     Qe{m} = sparse(0);
end
AQe{1} = sparse(0);                         % aligned error covariance


for i = 1:Nl
                       
    UYYU{i}  = sparse(0);                          % covariance for ReML

    % Time-window of interest
    %----------------------------------------------------------------------
    if isempty(woi)
        w{i} = 1000*[min(D{i}.time) max(D{i}.time)];
    else
        w{i} = woi;
    end
    It{i}  = (w{i}/1000 - D{i}.timeonset)*D{i}.fsample + 1;
    It{i}  = max(1,It{i}(1)):min(It{i}(end), length(D{i}.time));
    It{i}  = fix(It{i});
    
    % Peristimulus time
    %----------------------------------------------------------------------
    pst{i} = 1000*D{i}.time;
    pst{i} = pst{i}(It{i});                       % peristimulus time (ms)
    dur    = (pst{i}(end) - pst{i}(1))/1000;      % duration (s)
    dct{i} = (It{i} - It{i}(1))/2/dur;            % DCT frequencies (Hz)
    Nb(i)  = length(It{i});                       % number of time bins
    
    % Serial correlations
    %----------------------------------------------------------------------
    K      = exp(-(pst{i} - pst{i}(1)).^2/(2*sdv^2));
    K      = toeplitz(K);
    qV{i}  = sparse(K*K');
    
    % Confounds and temporal subspace
    %----------------------------------------------------------------------
    T      = spm_dctmtx(Nb(i),Nb(i));
    j      = find( (dct{i} >= lpf) & (dct{i} <= hpf) );
    T      = T(:,j);
    dct{i} = dct{i}(j);
    

    % Hanning operator (if requested)
    %----------------------------------------------------------------------
    if Han
        W   = sparse(1:Nb(i),1:Nb(i),spm_hanning(Nb(i)));
    else
        W   = 1;
    end
        
    % get trials
    %----------------------------------------------------------------------
    try
        trial = D{i}.inv{D{i}.val}.inverse.trials;
    catch
        trial = unique(D{i}.conditions);
    end
    Nt(i) = length(trial);
       
    % get temporal covariance (Y'*Y) to find temporal modes
    %======================================================================
        
    for m = 1:Nmod       
                   
        YTY   = sparse(0);
       
        % get data 
        %----------------------------------------------------------------------
        for j = 1:Nt(i)
            c     = D{i}.pickconditions(trial{j});
            for k = 1:length(c)  
                Yk     = U{i}{m}'*R{i}{m}*squeeze(D{i}(Ic{i}{m},It{i},c(k)));
                YTY    = YTY + Yk'*Yk;
            end
        end     
     
        % Scale data (in case fusing multiple modalities)
        %--------------------------------------------------------------------------
        Yscale(i,m) = 1/sqrt(trace(YTY)/Nm(i,m));
        YTY         = YTY*(Yscale(i,m)^2);

        % temporal projector (at most Nrmax modes) S = T*V
        %======================================================================
       
        YTY     = W'*YTY*W;                           % Apply any Hanning
        YTY     = T'*YTY*T;                           % Apply filtering
        [v E]   = spm_svd(YTY,exp(-8));
        E       = diag(E)/trace(YTY);
        Nr(i,m) = min(length(E),Nrmax);               % number of temporal modes
        V{i}{m} = v(:,1:Nr(i,m));                     % temporal modes
        VE(i,m) = sum(E(1:Nr(i,m)));                  % variance explained
        
     end
                 
    % combine temporal modes across modalities
    %--------------------------------------------------------------------------
    fV{i}     = full(spm_cat(V{i}));
    [fV{i} E] = spm_svd(fV{i}*fV{i}',-Inf);
    fNr(i)    = min(length(E),sum(Nr(i,:),2));
    fV{i}     = fV{i}(:,1:fNr(i));
    E         = diag(E);
    E         = E/sum(E);
    fVE(i)    = mean(VE(i,:),2)*full(sum(E(1:fNr(i))));

    fprintf('Using %i temporal modes for subject %i, ',fNr(i))
    fprintf('accounting for %0.2f percent average variance\n',full(100*fVE(i)))
    
    % projection and whitening
    %--------------------------------------------------------------------------
    S{i}  = T*fV{i};                                    % temporal projector
    qP    = inv(S{i}'*qV{i}*S{i});                   % precision (mode)
    Vq{i} = S{i}*qP*S{i}';                           % precision (time)
     
    
    % Get all data...
    %======================================================================  

    % get temporal covariance (Y*Y') for Gaussian process model
    %----------------------------------------------------------------------
    for j = 1:Nt(i)
                    
        UY{i,j} = sparse(0);
         
        c     = D{i}.pickconditions(trial{j});
        for k = 1:length(c)
            
            Uy    = [];
            Ay    = [];
            
            for m = 1:Nmod       
        
                % concatenate trial-specific data across modalities
                %--------------------------------------------------------------
                y        = squeeze(D{i}(Ic{i}{m},It{i},c(k)))*S{i};
                y        = y * Yscale(i,m);
                Uy       = [Uy; U{i}{m}'*y];
                Ay       = [Ay; A{i,m}  *y];
            end
            
            % accumulate first and second-order statistics
            %--------------------------------------------------------------
            UY{i,j}    = UY{i,j} + Uy;
            UYYU{i}    = UYYU{i} + Uy*Uy';
           
            % and repeat for aligned data
            %--------------------------------------------------------------
            AY{end + 1} = Ay;
            AYYA        = AYYA  + Ay*Ay';
                        
            Nn(i)       = Nn(i) + fNr(i);

        end
    end
   
    % generate sensor component (Qe) per modality; 
    %----------------------------------------------------------------------
    for m = 1:Nmod
        Qn{m,m}    = sparse(Nm(i,m),Nm(i,m));
    end
    for m = 1:Nmod
        Q         = Qn;
        Q{m,m}    = U{i}{m}'*QE{i}{m}*U{i}{m};
        Qe{m}     = spm_cat(Q);

        % assuming equal noise over subjects and modalities
        %----------------------------------------------------------------------
        Q{m,m}    = A{i,m}*QE{i}{m}*A{i,m}';
        AQe{1}    = AQe{1} + spm_cat(Q);
    end
end
         
% concatenate leadfields and aligned data 
%----------------------------------------------------------------------
UL       = spm_cat(UL(:));
AY       = spm_cat(AY);
    
Nr = fNr;
VE = fVE;
V  = fV;

%==========================================================================
% Step 1: Optimise spatial priors over subjects
%==========================================================================
    
    
% create source components (Qp)
%==========================================================================
switch(type)
    
    case {'MSP','GS','ARD'}
        
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        Ip    = ceil([1:Np]*Ns/Np);
        for i = 1:Np
            
            % left hemisphere
            %--------------------------------------------------------------
            q               = QG(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
            
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                vert(:,2) - vert(Ip(i),2), ...
                vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
            
            % bilateral
            %--------------------------------------------------------------
            q               = QG(:,Ip(i)) + QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
            
        end
        
    case {'LOR','COH'}
        
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
        
        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = UL*Qp{2}*UL';
        
        
    case {'IID','MMN'}
        
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
        
end

% augment with exogenous (eg, fMRI) source priors in pQ
%==========================================================================
for i = 1:length(pQ)
    
    switch(type)
        
        case {'MSP','GS','ARD'}
            %------------------------------------------------------------------
            if isvector(pQ{i}) && length(pQ{i}) == Ns
                Qp{end + 1}.q   = pQ{i}(:);
                LQpL{end + 1}.q = UL*Qp{end}.q;
            else
                error('Using MSP(GS/ARD) please supply spatial priors as vectors')
            end
            
        case {'LOR','COH','IID','MMN'}
            %------------------------------------------------------------------
            if isvector(pQ{i}) && length(pQ{i}) == Ns
                pQ{i}         = pQ{i}(:);     
                Qp{end + 1}   = sparse(diag(pQ{i}*pQ{i}'));
                LQpL{end + 1} = UL*Qp{end}*UL';
            elseif size(pQ{i},1) == Ns && size(pQ{i},2) == Ns
                Qp{end + 1}   = pQ{i};
                LQpL{end + 1} = UL*Qp{end}*UL';
            else
                error('spatial priors are the wrong size')
            end
    end
end

if length(pQ)
    fprintf('Using %d spatial source priors provided...\n',length(pQ));
end

% Inverse solution
%==========================================================================

QP = {};

% get source-level priors (using all subjects)
%--------------------------------------------------------------------------
switch(type)
    
    case {'MSP','GS'}
        
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Npm   = length(Qp);
        Q     = sparse(Ns,Npm);
        for i = 1:Npm
            Q(:,i) = Qp{i}.q;
        end
        
        % Multivariate Bayes
        %------------------------------------------------------------------
        MVB   = spm_mvb(AY,UL,[],Q,AQe,16);
        
        % Spatial priors (QP); eliminating minor patterns
        %------------------------------------------------------------------
        cp    = diag(MVB.cp);
        [i j] = sort(-diag(MVB.cp));
        try
            j = j(1:256);
        end
        qp    = Q(:,j)*MVB.cp(j,j)*Q(:,j)';
        
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;
        
end

switch(type)
    
    case {'MSP','ARD'}
        
        % or ReML - ARD
        %------------------------------------------------------------------
        qp          = sparse(0);
        Q           = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F] = spm_sp_reml(AYYA,[],Q,sum(Nn));
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Npm   = length(Qp);
        hp    = h([1:Npm] + Ne);
        for i = 1:Npm
            if hp(i) > max(hp)/128;
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
        
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;
        
        
    case {'IID','MMN','LOR','COH'}
        
        % or ReML - ARD
        %------------------------------------------------------------------
        qp          = sparse(0);
        Q           = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F] = spm_reml_sc(AYYA,[],Q,sum(Nn));
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Npm   = length(Qp);
        hp    = h([1:Npm] + Ne);
        for i = 1:Npm
            qp = qp + hp(i)*Qp{i};
        end
        
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;
        
end

 
%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================

Lscale = [];

for i = 1:Nl
 
    fprintf('Inverting subject %i\n',i)
    
    % Concatenate leadfields over modalities
    %--------------------------------------------------------------------------
    
    L = {};
    for m    = 1:Nmod
        L{m}        = U{i}{m}'*R{i}{m}* ...
                      spm_eeg_lgainmat(D{1},Is,D{1}.chanlabels(Ic{i}{m}));
        Lscale(i,m) = 1/sqrt(trace(L{m}*L{m}')/Nm(i,m));
        L{m}        = L{m} * Lscale(i,m);
    end
    L = spm_cat(L(:));
           
    % generate sensor component (Qe) per modality
    %----------------------------------------------------------------------
    Qe = {};
    for m = 1:Nmod
        Q         = Qn;
        Q{m,m}    = U{i}{m}'*QE{i}{m}*U{i}{m};
        Qe{m}     = spm_cat(Q);
    end  

    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    Ne    = length(Qe);
    Np    = length(QP);
    LQpL  = {};
    for j = 1:Np
        LQpL{j}  = L*QP{j}*L';
    end
    Q     = {Qe{:} LQpL{:}};
 
    % re-do ReML
    %----------------------------------------------------------------------
    [Cy,h,Ph,F] = spm_reml_sc(UYYU{i},[],Q,Nn(i));
    
    % Data ID
    %==========================================================================  
    ID  = spm_data_id(UYYU{i});
 
    % Covariances: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    Qp    = sparse(0);
    hp    = h([1:Np] + Ne);
    for j = 1:Np
        Qp = Qp + hp(j)*QP{j};
    end
    LCp   = L*Qp;
 
    % MAP estimates of instantaneous sources
    %======================================================================
    iC    = inv(Cy);
    M     = LCp'*iC;
 
    % conditional covariance (leading diagonal)
    % Cq    = Cp - Cp*L'*iC*L*Cp;
    %----------------------------------------------------------------------
    Cq    = diag(Qp) - sum(LCp.*M')';
 
    % evaluate conditional expectation (of the sum over trials)
    %----------------------------------------------------------------------
    SSR   = 0;
    SST   = 0;
    J     = {};
    for j = 1:Nt(i)
        
        % trial-type specific source reconstruction
        %------------------------------------------------------------------
        J{j} = M*UY{i,j};
 
        % sum of squares
        %------------------------------------------------------------------
        SSR  = SSR + sum(var((UY{i,j} - L*J{j}),0,2));
        SST  = SST + sum(var( UY{i,j},0,2));
 
    end
 
    % accuracy; signal to noise (over sources)
    %======================================================================
    R2   = 100*(SST - SSR)/SST;
    fprintf('Percent variance explained %.2f (%.2f)\n',R2,R2*VE(i))
    
    % Save results
    %======================================================================
    inverse.type   = type;                 % inverse model
    inverse.smooth = s;                    % smoothness (0 - 1)
    inverse.xyz    = xyz;                  % VOI (XYZ)
    inverse.rad    = rad;                  % VOI (radius)
    inverse.Lscale = Lscale(i,:);          % leadfield scale-factor
    inverse.scale  = Yscale(i,:);          % data scale-factor
    inverse.M      = M;                    % MAP projector (reduced)
    inverse.J      = J;                    % Conditional expectation
    inverse.Y      = UY(i,:);              % ERP data (reduced)
    inverse.L      = L;                    % Lead-field (reduced)    
    inverse.qC     = Cq;                   % spatial covariance
    inverse.qV     = Vq{i};                % temporal correlations
    inverse.T      = S{i};                 % temporal subspace  
    inverse.U      = spm_cat(diag(U{i}));  % spatial subspace
    inverse.R      = spm_cat(diag(R{i}));  % Re-referencing matrix
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It{i};                % Indices of time bins
    inverse.Ic     = Ic{i};                % Indices of good channels
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.pst    = pst{i};               % peristimulus time
    inverse.dct    = dct{i};               % frequency range
    inverse.F      = F;                    % log-evidence
    inverse.ID     = ID;                   % data ID
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.VE     = VE(i);                % variance explained
    inverse.woi    = w{i};                 % time-window inverted
    
    inverse.modality = modalities;         % modalities inverted
    
    % save in struct
    %----------------------------------------------------------------------
    D{i}.inv{D{i}.val}.inverse = inverse;
    D{i}.inv{D{i}.val}.method  = 'Imaging';
    
    % and delete old contrasts
    %----------------------------------------------------------------------
    try
        D{i}.inv{D{i}.val} = rmfield(D{i}.inv{D{i}.val},'contrast');
    end
 
    % display
    %======================================================================
    spm_eeg_invert_display(D{i});
    drawnow
 
end
 
if length(D) == 1, D = D{1}; end
return
