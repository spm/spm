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
%     inverse.J{i}   - Conditional expectation (i conditions) J = M*U*Y
%     inverse.L      - Lead field (reduced UL := U*L)
%     inverse.qC     - spatial covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal projector
%     inverse.U(j)   - spatial projector (j modalities)
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
% 1. This routine implements "group-based" inversion, corresponding to
% ill-posed linear models of the following form:
%
% [AY{1}...AY{n}] = L(1} * [J{1}...J{n}]   +  [e{1}...e{n}]
%
% where AY{i} are the spatially normalized or adjusted data from subject i
% that would have been seen if the lead-field L{i} = L{1}. The ensuing
% Gaussian process priors on sources are then used to estimate subject-
% specific MAP estimates of J{i} using
%
% AY{i} = L(1} * J{i}  +  e{i}
%
% using spatial priors from the group model above.
%
% Here, A{i}  = L{1}*pinv(L{i}) =>
%       AY{i} = A(i}*L(i}*J{i}
%             = L(1}*J{i}
%
% Potential scaling differences between the lead-fields are handled by
% scaling L{1} such that trace(L{1}*L{1}') = constant (number of spatial 
% modes or channels), while scaling the data such that trace(AY{n}*AY{n}') =
% constant over subjects (and modalities; see below).
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
%             AY{1}{1,...,t}  = L(1} * J{1,...,t}   +  e{{1,...,t}}
%             AY{2}{1,...,t}  = L(2}                   e{{2,...,t}}
%                  .
%                  .
%                  .
%             AY{m}{1,...,t}  = L(n}                   e{{n,...,t}}
%
% Under empirical priors on J{1,...,t} for m modalities with t trial types.
%
% See: MEG and EEG data fusion: Simultaneous localisation of face-evoked
% responses.
% Henson R, Mouchlianitis E & Friston K.
% Neuroimage. 2009. 47:581-9.
%__________________________________________________________________________
%
% 3. It also allows incorporation of spatial source priors, eg, from fMRI
% (see spm_eeg_inv_fmripriors.m). Note that if a vector is passed in 
% inverse.pQ, then variance components used (pass a matrix if a covariance
% component is desired).
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
% $Id: spm_eeg_invert.m 3633 2009-12-11 10:34:15Z guillaume $
 
% check whether this is a group inversion
%--------------------------------------------------------------------------
if ~iscell(D), D = {D}; end
Nl = length(D);                                  % number of subjects
 
 
% D - SPM data structure
%==========================================================================
if nargin > 1
    D{1}.val = val;
elseif ~isfield(D{1}, 'val')
    D{1}.val = 1;
end

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
 
 
% get specified modalities to invert (default to all)
%--------------------------------------------------------------------------
try
    modalities     = inverse.modality;
    if ~iscell(modalities)
        modalities = {modalities};
    end
catch
    for m = 1:length(D{1}.inv{D{1}.val}.forward)
        modalities{m} = D{1}.inv{D{1}.val}.forward(m).modality;
    end
end
Nmod      = numel(modalities);               % number of modalities
Nmmax     = Nm;                              % max number of spatial modes
Nrmax     = Nr;                              % max number of temporal modes

 
%==========================================================================
for i = 1:Nl
    fprintf('Checking leadfields for subject %i\n',i)
    [L D{i}] = spm_eeg_lgainmat(D{i});
    
    for m = 1:Nmod
        
        % Check gain or lead-field matrices
        %------------------------------------------------------------------
        Ic{i,m}  = setdiff(meegchannels(D{i}, modalities{m}), badchannels(D{i}));
        Nd(i)    = size(L,2);                          % number of dipoles
        Nc(i,m)  = length(Ic{i,m});                    % number of channels
 
        if isempty(Ic{i,m})
            errordlg(['Modality ' modalities{m} 'is missing from file ' D{i}.fname]);
            return
        end
 
        if any(diff(Nd))
            errordlg('Please ensure subjects have the same number of dipoles')
            return
        end
 
        % Check for null space over sensors and remove it
        %------------------------------------------------------------------
        try
            SX     = D{i}.sconfounds{m};
            R{i,m} = speye(Nc(i,m),Nc(i,m)) - SX*pinv(SX);
        catch
            R{i,m} = speye(Nc(i,m),Nc(i,m));
        end
    end
end
fprintf(' - done\n')
 
 
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================
 
 
fprintf('Computing Green function from graph Laplacian:')
%--------------------------------------------------------------------------
Nd    = Nd(1);                                     % number of dipoles
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
 
 
% check restriction (assume radii are the same for all VOI)
%==========================================================================
Nv    = size(xyz,1);                               % number of VOI
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end
 
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
Ns    = length(Is);                               % number of sources
 
 
% check for (e.g., empty-room) sensor components (in Qe{1})
%==========================================================================
for i = 1:Nl
    for m = 1:Nmod
        try
            Qe      = D{i}.inv{D{i}.val}.inverse.Qe{m};
            QE{i,m} = Nc(i,m)*Qe/trace(Qe);
            if length(QE{i,m}) ~= Nc(i,m)
                errordlg('error component (modality %s; subject %d) does not match number of channels (%d)\n',modalities{m},i,Nc(i,m))
                return
            end
            
        % assume i.i.d. if not specified
        %------------------------------------------------------------------
        catch
            Qe      = [];
            QE{i,m} = 1;
        end
    end
end
if ~isempty(Qe)
    fprintf('Using sensor error component provided...\n');
end
 
    
%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================
 
for m = 1:Nmod
 
    % Project to channel modes (U := A{1,m}) (with null space over sensors)
    %----------------------------------------------------------------------
    L      = R{1,m}*spm_eeg_lgainmat(D{1},Is,D{1}.chanlabels(Ic{1,m}));
    U      = spm_svd(L*L',exp(-16));
    Nm(m)  = min(size(U,2),Nmmax);
    A{1,m} = U(:,1:Nm(m))';
    UL{m}  = A{1,m}*L;
 
    % Scale projected lead-field
    %----------------------------------------------------------------------
    UL{m}     = UL{m}/sqrt(trace(UL{m}*UL{m}')/Nm(m));
 
    % Spatial projectors A{i,m) for subsequent subjects (i)
    %----------------------------------------------------------------------
    for i = 2:Nl
        L      = R{i,m}*spm_eeg_lgainmat(D{i},Is,D{i}.chanlabels(Ic{i,m}));
        A{i,m} = UL{m}*pinv(full(L));
    end
 
    % Report
    %----------------------------------------------------------------------
    fprintf('Using %d spatial modes for modality %s\n',Nm(m),modalities{m})
end
 
 
 
%==========================================================================
% Temporal projector
%==========================================================================
 
Nn    = sparse(1,Nl);                      % total number of samples
AY    = {};                                % group response for MVB
AYYA  = sparse(0);                         % group response for ReML
 
% loop over Nl lead-fields (subjects)
%----------------------------------------------------------------------
for i = 1:Nl
 
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
    pst{i} = 1000*D{i}.time;                      % peristimulus time (ms)
    pst{i} = pst{i}(It{i});                       % windowed time (ms)
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
 
    % get trials or conditions
    %----------------------------------------------------------------------
    try
        trial = D{i}.inv{D{i}.val}.inverse.trials;
    catch
        trial = unique(D{i}.conditions);
    end
    Nt(i) = length(trial);
    
 
    % get temporal covariance (Y'*Y) to find temporal modes
    %======================================================================
    Ay    = cell(Nmod,1);                            % modality-specifc
    YTY   = sparse(0);                               % accumulator
    for m = 1:Nmod                                   % loop over modalities
 
        % get (spatially aligned) data
        %------------------------------------------------------------------
        YY    = 0;
        Ay{m} = 0;
        for j = 1:Nt(i)                              % pool over conditions
            c     = D{i}.pickconditions(trial{j});   % and trials
            for k = 1:length(c)
                y     = A{i,m}*D{i}(Ic{i,m},It{i},c(k));
                Ay{m} = Ay{m} + y;
                YY    = YY + y'*y;
            end
        end
        
        % Scale data (in case fusing multiple modalities)
        %------------------------------------------------------------------
        scale(i,m) = sign(trace(Ay{m}'*(UL{m}*UL{1}')*Ay{1}));
        scale(i,m) = scale(i,m)/sqrt(trace(YY)/Nm(m));
        YTY        = YTY + YY*(scale(i,m)^2);
        
    end
 
    % temporal projector (at most Nrmax modes) S = T*V
    %======================================================================
    YTY     = W'*YTY*W;                           % Apply any Hanning
    YTY     = T'*YTY*T;                           % Apply filtering
    [VT E]  = spm_svd(YTY,exp(-8));               % get temporal modes
    E       = diag(E)/trace(YTY);                 % normalise variance
    Nr(i)   = min(length(E),Nrmax);               % number of temporal modes
    V{i}    = VT(:,1:Nr(i));                      % temporal modes
    VE(i)   = sum(E(1:Nr(i)));                    % variance explained
 
    fprintf('Using %i temporal modes for subject %i, ',Nr(i),i)
    fprintf('accounting for %0.2f percent average variance\n',full(100*VE(i)))
 
    % projection and whitening
    %----------------------------------------------------------------------
    S{i}  = T*V{i};                               % temporal projector
    qP    = inv(S{i}'*qV{i}*S{i});                % precision (mode)
    Vq{i} = S{i}*qP*S{i}';                        % precision (time)
 
 
    % get spatial covariance (Y*Y') for Gaussian process model.
    %======================================================================
    UYYU{i}  = sparse(0);
 
    % loop over Nt trial types
    %----------------------------------------------------------------------
    for j = 1:Nt(i)
 
        UY{i,j} = sparse(0);
        c     = D{i}.pickconditions(trial{j});
        for k = 1:length(c)
 
            % stack (scaled aligned data) over modalities
            %--------------------------------------------------------------
            Ay    = cell(Nmod,1);                % modality accumulator
            for m = 1:Nmod
                y       = D{i}(Ic{i,m},It{i},c(k))*S{i};
                Ay{m}   = A{i,m}*y*scale(i,m);
            end
            Ay          = spm_cat(Ay)/length(c); % contribution to ERP
            AyyA        = Ay*Ay';                % and covariance
            Nn(i)       = Nn(i) + Nr(i);         % number of samples
 
            % accumulate first & second-order statistics (subject-specific)
            %--------------------------------------------------------------
            UY{i,j}     = UY{i,j} + Ay;          % condition-specific ERP 
            UYYU{i}     = UYYU{i} + AyyA;        % subject-specific covariance
 
            % and pool for optimisation of spatial priors over subjects
            %--------------------------------------------------------------
            AY{end + 1} = Ay;                    % pooled response for MVB
            AYYA        = AYYA  + AyyA;          % pooled response for ReML
 
        end
    end
 
end
 
% concatenate responses over realisations and lead fields over modalities
%--------------------------------------------------------------------------
AY   = spm_cat(AY);
UL   = spm_cat(UL(:));
 
 
% generate sensor error components (Qe)
%==========================================================================
AQe   = {sparse(0)};
for m = 1:Nmod
    Qe{m} = sparse(0);
end
 
% assuming equal noise over subjects (Qe{m}) and modalities AQe
%----------------------------------------------------------------------
for i = 1:Nl
    for m = 1:Nmod
        N{m,m}    = sparse(Nm(m),Nm(m));
    end
    for m = 1:Nmod
        Q         = N;
        Q{m,m}    = A{i,m}*QE{i,m}*A{i,m}';
        Qe{m}     = Qe{m}  + spm_cat(Q);
        AQe{1}    = AQe{1} + Qe{m};
    end
end
 
 
 
 
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
 
 
% augment with exogenous (e.g., fMRI) source priors in pQ
%==========================================================================
for i = 1:length(pQ)
 
    switch(type)
 
        case {'MSP','GS','ARD'}
            %--------------------------------------------------------------
            if isvector(pQ{i}) && length(pQ{i}) == Ns
                Qp{end + 1}.q   = pQ{i}(:);
                LQpL{end + 1}.q = UL*Qp{end}.q;
            else
                errordlg('Using MSP(GS/ARD) please supply spatial priors as vectors')
                return
            end
 
        case {'LOR','COH','IID','MMN'}
            %--------------------------------------------------------------
            if isvector(pQ{i}) && length(pQ{i}) == Ns
                pQ{i}         = pQ{i}(:);
                Qp{end + 1}   = sparse(diag(pQ{i}.^2));
                LQpL{end + 1} = UL*Qp{end}*UL';
            elseif size(pQ{i},1) == Ns && size(pQ{i},2) == Ns
                Qp{end + 1}   = pQ{i};
                LQpL{end + 1} = UL*Qp{end}*UL';
            else
                errordlg('spatial priors are the wrong size')
                return
            end
    end
end
if ~isempty(pQ)
    fprintf('Using %d spatial source priors provided...\n',length(pQ));
end
 
 
% Inverse solution
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};
 
% Get source-level priors (using all subjects)
%--------------------------------------------------------------------------
switch(type)
 
    case {'MSP','GS'}
 
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = sparse(Ns,Np);
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
 
        % Multivariate Bayes
        %------------------------------------------------------------------
        MVB   = spm_mvb(AY,UL,[],Q,AQe,16);
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
 
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
        Np    = length(Qp);
        hp    = h([1:Np] + Ne);
        for i = 1:Np
            if hp(i) > max(hp)/128;
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
 
end
 
switch(type)
 
    case {'IID','MMN','LOR','COH'}
 
        % or ReML - ARD
        %------------------------------------------------------------------
        qp          = sparse(0);
        Q           = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F] = spm_reml_sc(AYYA,[],Q,sum(Nn));
 
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h([1:Np] + Ne);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
 
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
 
end
 
 
 
%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================
 
for i = 1:Nl
 
    fprintf('Inverting subject %i\n',i)
    
    % generate sensor component (Qe) per modality
    %----------------------------------------------------------------------
    Qe    = {};
    for m = 1:Nmod
        N{m,m}    = sparse(Nm(m),Nm(m));
    end
    for m = 1:Nmod
        Q         = N;
        Q{m,m}    = A{i,m}*QE{i,m}*A{i,m}';
        Qe{m}     = spm_cat(Q);
    end
    
    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    Np    = length(LQPL);
    Ne    = length(Qe);
    Q     = {Qe{:} LQPL{:}};
 
    % re-do ReML
    %----------------------------------------------------------------------
    [Cy,h,Ph,F] = spm_reml_sc(UYYU{i},[],Q,Nn(i));
 
    % Data ID
    %==========================================================================
    ID  = spm_data_id(UYYU{i});
 
 
    % Covariance: sensor space - Ce and source space - L*Cp
    %--------------------------------------------------------------------------
    Cp    = sparse(0);
    LCp   = sparse(0);
    hp    = h([1:Np] + Ne);
    for j = 1:Np
        Cp  =  Cp + hp(j)*QP{j};
        LCp = LCp + hp(j)*LQP{j};
    end
 
    % MAP estimates of instantaneous sources
    %======================================================================
    iC    = inv(Cy);
    M     = LCp'*iC;
 
    % conditional variance (leading diagonal)
    % Cq    = Cp - Cp*L'*iC*L*Cp;
    %----------------------------------------------------------------------
    Cq    = Cp - sum(LCp.*M')';
 
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
        SSR  = SSR + sum(var((UY{i,j} - UL*J{j}),0,2));
        SST  = SST + sum(var( UY{i,j},0,2));
 
    end
 
    % accuracy; signal to noise (over sources)
    %======================================================================
    R2   = 100*(SST - SSR)/SST;
    fprintf('Percent variance explained %.2f (%.2f)\n',R2,R2*VE(i))
 
    % Save results (for first modality)
    %======================================================================
    inverse.type   = type;                 % inverse model
    inverse.smooth = s;                    % smoothness (0 - 1)
    inverse.xyz    = xyz;                  % VOI (XYZ)
    inverse.rad    = rad;                  % VOI (radius)
    inverse.scale  = scale(i,:);           % data scale-factor
    inverse.M      = M;                    % MAP projector (reduced)
    inverse.J      = J;                    % Conditional expectation
    inverse.Y      = UY(i,:);              % ERP data (reduced)
    inverse.L      = UL;                   % Lead-field (reduced)
    inverse.qC     = Cq;                   % spatial covariance
    inverse.qV     = Vq{i};                % temporal correlations
    inverse.T      = S{i};                 % temporal projector
    inverse.U      = A(i,:);               % spatial projector
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It{i};                % Indices of time bins
    inverse.Ic     = Ic(i,:);              % Indices of good channels
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
