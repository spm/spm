function [D] = spm_eeg_invert(D, val)
% ReML inversion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert(D, val)
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
%     inverse.pQ     - any source priors (e.g. from fMRI); vector or matrix
%     inverse.Qe     - any sensor error components (e.g. empty-room data)
%     inverse.dplot  - make diagnostics plots (0 or 1)
%     inverse.STAT   - flag for stationarity assumption, which invokes a 
%                      full DCT temporal projector (from lpf to hpf Hz)
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
% Human Brain Mapping. 2010. 1(10):1512-31.
%__________________________________________________________________________
%
% The routine essentially consists of two steps:
%
%   1. Optimisation of spatial source priors over subjects
%   2. Re-inversion of each subject, fusing across all modalities
%__________________________________________________________________________
% Copyright (C) 2006-2017 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_invert.m 7118 2017-06-20 10:33:27Z guillaume $


SVNid = '$Rev: 7118 $';

%-Say hello
%--------------------------------------------------------------------------
spm('FnBanner',mfilename,SVNid);

% check whether this is a group inversion for (Nl) number of subjects
%--------------------------------------------------------------------------
if ~iscell(D), D = {D}; end
Nl = length(D);
 
 
% D - SPM data structure
%==========================================================================
if nargin > 1
    D{1}.val = val;
elseif ~isfield(D{1}, 'val')
    D{1}.val = 1;
end
for i = 2:Nl
    D{i}.val = D{1}.val;
end
 
% forward model
%--------------------------------------------------------------------------
inverse   = D{1}.inv{D{1}.val}.inverse;
 
% defaults
%--------------------------------------------------------------------------
try, STAT = inverse.STAT;   catch, STAT = 0;        end
try, type = inverse.type;   catch, type = 'GS';     end
try, s    = inverse.smooth; catch, s    = 0.6;      end
try, Np   = inverse.Np;     catch, Np   = 256;      end
try, Nr   = inverse.Nr;     catch, Nr   = 16;       end
try, xyz  = inverse.xyz;    catch, xyz  = [0 0 0];  end
try, rad  = inverse.rad;    catch, rad  = 128;      end
try, mask = inverse.mask;   catch, mask = [];       end
try, lpf  = inverse.lpf;    catch, lpf  = 0;        end
try, hpf  = inverse.hpf;    catch, hpf  = 48;       end
try, sdv  = inverse.sdv;    catch, sdv  = 4;        end
try, Han  = inverse.Han;    catch, Han  = 1;        end
try, woi  = inverse.woi;    catch, woi  = [];       end
try, pQ   = inverse.pQ;     catch, pQ   = [];       end
try, dp   = inverse.dplot;  catch, dp   = 0;        end


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
Nmod  = numel(modalities);                  % number of modalities
Nmax  = Nr;                                 % max number of temporal modes
 
 
% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================
for i = 1:Nl
    
    fprintf('Checking lead fields for subject %i\n',i)
    [L,D{i}] = spm_eeg_lgainmat(D{i});
    
    for m = 1:Nmod
        
        % Check gain or lead-field matrices
        %------------------------------------------------------------------
        Ic{i,m}  = indchantype(D{i}, modalities{m}, 'GOOD');
        Nd(i)    = size(L,2);
        Nc(i,m)  = length(Ic{i,m});
        
        if isempty(Ic{i,m})
            errordlg(['Modality ' modalities{m} ' is missing from file ' D{i}.fname]);
            return
        end
        
        if any(diff(Nd))
            errordlg('Please ensure subjects have the same number of dipoles.')
            return
        end
        
        % Check for null space over sensors (SX) and remove it
        %------------------------------------------------------------------
        try
            SX     = D{i}.sconfounds{m};
            R{i,m} = speye(Nc(i,m),Nc(i,m)) - SX*spm_pinv(SX);
        catch
            R{i,m} = speye(Nc(i,m),Nc(i,m));
        end
    end
end
fprintf(' - done\n')
 
 
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================
 
fprintf('%-40s: %30s','Green function from graph Laplacian','...computing'); %-#

Nd    = Nd(1);                                          % number of dipoles
vert  = D{1}.inv{D{1}.val}.mesh.tess_mni.vert;
face  = D{1}.inv{D{1}.val}.mesh.tess_mni.face;
A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
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

fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
 
 
% Check for (e.g., empty-room) sensor components (in Qe{1})
%==========================================================================
QE    = cell(Nl,Nmod);
for i = 1:Nl
    for m = 1:Nmod
        try
            QE{i,m} = D{i}.inv{D{i}.val}.inverse.Qe{m};
            QE{i,m} = Nc(i,m)*QE{i,m}/trace(QE{i,m});
            if length(QE{i,m}) ~= Nc(i,m)
                errordlg('error component (modality %s; subject %d) does not match number of channels (%d)\n',modalities{m},i,Nc(i,m))
                return
            end
            fprintf('Using sensor error component provided...\n');
            
        % assume i.i.d. if not specified
        %------------------------------------------------------------------
        catch
            QE{i,m} = 1; 
        end
    end
end
 
 
%==========================================================================
% Spatial projectors (adjusting for different Lead-fields)
%==========================================================================
 
fprintf('Optimising and aligning spatial modes ...\n')
 
% Use recursive (regularised) least squares to find average lead-field
%--------------------------------------------------------------------------
Is    = 1:Nd;
UL    = cell(Nmod,1);
for m = 1:Nmod
    
    % Initialise average lead-field with L{1}
    %----------------------------------------------------------------------
    UL{m} = R{1,m}*spm_eeg_lgainmat(D{1},Is,D{1}.chanlabels(Ic{1,m}));
    AA    = 1;
    
    % pre-compute regularised inverses (for speed)
    %----------------------------------------------------------------------
    for i = 1:Nl
        L     = R{i,m}*spm_eeg_lgainmat(D{i},Is,D{i}.chanlabels(Ic{i,m}));
        iL{i} = spm_inv(L*L');
    end
    
    % Optimise alignment matrices A such that A{m}*L{m} = <A{m}*L{m}>m
    %----------------------------------------------------------------------
    for j = 1:8
        
        % eliminate redundant virtual channels
        %------------------------------------------------------------------
        fprintf('Aligning - iteration: %i\n',j)
        UL{m} = spm_sqrtm(spm_inv(AA))*UL{m};
        
        % eliminate low SNR spatial modes
        %------------------------------------------------------------------
        U     = spm_svd((UL{m}*UL{m}'));
        UL{m} = U'*UL{m};
        Nm(m) = size(UL{m},1);
        
        % normalise lead-field
        %------------------------------------------------------------------
        Scale = sqrt(trace(UL{m}*UL{m}')/Nm(m));
        UL{m} = UL{m}/Scale;
        
        % spatial projectors A{i,m) for i = 1,...,Nl subjects
        %------------------------------------------------------------------
        AL    = 0;
        AA    = 0;
        for i = 1:Nl
            L      = R{i,m}*spm_eeg_lgainmat(D{i},Is,D{i}.chanlabels(Ic{i,m}));
            A{i,m} = UL{m}*L'*iL{i};
            AL     = AL + A{i,m}*L;
            AA     = AA + A{i,m}*A{i,m}';
        end
        
        % re-compute average
        %------------------------------------------------------------------
        UL{m} = AL/Nl;
        if Nl == 1, break, end
        
    end
    
    % Report
    %----------------------------------------------------------------------
    fprintf('Using %d spatial modes for modality %s\n',Nm(m),modalities{m})
    
    if dp
        spm_figure('GetWin','Lead fields');
        for i = 1:Nl
            L  = R{i,m}*spm_eeg_lgainmat(D{i},Is(1:8:end),D{i}.chanlabels(Ic{i,m}));
            subplot(Nl,4,(i - 1)*4 + 1)
            imagesc(A{i,m}*A{i,m}'), axis square
            subplot(Nl,4,(i - 1)*4 + 2)
            imagesc(A{i,m}'*A{i,m}), axis square
            subplot(Nl,4,(i - 1)*4 + 3)
            plot(sum(A{i,m}.^2,2)),  axis square
            subplot(Nl,4,(i - 1)*4 + 4)
            plot(L'*A{i,m}'),        axis square tight
        end
        drawnow
    end
end
 
 
% check restriction: assume radii are the same for all (Nv) VOI
%==========================================================================
Nv  = size(xyz,1);
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end
 
% Restrict source space to Ns sources by eliminating dipoles
%--------------------------------------------------------------------------
if any(any(xyz)) || ~isempty(mask)
    
    Is    = sparse(Nd,1);
    
    if any(any(xyz))
        for i = 1:Nv
            Iv = sum([vert(:,1) - xyz(i,1), ...
                vert(:,2) - xyz(i,2), ...
                vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
            Is = Is | Iv;
        end
    end
    
    if ~isempty(mask)
        Iv = spm_mesh_project(struct('vertices',vert,'faces',face), mask);
        Is = Is | Iv(:);
    end
    
    Is    = find(Is);
else
    Is    = 1:Nd;
end


vert  = vert(Is,:);
QG    = QG(Is,Is);
for m = 1:Nmod
    UL{m} = UL{m}(:,Is);
end
Ns    = length(Is);
 
 
 
 
%==========================================================================
% Temporal projector
%==========================================================================
 
% loop over Nl lead-fields (subjects)
%--------------------------------------------------------------------------
Nn    = zeros(1,Nl);                             % number of samples
AY    = {};                                      % pooled response for MVB
AYYA  = 0;                                       % pooled response for ReML
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
    pst{i} = 1000*D{i}.time;                     % peristimulus time (ms)
    pst{i} = pst{i}(It{i});                      % windowed time (ms)
    dur    = (pst{i}(end) - pst{i}(1))/1000;     % duration (s)
    dct{i} = (It{i} - It{i}(1))/2/dur;           % DCT frequencies (Hz)
    Nb(i)  = length(It{i});                      % number of time bins
    
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
    
    % notch filter nf (Hz)
    %----------------------------------------------------------------------
    % nf   = 10.2/1000;
    try
        T0 = [sin(2*pi*nf*pst{i}(:)) cos(2*pi*nf*pst{i}(:))];
        T  = T - T0*pinv(T0)*T;
    end
    
    % Hanning operator (if requested)
    %----------------------------------------------------------------------
    if Han
        W  = sparse(1:Nb(i),1:Nb(i),spm_hanning(Nb(i)));
    else
        W  = 1;
    end
    
    % get trials or conditions
    %----------------------------------------------------------------------
    try
        trial = D{i}.inv{D{i}.val}.inverse.trials;
    catch
        trial = D{i}.condlist;
    end
    Nt(i) = length(trial);
    
    
    % get temporal covariance (Y'*Y) to find temporal modes
    %======================================================================
    MY    = cell(Nmod,1);                        % mean response
    YTY   = sparse(0);                           % accumulator
    for m = 1:Nmod                               % loop over modalities
        
        % get (spatially aligned) data
        %------------------------------------------------------------------
        N     = 0;
        YY    = 0;
        MY{m} = 0;
        for j = 1:Nt(i)                          % pool over conditions
            c     = D{i}.indtrial(trial{j});     % and trials
            Nk    = length(c);
            for k = 1:Nk
                Y     = A{i,m}*D{i}(Ic{i,m},It{i},c(k));
                MY{m} = MY{m} + Y;
                YY    = YY + Y'*Y;
                N     = N + Nb(i);
            end
        end
        
        % Apply any Hanning and filtering
        %------------------------------------------------------------------
        YY         = W'*YY*W;
        YY         = T'*YY*T;
        
        % Scale data (to remove subject and modality scaling differences)
        %------------------------------------------------------------------
        scale(i,m) = sign(trace(MY{m}'*(UL{m}*UL{1}')*MY{1}));
        scale(i,m) = scale(i,m)/sqrt(trace(YY)/(Nm(m)*N));
        YTY        = YTY + YY*(scale(i,m)^2);
        
    end
    
    % temporal projector (at most Nmax modes) S = T*V
    %======================================================================
    if STAT % Stationarity assumption
        
        S{i}  = T;                               % temporal modes
        Nr(i) = size(T,2);                       % number of temporal modes
        VE(i) = 1;                               % variance explained
        
    else
        [U,E] = spm_svd(YTY,exp(-8));            % get temporal modes
        E     = diag(E)/trace(YTY);              % normalise variance
        Nr(i) = min(length(E),Nmax);             % number of temporal modes
        S{i}  = T*U(:,1:Nr(i));                  % temporal modes
        VE(i) = sum(E(1:Nr(i)));                 % variance explained
    end
 
    fprintf('Using %i temporal modes for subject %i, ',Nr(i),i)
    fprintf('accounting for %0.2f percent average variance\n',full(100*VE(i)))
    
    % whitening
    %----------------------------------------------------------------------
    Vq{i}  = S{i}*inv(S{i}'*qV{i}*S{i})*S{i}';   % temporal precision
    
    
    % get spatial covariance (Y*Y') for Gaussian process model.
    %======================================================================
    
    % loop over Nt trial types
    %----------------------------------------------------------------------
    UYYU{i} = 0;
    for j = 1:Nt(i)
        
        UY{i,j} = sparse(0);
        c       = D{i}.indtrial(trial{j});
        Nk      = length(c);
        
        % loop over epochs
        %------------------------------------------------------------------
        for k = 1:Nk
            
            % stack (scaled aligned data) over modalities
            %--------------------------------------------------------------
            for m = 1:Nmod
                Y       = D{i}(Ic{i,m},It{i},c(k))*S{i};
                MY{m}   = A{i,m}*Y*scale(i,m)/Nk;
            end
            
            % accumulate first & second-order responses
            %--------------------------------------------------------------
            Nn(i)       = Nn(i) + Nr(i);         % number of samples
            Y           = spm_cat(MY);           % contribution to ERP
            YY          = Y*Y';                  % and covariance
            
            % accumulate statistics (subject-specific)
            %--------------------------------------------------------------
            UY{i,j}     = UY{i,j} + Y;           % condition-specific ERP
            UYYU{i}     = UYYU{i} + YY;          % subject-specific covariance
            
            % and pool for optimisation of spatial priors over subjects
            %--------------------------------------------------------------
            AY{end + 1} = Y;                     % pooled response for MVB
            AYYA        = AYYA    + YY;          % pooled response for ReML
            
        end
    end
end
 
% and concatenate for optimisation of spatial priors over subjects
%--------------------------------------------------------------------------
AY    = spm_cat(AY);                             % pooled response for MVB
UL    = spm_cat(UL);                             % pooled lead fields
 
 
% generate sensor error components (Qe)
%==========================================================================
AQ{1} = sparse(0);
for m = 1:Nmod
    Qe{m} = sparse(0);
end
 
% assuming equal noise over subjects (Qe{m}) and modalities AQ
%--------------------------------------------------------------------------
N     = cell(Nmod,Nmod);
for i = 1:Nl
    for m = 1:Nmod
        N{m,m} = sparse(Nm(m),Nm(m));
    end
    for m = 1:Nmod
        Q      = N;
        AQeA   = A{i,m}*QE{i,m}*A{i,m}';
        Q{m,m} = AQeA/(trace(AQeA)/Nm(m));
        Q      = spm_cat(Q)/Nl;
        Qe{m}  = Qe{m} + Q;
        AQ{1}  = AQ{1} + Q;
    end
end
 
 
%==========================================================================
% Step 1: Optimise spatial priors over subjects
%==========================================================================
 
% create source components (Qp)
%==========================================================================
switch(type)
    
    case {'MSP','GS','ARD','BMR'}
        
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        LQL   = {};
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
        
    case {'EBB'}
        % create beamforming prior. See:
        % Source reconstruction accuracy of MEG and EEG Bayesian inversion approaches. 
        % Belardinelli P, Ortiz E, Barnes G, Noppeney U, Preissl H. PLoS One. 2012;7(12):e51985. 
        %------------------------------------------------------------------
        InvCov      = spm_inv(YY);
        allsource   = zeros(Ns,1);
        Sourcepower = zeros(Ns,1);
        for bk = 1:Ns
            normpower       = 1/(UL(:,bk)'*UL(:,bk));
            Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
            allsource(bk)   = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);   % Normalise
        
        Qp{1}     = diag(allsource);
        LQpL{1}   = UL*diag(allsource)*UL';
        
end
 
 
% augment with exogenous (e.g., fMRI) source priors in pQ
%==========================================================================
for i = 1:length(pQ)
    
    switch(type)
        
        case {'MSP','GS','ARD','BMR'}
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
        MVB   = spm_mvb(AY,UL,[],Q,AQ,16);
        
        % Accumulate empirical priors
        %------------------------------------------------------------------
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
        
end
 
switch(type)
    
    case {'BMR'}
        
        % convert patterns into covariance components
        %------------------------------------------------------------------
        Np    = length(Qp);
        for i = 1:Np
            LQpL{i} = LQpL{i}.q*LQpL{i}.q';
        end
        
        % hyperparameter estimation
        %------------------------------------------------------------------
        [C,h,Ph,F,Fa,Fc,Eh,Ch,hE,hC] = spm_reml_sc(AYYA,[],[Qe LQpL],sum(Nn),-16,32);
        
        
        % Bayesian model reduction
        %------------------------------------------------------------------
        DCM.M.pE = hE;
        DCM.M.pC = hC;
        DCM.Ep   = Eh;
        DCM.Cp   = Ch;
        h        = spm_dcm_sparse(DCM);
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h((1:Np) + Ne);
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i}.q*Qp{i}.q';
        end
        
        % Accumulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
        
end


switch(type)
    
    case {'MSP','ARD'}
        
        % ReML - ARD
        %------------------------------------------------------------------
        [C,h] = spm_sp_reml(AYYA,[],[Qe LQpL],sum(Nn));
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h((1:Np) + Ne);
        qp    = sparse(0);
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
    
    case {'IID','MMN','LOR','COH','EBB'}
        
        % or ReML - ARD
        %------------------------------------------------------------------
        Q0    = exp(-2)*trace(AYYA)/sum(Nn)*AQ{1}/trace(AQ{1});
        [C,h] = spm_reml_sc(AYYA,[],[Qe LQpL],sum(Nn),-4,16,Q0);
        
        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h((1:Np) + Ne);
        qp    = sparse(0);
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
    AQ    = 0;
    Qe    = {};
    for m = 1:Nmod
        N{m,m} = sparse(Nm(m),Nm(m));
    end
    for m = 1:Nmod
        Q      = N;
        AQeA   = A{i,m}*QE{i,m}*A{i,m}';
        Q{m,m} = AQeA/(trace(AQeA)/Nm(m));
        Qe{m}  = spm_cat(Q);
        AQ     = AQ + Qe{m};
    end
    
    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    Np    = length(LQPL);
    Ne    = length(Qe);
    Q     = [Qe LQPL];
    
    % re-do ReML (with informative hyperpriors)
    %======================================================================
    Q0          = exp(-2)*trace(UYYU{i})/Nn(i)*AQ/trace(AQ);
    [Cy,h,Ph,F] = spm_reml_sc(UYYU{i},[],Q,Nn(i),-4,16,Q0);
    
    % Data ID
    %----------------------------------------------------------------------
    ID    = spm_data_id(AYYA);
    
    
    % Covariance: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    Cp    = sparse(0);
    LCp   = sparse(0);
    hp    = h(Ne + (1:Np));
    for j = 1:Np
        Cp  =  Cp + hp(j)*QP{j};
        LCp = LCp + hp(j)*LQP{j};
    end
    
    % MAP estimates of instantaneous sources
    %======================================================================
    M     = LCp'/Cy;
    
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
        SSR  = SSR + sum(var(full(UY{i,j} - UL*J{j}),0,2));
        SST  = SST + sum(var(full(UY{i,j}),0,2));
        
    end
    
    % accuracy; signal to noise (over sources)
    %======================================================================
    R2   = 100*(SST - SSR)/SST;
    fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE(i)));
    
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
    inverse.R2     = R2;                   % variance explained (reduced)
    inverse.VE     = R2*VE(i);             % variance explained
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
    
    % Display
    %======================================================================
    if ~spm('CmdLine'), spm_eeg_invert_display(D{i}); end
    
end
 
if length(D) == 1, D = D{1}; end

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
