function [D] = spm_eeg_invert(D, val)
% ReML inversion of multiple forward models for EEG-MEG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG data
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.trials - D.events.types to invert
%     inverse.modality modality to use in case of multimodal datasets
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
%     inverse.Is     - Indices of active dipoles
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.R2     - variance accounted for (%)
%__________________________________________________________________________
% 
% This routine solves ill-posed linear models of the following form
%
% [AY{1}...AY{n}] = L(1} * [J{1}...J{1}]   +  [e{1}...e{n}]
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_eeg_invert.m 3173 2009-06-02 14:50:19Z karl $
 
% check whether this is a group inversion
%--------------------------------------------------------------------------
if ~iscell(D), D = {D}; end
Nl         = length(D);                          % number of forward models
 
% D - SPM data structure
%==========================================================================
if nargin == 1
    val = 1;
end
 
inverse    = D{1}.inv{D{1}.val}.inverse;
 
% defaults
%--------------------------------------------------------------------------
try, type  = inverse.type;   catch, type  = 'GS';              end
try, s     = inverse.smooth; catch, s     = 0.6;               end
try, Np    = inverse.Np;     catch, Np    = 256;               end
try, Nm    = inverse.Nm;     catch, Nm    = 128;               end
try, Nr    = inverse.Nr;     catch, Nr    = 8;                 end
try, xyz   = inverse.xyz;    catch, xyz   = [0 0 0];           end
try, rad   = inverse.rad;    catch, rad   = 128;               end
try, lpf   = inverse.lpf;    catch, lpf   = 1;                 end
try, hpf   = inverse.hpf;    catch, hpf   = 256;               end
try, sdv   = inverse.sdv;    catch, sdv   = 4;                 end
try, Han   = inverse.Han;    catch, Han   = 1;                 end
try, Na    = inverse.Na;     catch, Na    = 1024;              end
try, woi   = inverse.woi;    catch, woi   = [];                end

try
    modality = inverse.modality;  
catch
    modality = D{1}.inv{D{1}.val}.forward(1).modality;
end
 
%==========================================================================
% Spatial parameters
%==========================================================================
 
% Check gain or lead-field matrices
%--------------------------------------------------------------------------
for i = 1:Nl
    [L D{i}] = spm_eeg_lgainmat(D{i});
    
    Nd(i)    = size(L,2);                             % number of dipoles

    Ic{i} = setdiff(meegchannels(D{i}, modality), badchannels(D{i}));
    
    if isempty(Ic{i})
        error(['The specified modality is missing from file ' D{i}.fname]);
    end

    Nc(i)    = length(Ic{i});                         % number of channels
    
    if any(diff(Nd))
        warndlg('please ensure the number of dipoles is the same')
        return
    end
end

% chaeck restriction; assume radii are the same for all VOI
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
clear Qi
QG    = QG.*(QG > exp(-8));
QG    = QG*QG;
fprintf(' - done\n')
 
 
% Restrict source space
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

% Check for null space over sensors and remove it
%--------------------------------------------------------------------------
for i = 1:Nl
    try
        SX   = D{i}.sconfounds;
        R{i} = speye(Nc(i),Nc(i)) - SX*pinv(SX);
    catch
        R{i} = speye(Nc(i),Nc(i));
    end
end
 
% Project to channel modes (U); Checking for null space over sensors
%--------------------------------------------------------------------------
TOL   = 16;
Nmax  = Nm;
for i = 1:Nl
    L     = R{i}*spm_eeg_lgainmat(D{i},Is,D{i}.chanlabels(Ic{i}));
    U{i}  = spm_svd(L*L',exp(-TOL));
    Nm(i) = min(size(U{i},2),Nmax);
    U{i}  = U{i}(:,1:Nm(i));
    UL{i} = U{i}'*L;
end
fprintf('Using %i spatial modes\n',Nm)
 
 
 
%==========================================================================
% Temporal parameters
%==========================================================================

Nrmax = Nr;
A     = {};
AY    = {};
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
 
    % Peri-stimulus time
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
    j      = find( (dct{i} > lpf) & (dct{i} < hpf) );
    T      = T(:,j);
    dct{i} = dct{i}(j);
 
    % get data (with temporal filtering)
    %======================================================================
    
    % get trials and channels
    %----------------------------------------------------------------------
    try
        trial = D{i}.inv{D{i}.val}.inverse.trials; 
    catch
        trial = unique(D{i}.conditions);
    end
    Nt(i) = length(trial);
      
    for j = 1:Nt(i)
        Y{i,j} = sparse(0);
        c     = D{i}.pickconditions(trial{j});
        Ne    = length(c);
        for k = 1:Ne
            Y{i,j} = Y{i,j} + R{i}*squeeze(D{i}(Ic{i},It{i},c(k)))*T/Ne;
        end
    end
 
    % Hanning operator (if requested)
    %----------------------------------------------------------------------
    if Han
        W  = T'*sparse(1:Nb(i),1:Nb(i),spm_hanning(Nb(i)))*T;
        WY = W'*spm_cat(Y(i,:)')';
    else
        WY = spm_cat(Y(i,:)')';
    end
 
    % temporal projector (at most 8 modes) S = T*v
    %======================================================================
    [v u]  = spm_svd(WY,1);                      % temporal modes
    Nr(i)  = min(size(v,2),Nrmax);               % number of temporal modes
    v      = v(:,      1:Nr(i));
    u      = u(1:Nr(i),1:Nr(i));
    VE(i)  = sum(sum(u.^2))/sum(sum(WY.^2));     % variance explained
    S{i}   = T*v;                                % temporal projector
    iV{i}  = inv(S{i}'*qV{i}*S{i});              % precision (mode)
    Vq{i}  = S{i}*iV{i}*S{i}';                   % precision (time)
 
    % spatial projector (adjusting for different Lead-fields)
    %======================================================================
    A{i}   = UL{1}*pinv(full(UL{i}))*U{i}';
 
    % spatially adjust, temporally whiten and scale under i.i.d priors
    %----------------------------------------------------------------------
    for j = 1:Nt(i)
        Y{i,j}      = Y{i,j}*v;
        AY{end + 1} = A{i}*Y{i,j}*sqrtm(iV{i});
        AY{end}     = AY{end}*sqrt(trace(A{i}*A{i}'))/trace(AY{end}'*AY{end});
    end
 
    % create sensor components (Qe)
    %----------------------------------------------------------------------
    Qe{i} = A{i}*A{i}';
 
end
 
% adjusted data and sample covariance
%--------------------------------------------------------------------------
[Y scale] = spm_cond_units(Y);
AY        = spm_cat(AY);
YY        = AY*AY';
G         = UL{1};
 
fprintf('Using %i temporal modes\n',Nr)
fprintf('accounting for %0.2f percent variance\n',full(100*VE))
    
 
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
            LQpL{end + 1}.q = G*q;
 
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                             vert(:,2) - vert(Ip(i),2), ...
                             vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = G*q;
 
            % bilateral
            %--------------------------------------------------------------
            q               = QG(:,Ip(i)) + QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = G*q;
 
        end
 
    case {'LOR','COH'}
 
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';
 
        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = G*Qp{2}*G';
 
 
    case {'IID','MMN'}
 
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';
 
end
 
% Inverse solution
%==========================================================================
QP    = {};
 
% get source-level priors (using all subjects)
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
        MVB   = spm_mvb(AY,G,[],Q,Qe,16);
 
        % Spatial priors (QP); eliminating minor patterns
        %------------------------------------------------------------------
        cp    = diag(MVB.cp);
        [i j] = sort(-diag(MVB.cp));
        try
            j = j(1:256);
        end
        qp    = Q(:,j)*MVB.cp(j,j)*Q(:,j)';
 
        % Accmulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;
 
end
 
switch(type)

    case {'MSP','ARD'}

        % or ReML - ARD
        %------------------------------------------------------------------
        qp          = sparse(0);
        Q           = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,sum(Nr)*sum(Nt));

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

        % Accmulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;


    case {'IID','MMN','LOR','COH'}

        % or ReML - ARD
        %------------------------------------------------------------------
        qp          = sparse(0);
        Q           = {Qe{:} LQpL{:}};
        [Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,sum(Nr)*sum(Nt));

        % Spatial priors (QP)
        %------------------------------------------------------------------
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h([1:Np] + Ne);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end

        % Accmulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;

end


 
% re-estimate (one subject at a time)
%==========================================================================
for i = 1:Nl
 
    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    L     = R{i}*spm_eeg_lgainmat(D{i}, Is, D{i}.chanlabels(Ic{i}));
    Qe    = R(i);
    Ne    = length(Qe);
    Np    = length(QP);
    LQpL  = {};
    for j = 1:Np
        LQpL{j}  = L*QP{j}*L';
    end
    Q     = {Qe{:} LQpL{:}};
    YY    = spm_cat(Y(i,:))*kron(speye(Nt(i),Nt(i)),iV{i})*spm_cat(Y(i,:))';
 
    % re-do ReML
    %----------------------------------------------------------------------
    [Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,Nr(i)*Nt(i));
 
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
    for j = 1:Nt(i)
        
        % trial-type specific source reconstruction
        %------------------------------------------------------------------
        J{j} = M*Y{i,j};
 
        % sum of squares
        %------------------------------------------------------------------
        SSR  = SSR + sum(var((Y{i,j} - L*J{j}),0,2));
        SST  = SST + sum(var(Y{i,j},0,2));
 
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
    inverse.scale  = scale;                % scale-factor
 
    inverse.M      = M;                    % MAP projector (reduced)
    inverse.J      = J;                    % Conditional expectation
    inverse.Y      = Y(i,:);               % ERP data (reduced)
    inverse.L      = L;                    % Lead-field (reduced)
    inverse.R      = speye(Nc(i),Nc(i));   % Re-referencing matrix
    inverse.qC     = Cq;                   % spatial covariance
    inverse.qV     = Vq{i};                % temporal correlations
    inverse.T      = S{i};                 % temporal subspace
    inverse.U      = speye(Nc(i),Nc(i));   % spatial subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It{i};                % Indices of time bins
    inverse.Ic     = Ic{i};                % Indices of good channels
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.pst    = pst{i};               % peristimulus time
    inverse.dct    = dct{i};               % frequency range
    inverse.F      = F;                    % log-evidence
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.VE     = VE(i);                % variance explained
    inverse.woi    = w{i};                 % time-window inverted
    
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
