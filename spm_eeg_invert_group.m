function [D] = spm_eeg_invert_group(D)
% ReML inversion of multiple forward models for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the
% spatiotemporal hierarchy implicit in EEG data
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.trials - D.events.types to invert
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.Nm     - maximum number of channel modes
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'MSP' multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' LORETA and minimum norm
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%     inverse.lpf    - band-pass filter - low  frequency cutoff (Hz)
%     inverse.hpf    - band-pass filter - high frequency cutoff (Hz)
%     inverse.Lap    - switch for Laplace transform
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Han    - switch for Hanning window
%     inverse.Na     - number of most energetic dipoles
%     inverse.woi    - time window of interest ([start stop] in ms)
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J      - Conditional expectation
%     inverse.L      - Lead field (reduced)
%     inverse.R      - Re-referencing matrix
%     inverse.qC     - spatial  covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal subspace
%     inverse.U      - spatial  subspace
%     inverse.Is     - Indices of active dipoles
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.R2     - variance accounted for (%)
%__________________________________________________________________________
% Karl Friston: Modified by Rik Henson to allow more arguments to be passed       4/6/07
 
% check whether this is a group inversion
%--------------------------------------------------------------------------
if ~iscell(D), D = {D}; end
Nl         = length(D);                          % number of forward models
 
% D - SPM data structure
%==========================================================================
inverse    = D{1}.inv{D{1}.val}.inverse;
 
% defaults
%--------------------------------------------------------------------------
try, trial = inverse.trials; catch, trial = D{1}.events.types; end
try, type  = inverse.type;   catch, type  = 'GS';              end
try, s     = inverse.smooth; catch, s     = 0.6;               end
try, Np    = inverse.Np;     catch, Np    = 256;               end
try, Nm    = inverse.Nm;     catch, Nm    = 128;                end
try, xyz   = inverse.xyz;    catch, xyz   = [0 0 0];           end
try, rad   = inverse.rad;    catch, rad   = 128;               end
try, lpf   = inverse.lpf;    catch, lpf   = 1;                 end
try, hpf   = inverse.hpf;    catch, hpf   = 256;               end
try, sdv   = inverse.sdv;    catch, sdv   = 4;                 end
try, Han   = inverse.Han;    catch, Han   = 1;                 end
try, Na    = inverse.Na;     catch, Na    = 1024;              end
try, woi   = inverse.woi;    catch, woi   = [];                end
 
% Load Gain or Lead field matrix
%--------------------------------------------------------------------------
for i = 1:Nl
    gainmat   = D{i}.inv{D{i}.val}.forward.gainmat;
    try
        G     = load(gainmat);
    catch
        [p f] = fileparts(gainmat);
        G     = load(f);
    end
    name   = fieldnames(G);
    L{i}   = sparse(getfield(G, name{1}));
    Nc(i)  = size(L{i},1);                         % number of channels
end
 
% Time-window of interest
%--------------------------------------------------------------------------
if isempty(woi)
    woi = round([-D{1}.events.start D{1}.events.stop]*1000/D{1}.Radc);
end
It    = round(woi*(D{1}.Radc/1000)) + D{1}.events.start;
It    = max(1,It(1)):min(It(end),size(D{1}.data,2));
 
% Peri-stimulus time
%--------------------------------------------------------------------------
pst   = (It - D{1}.events.start - 1);
pst   = pst/D{1}.Radc*1000;                      % peristimulus time (ms)
dur   = (pst(end) - pst(1))/1000;                % duration (s)
dct   = (It - It(1))/2/dur;                      % DCT frequenices (Hz)
 
% parameters
%==========================================================================
Nd    = size(L{1},2);                            % number of dipoles
Nb    = length(It);                              % number of time bins
Nv    = size(xyz,1);                             % number of VOI
 
% assume radii are the same for all VBOI
%--------------------------------------------------------------------------
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end
 
 
% Serial correlations
%--------------------------------------------------------------------------
K     = exp(-(pst - pst(1)).^2/(2*sdv^2));
K     = toeplitz(K);
qV    = sparse(K*K');
 
% Confounds and temporal subspace
%--------------------------------------------------------------------------
T     = spm_dctmtx(Nb,Nb);
i     = (dct > lpf) & (dct < hpf);
T     = T(:,i);
dct   = dct(i);
 
% get data (with temporal filtering)
%==========================================================================
Nt    = length(trial);
for k = 1:Nl
    Ic{k} = setdiff(D{k}.channels.eeg, D{k}.channels.Bad);
    for i = 1:Nt
        Y{k,i}  = sparse(0);
        if isfield(D{k}.events,'reject')
            c = find(D{k}.events.code == trial(i) & ~D{k}.events.reject);
        else
            c = find(D{k}.events.code == trial(i));
        end
        for j = 1:length(c)
            Y{k,i} = Y{k,i} + squeeze(D{k}.data(Ic{k},It,c(j)))*T;
        end
    end
end
 
% temporal covariance (with Hanning if requested)
%--------------------------------------------------------------------------
if Han
    W = T'*diag(spm_hanning(Nb))*T;
else
    W = T'*T;
end
YY    = W'*spm_cat(Y(:))'*spm_cat(Y(:))*W;
 
% project onto temporal modes
%--------------------------------------------------------------------------
S     = spm_svd(YY,exp(-4));
T     = T*S;                                     % temporal projector
iV    = inv(T'*qV*T);                            % precision (mode)
sV    = sqrtm(iV);                               % sqrt (mode)
Vq    = T*iV*T';                                 % precision (time)
Nr    = size(T,2);                               % number of temporal modes
for i = 1:Nl
    for j = 1:Nt
        Y{i,j} = Y{i,j}*S;
    end
end
fprintf('Using %i temporal modes\n',Nr)
 
% sample covariance over time (adjusting for different Lead-fields)
%--------------------------------------------------------------------------
G     = L{1};
AY    = {};
YY    = sparse(Nc(1),Nc(1));
for i = 1:Nl
    R = speye(Nc(i),Nc(i))*exp(-32)*norm(L{i}*L{i}',1);
    A = (G*L{i}')*inv((L{i}*L{i}' + R));
    for j = 1:Nt
        AY{end + 1} = A*Y{i,j}*sV;
        YY          = YY + AY{end}*AY{end}';
    end
end
AY    = spm_cat(AY);
 
 
% Project to channel modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(G*G',exp(-32));
try
    U = U(:,1:Nm);
end
Nm    = size(U,2);
AY    = U'*AY;
YY    = U'*YY*U;
G     = U'*G;
fprintf('Using %i spatial modes\n',Nm)
 
 
% Restrict source space
%==========================================================================
vert  = D{1}.inv{D{1}.val}.mesh.tess_mni.vert;
face  = D{1}.inv{D{1}.val}.mesh.tess_mni.face;
Is    = sparse(Nd,1);
for i = 1:Nv
    Iv = sum([vert(:,1) - xyz(i,1), ...
        vert(:,2) - xyz(i,2), ...
        vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
    Is = Is | Iv;
end
Is    = find(Is);
vert  = vert(Is,:);
for i = 1:Nl
    L{i} = L{i}(:,Is);
end
G     = G(:,Is);
Ns    = length(Is);
 
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================
if ~strcmp(type,'IID')
 
    fprintf('Computing Green''s function from graph Laplacian:')
    %----------------------------------------------------------------------
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
    QG    = QG(Is,Is);
    fprintf(' - done\n')
    
end
 
 
% covariance components
%==========================================================================
 
% sensor noise (accommodating spatial projector)
%--------------------------------------------------------------------------
Qe    = {U'*U};
 
% create source components
%--------------------------------------------------------------------------
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
 
% get source-level hyperparameters (using all subjects)
%--------------------------------------------------------------------------
if strcmp(type,'GS')
 
    % Greedy search over MSPs
    %----------------------------------------------------------------------
    Np    = length(Qp);
    QP    = sparse(0);
    Q     = sparse(Ns,Np);
    for i = 1:Np
        Q(:,i) = Qp{i}.q;
    end
 
    % Multivariate Bayes
    %----------------------------------------------------------------------
    MVB   = spm_mvb(AY,G,[],Q,Qe,16,1/2);
 
    % Spatial priors (QP); eliminating minor patterns
    %----------------------------------------------------------------------
    QP    = MVB.Cp;
    pV    = diag(QP);
    i     = find(pV > max(pV)/128);
    QP    = Q(:,i)*QP(i,i)*Q(:,i)';
 
else
 
    % or ReML - ARD
    %----------------------------------------------------------------------
    QP          = sparse(0);
    Q           = {Qe{:} LQpL{:}};
    [Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,Nr*Nt);
 
    % Spatial priors (QP)
    %----------------------------------------------------------------------
    Ne    = length(Qe);
    Np    = length(Qp);
    hp    = h([1:Np] + Ne);
    for i = 1:Np
        if hp(i) > max(hp)/128;
            try
                QP  = QP + hp(i)*Qp{i}.q*Qp{i}.q';
            catch
                QP  = QP + hp(i)*Qp{i};
            end
        end
    end
end
 
% eliminate sources with low [empirical] prior variance
%--------------------------------------------------------------------------
pV    = diag(QP);
i     = find(pV > max(pV)/exp(16));
Is    = Is(i);
QP    = QP(i,i);
for j = 1:Nl
    L{j} = L{j}(:,i);
end
Ns    = length(Is);
 
% re-estimate (one subject at a time)
%==========================================================================
for i = 1:Nl
 
    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    Qe    = speye(Nc(i),Nc(i));
    LQpL  = L{i}*QP*L{i}';
    Q     = {Qe LQpL};
    YY    = spm_cat(Y(i,:))*kron(speye(Nt,Nt),iV)*spm_cat(Y(i,:))';
 
    % re-do ReML
    %----------------------------------------------------------------------
    [Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,Nr*Nt);
 
    % Covariances: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    hp    = h(2);
    LCp   = hp*L{i}*QP;
 
    % MAP estimates of instantaneous sources
    %======================================================================
    iC    = inv(Cy);
    M     = LCp'*iC;
 
    % conditional covariance (leading diagonal)
    % Cq    = Cp - Cp*L'*iC*L*Cp;
    %----------------------------------------------------------------------
    Cq    = hp*diag(QP) - sum(LCp.*M')';
 
    % evaluate conditional expectation (of the sum over trials)
    %----------------------------------------------------------------------
    SSR   = 0;
    SST   = 0;
    for j = 1:Nt
        
        % trial-type specific source reconstruction
        %------------------------------------------------------------------
        J{j}    = M*Y{i,j};
 
        % sum of squares
        %------------------------------------------------------------------
        SSR   = SSR + sum(var((Y{i,j} - L{i}*J{j}),0,2));
        SST   = SST + sum(var(Y{i,j},0,2));
 
    end
 
    % assess accuracy; signal to noise (over sources)
    %======================================================================
    R2   = 100*(SST - SSR)/SST;
    fprintf('Variance explained %.2f (percent)\n',R2)
 
    % Save results
    %======================================================================
    U              = speye(Nc(i),Nc(i));
    inverse.type   = type;                 % inverse model
    inverse.smooth = s;                    % smoothness (0 - 1)
    inverse.xyz    = xyz;                  % VOI (XYZ)
    inverse.rad    = rad;                  % VOI (radius)
 
    inverse.M      = M;                    % MAP projector (reduced)
    inverse.J      = J;                    % Conditional expectation
    inverse.Y      = Y(i,:);               % ERP data (reduced)
    inverse.L      = L{i};                 % Lead-field (reduced)
    inverse.R      = U;                    % Re-referencing matrix
    inverse.qC     = Cq;                   % spatial covariance
    inverse.qV     = Vq;                   % temporal correlations
    inverse.T      = T;                    % temporal subspace
    inverse.U      = U;                    % spatial subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It;                   % Indices of time bins
    inverse.Ic     = Ic{i};                % Indices of good channels
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.pst    = pst;                  % peristimulus time
    inverse.dct    = dct;                  % frequency range
    inverse.F      = F;                    % log-evidence
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.woi    = woi;		   % timewindow inverted
 
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
