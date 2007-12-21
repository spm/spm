function [D] = spm_eeg_invert(D)
% ReML fusion of multiple forward models for EEG-EMG
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
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
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
try, Nm    = inverse.Nm;     catch, Nm    = 128;               end
try, xyz   = inverse.xyz;    catch, xyz   = [0 0 0];           end
try, rad   = inverse.rad;    catch, rad   = 128;               end
try, lpf   = inverse.lpf;    catch, lpf   = 1;                 end
try, hpf   = inverse.hpf;    catch, hpf   = 256;               end
try, sdv   = inverse.sdv;    catch, sdv   = 4;                 end
try, Han   = inverse.Han;    catch, Han   = 1;                 end
try, Na    = inverse.Na;     catch, Na    = 1024;              end
try, woi   = inverse.woi;    catch, woi   = [];                end


%==========================================================================
% Spatial parameters
%==========================================================================

% Load Gain or Lead-field matrices
%--------------------------------------------------------------------------
for i = 1:Nl
    gainmat   = D{i}.inv{D{i}.val}.forward.gainmat;
    try
        G     = load(gainmat);
    catch
        [p f] = fileparts(gainmat);
        try
            G = load(f);
        catch
            G = load(fullfile(D{i}.path,f));
        end
    end
    name   = fieldnames(G);
    L{i}   = sparse(getfield(G, name{1}));
    L{i}   = spm_cond_units(L{i});
    Nc(i)  = size(L{i},1);                         % number of channels
end


% assume radii are the same for all VOI
%--------------------------------------------------------------------------
Nd    = size(L{1},2);                              % number of dipoles
Nv    = size(xyz,1);                               % number of VOI
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end

% Restrict source space
%--------------------------------------------------------------------------
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
Ns    = length(Is);

 
% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================

fprintf('Computing Green''s function from graph Laplacian:')
%--------------------------------------------------------------------------
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


% Project to channel modes (U)
%--------------------------------------------------------------------------
Nmax  = Nm;
U{1}  = spm_svd(L{1}*L{1}',exp(-16));
Nm    = min(size(U{1},2),Nmax);
U{1}  = U{1}(:,1:Nm);
UL{1} = U{1}'*L{1};
for i = 2:Nl
    U{i}  = spm_svd(L{i}*L{i}',exp(-16));
    Nm(i) = min(size(U{i},2),Nmax);
    U{i}  = U{i}(:,1:Nm(i));
    UL{i} = U{i}'*L{i};
end
fprintf('Using %i spatial modes\n',Nm)

%==========================================================================
% Temporal parameters
%==========================================================================
Nt    = length(trial);
G     = UL{1};
AY    = {};
for i = 1:Nl

    % Time-window of interest
    %----------------------------------------------------------------------
    if isempty(woi)
        w{i} = round([-D{i}.events.start D{i}.events.stop]*1000/D{i}.Radc);
    else
        w{i} = woi;
    end
    It{i}   = round(w{i}*(D{i}.Radc/1000)) + D{i}.events.start;
    It{i}   = max(1,It{i}(1)):min(It{i}(end),size(D{i}.data,2));

    % Peri-stimulus time
    %----------------------------------------------------------------------
    pst{i}  = (It{i} - D{i}.events.start - 1);
    pst{i}  = pst{i}/D{i}.Radc*1000;               % peristimulus time (ms)
    dur     = (pst{i}(end) - pst{i}(1))/1000;      % duration (s)
    dct{i}  = (It{i} - It{i}(1))/2/dur;            % DCT frequenices (Hz)
    Nb(i)   = length(It{i});                       % number of time bins

    % Serial correlations
    %----------------------------------------------------------------------
    K       = exp(-(pst{i} - pst{i}(1)).^2/(2*sdv^2));
    K       = toeplitz(K);
    qV{i}   = sparse(K*K');

    % Confounds and temporal subspace
    %----------------------------------------------------------------------
    T       = spm_dctmtx(Nb(i),Nb(i));
    j       = (dct{i} > lpf) & (dct{i} < hpf);
    T       = T(:,j);
    dct{i}  = dct{i}(j);

    % get data (with temporal filtering)
    %======================================================================
    Ic{i} = setdiff(D{i}.channels.eeg, D{i}.channels.Bad);
    for j = 1:Nt
        Y{i,j} = sparse(0);
        if isfield(D{i}.events,'reject')
            c = find(D{i}.events.code == trial(j) & ~D{i}.events.reject);
        else
            c = find(D{j}.events.code == trial(j));
        end
        Ne    = length(c);
        for k = 1:Ne
            Y{i,j} = Y{i,j} + squeeze(D{i}.data(Ic{i},It{i},c(k)))*T/Ne;
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
    [V u]  = spm_svd(WY,1);                      % temporal modes
    Nr(i)  = min(size(v,2),8);                   % number of temporal modes
    V      = V(:,      1:Nr(i));
    u      = u(1:Nr(i),1:Nr(i));
    VE(i)  = sum(sum(u.^2))/sum(sum(WY.^2));     % variance explained
    S{i}   = T*V;                                % temporal projector
    iV{i}  = inv(S{i}'*qV{i}*S{i});              % precision (mode)
    Vq{i}  = S{i}*iV{i}*S{i}';                   % precision (time)

    % temproal projectors (adjusting for different PSTs)
    %======================================================================
    B{i}  = V*S{i}'*pinv(full(S{i}'))*S{1}';
    for j = 1:Nt
        AY{i,j}  = U{i}'*Y{i,j}*B{i};
    end

    % spatially adjusted and temporally whitened data
    %----------------------------------------------------------------------
    for j = 1:Nt
        Y{i,j}      = Y{i,j}*v;
        AY{end + 1} = A{i}*Y{i,j}*sqrtm(iV{i});
    end
end

fprintf('Using %i temporal modes\n',Nr)
fprintf('accounting for %0.2f percent variance\n',full(100*VE))

% adjusted and sample covariance (whitened)
%----------------------------------------------------------------------
for j = 1:Nt
    UY{1,j}  = spm_cat(AY(:,j))*sqrtm(Vq{1});
end
[Y scale] = spm_cond_units(Y);
UY        = spm_cat(UY)*scale;
YY        = UY*UY';
G         = spm_cat(UL(:));

% create sensor components (Qe)
%==========================================================================
Qe    = {};
for i = 1:Nl
    QE{i,i} = sparse(Nm(i),Nm(i));
end
for i = 1:Nl
    Q      = QE;
    Q{i,i} = U{i}'*U{i};
    Qe{end + 1} = spm_cat(Q);
end
Ne    = length(Qe);

 
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
        MVB   = spm_mvb(UY,G,[],Q,Qe,16);

        % Spatial priors (QP); eliminating minor patterns
        %------------------------------------------------------------------
        pV    = diag(MVB.Cp);
        for i = 1:8
            j = find(pV > 2^i*(max(pV)/256));
            if length(j) < 128
                break
            end
        end
        qp    = Q(:,j)*MVB.Cp(j,j)*Q(:,j)';

        % Accmulate empirical priors
        %------------------------------------------------------------------
        QP{end + 1} = qp;

end
 
switch(type)
    
    case {'MSP','ARD','IID','MMN','LOR','COH'}
 
    % or ReML - ARD
    %----------------------------------------------------------------------
    qp          = sparse(0);
    Q           = {Qe{:} LQpL{:}};
    [Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,sum(Nr)*Nt);
 
    % Spatial priors (QP)
    %----------------------------------------------------------------------
    Ne    = length(Qe);
    Np    = length(Qp);
    hp    = h([1:Np] + Ne);
    for i = 1:Np
        if hp(i) > max(hp)/128;
            try
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            catch
                qp  = qp + hp(i)*Qp{i};
            end
        end
    end
    
    % Accmulate empirical priors
    %----------------------------------------------------------------------
    QP{end + 1} = qp;
    
end

% re-estimate (all together)
%==========================================================================

 



        
    % using spatial priors from group analysis
    %----------------------------------------------------------------------
    LQpL  = {};
    Np    = length(QP);
    for j = 1:Np
        LQpL{j}  = G*QP{j}*G';
    end
    Q     = {Qe{:} LQpL{:}};
 
    % re-do ReML
    %----------------------------------------------------------------------
    [Cy,h,Ph,F] = spm_reml_sc(YY,[],Q,Nb(1)*Nt);
 
    % Covariances: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    Qp    = sparse(0);
    hp    = h([1:Np] + Ne);
    for j = 1:Np
        Qp = Qp + hp(j)*QP{j};
    end
    LCp   = G*Qp;
 
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
    for j = 1:Nt
        
        % trial-type specific source reconstruction
        %------------------------------------------------------------------
        J{j}    = M*UY{j};
 
        % sum of squares
        %------------------------------------------------------------------
        SSR   = SSR + sum(var((UY{j} - G*J{j}),0,2));
        SST   = SST + sum(var(UY{j},0,2));
 
    end
 
    % assess accuracy; signal to noise (over sources)
    %======================================================================
    R2   = 100*(SST - SSR)/SST;
    fprintf('Percent variance explained %.2f (%.2f)\n',R2,R2*mean(VE(i)))
 
    % Save results
    %======================================================================
    inverse.type   = type;                 % inverse model
    inverse.smooth = s;                    % smoothness (0 - 1)
    inverse.xyz    = xyz;                  % VOI (XYZ)
    inverse.rad    = rad;                  % VOI (radius)
    inverse.scale  = scale;                % scalefactor

    inverse.M      = M;                    % MAP projector (reduced)
    inverse.J      = J;                    % Conditional expectation
    inverse.Y      = UY;                   % ERP data (reduced)
    inverse.L      = G;                    % Lead-field (reduced)
    inverse.R      = speye(Nc(1),Nc(1));   % Re-referencing matrix
    inverse.qC     = Cq;                   % spatial covariance
    inverse.qV     = Vq{1};                % temporal correlations
    inverse.T      = speye(Nb(1),Nb(1));   % temporal subspace
    inverse.U      = spm_cat(diag(U));     % spatial subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It{1};                % Indices of time bins
    inverse.Ic     = Ic{1};                % Indices of good channels
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.pst    = pst{1};               % peristimulus time
    inverse.dct    = dct{1};               % frequency range
    inverse.F      = F;                    % log-evidence
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.VE     = mean(VE(i));          % variance explained
    inverse.woi    = w{1};                 % time-window inverted
    
    % save in struct
    %----------------------------------------------------------------------
    DD     = D{1};
    DD.val = 1;
    DD.inv{1}.inverse = inverse;
    DD.inv{1}.method  = 'Imaging';
    
    % and delete old contrasts
    %----------------------------------------------------------------------
    try
        D{i}.inv{D{i}.val} = rmfield(D{i}.inv{D{i}.val},'contrast');
    end
 
    % display
    %======================================================================
    spm_eeg_invert_display(D{i});
    drawnow
 

 
if length(D) == 1, D = D{1}; end
return
