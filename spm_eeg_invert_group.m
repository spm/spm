function [D] = spm_eeg_invert(D)
% ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the
% spatio-temporal hierarchy implicit in EEG data
% Requires:
% D{i}.inv{val}.inverse:
%
%     inverse.trials - D.events.types to invert
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.Nm     - maximum number of channel modes
%     inverse.type   - 'GS' Greedu serach on MSPs
%                      'MSP' multiple sparse priors
%                      'LOR' LORETA-like model
%                      'IID' LORETA and minimum norm
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%     inverse.lpf    - band-pass filter - low  frequency cutoff (Hz)
%     inverse.hpf    - band-pass filter - high frequency cutoff (Hz)
%     inverse.Lap    - switch for Laplace transform
%     inverse.sdv    - standard devations of Gaussian temporal correlation
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
%     inverse.pst    - pers-stimulus time
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
try, Nm    = inverse.Nm;     catch, Nm    = 96;                end
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
G     = spm_cat(L(:));



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

% get data (with temporal reduction)
%==========================================================================
Nt    = length(trial);
for k = 1:Nl
    Ic    = setdiff(D{k}.channels.eeg, D{k}.channels.Bad);
    for i = 1:Nt
        Y{k,i}  = sparse(0);
        if isfield(D{k}.events,'reject')
            c = find(D{k}.events.code == trial(i) & ~D{k}.events.reject);
        else
            c = find(D{k}.events.code == trial(i));
        end
        for j = 1:length(c)
            Y{k,i} = Y{k,i} + squeeze(D{k}.data(Ic,It,c(j)))*T;
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
YY = W'*spm_cat(Y(:))'*spm_cat(Y(:))*W;

% eliminate unnecessary temporal modes
%--------------------------------------------------------------------------
S     = spm_svd(YY,1/512);
T     = T*S;
Nr    = size(T,2);                               % number of temporal modes
iV    = inv(T'*qV*T);                            % precision (mode space)
Vq    = T*iV*T';
fprintf('Using %i temporal modes\n',Nr)


% Project onto temporal modes (S)
%--------------------------------------------------------------------------
YY    = spm_cat(Y)*kron(speye(Nt,Nt),S*iV*S')*spm_cat(Y)';


% Re-reference matrix (R) (to minimum variance channel)
%--------------------------------------------------------------------------
for i = 1:Nl
    if strcmp(D{i}.modality,'EEG')

        [m j] = min(sum(Y{i}.^2,2));
        R{i}  = speye(Nc(i),Nc(i)) - sparse(1:Nc(i),j,1,Nc(i),Nc(i));
        L{i}  = R{i}*L{i};

    else
        R{i}  = speye(Nc(i),Nc(i));
    end
end
YY    = spm_cat(diag(R))*YY*spm_cat(diag(R))';


% Project to channel modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(G*G');
try
    U = U(:,1:Nm);
end
Nm    = size(U,2);
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
G     = G(:,Is);
Ns    = length(Is);
Na    = min(Na,Ns);

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

% sensor noise (accommodating re-reference)
%--------------------------------------------------------------------------
for i = 1:Nl
    Q{i} = R{i}*R{i}';
end
Qe    = {U'*spm_cat(diag(Q))*U};

% create source compeonts
%--------------------------------------------------------------------------
switch(type)

    case {'MSP','GS'}

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

    case {'LOR'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';

        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = G*Qp{2}*G';


    case {'IID'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = G*G';

end

% Inverse solution
%==========================================================================

% Greedy search over MSPs
%--------------------------------------------------------------------------
if strcmp(type,'GS')

    % extract patterns
    %----------------------------------------------------------------------
    Np    = length(Qp);
    Q     = sparse(Ns,Np);
    for i = 1:Np
        Q(:,i) = Qp{i}.q;
    end

    % Multivariate Bayes
    %----------------------------------------------------------------------
    MVB   = spm_mvb(U'*R*[Y{:}],L,[],Q,Qe,8,1/4);
    M     = MVB.M;
    Cq    = MVB.qC;
    F     = max(MVB.F);

else

    % or ReML - ARD
    %----------------------------------------------------------------------
    Q           = {Qe{:} LQpL{:}};
    [Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,Nr*Nt);
    
end

for i = 1:Nl
    
    % Covariances: sensor space - Ce and source space - L*Cp
    %----------------------------------------------------------------------
    Ne    = length(Qe);
    Np    = length(Qp);
    Ce    = sparse(Nc(i),Nc(i));
    LCp   = sparse(Nc(i),Ns);
    LCpL  = sparse(Nc(i),Nc(i));
    he    = h([1:Ne]);
    hp    = h([1:Np] + Ne);
    for j = 1:Ne
       Ce = Ce + he(j)*speye(Nc(i),Nc(i));
    end
    for j = 1:Np
        try
            LCp  = LCp  + hp(j)*L{i}*Qp{j}.q*Qp{j}.q';
            LCpL = LCpL + hp(j)*L{i}*Qp{j}.q*Qp{j}.q'*L{i}';
        catch
            LCp  = LCp  + hp(j)*L{i}*Qp{j};
            LCpL = LCpL + hp(j)*L{i}*Qp{j}*L{i}';

        end
    end

    % MAP estimates of instantaneous sources
    %======================================================================
    iC    = inv(Ce + LCpL);
    M{i}  = LCp'*iC;

    % conditional covariance (leading diagonal)
    % Cq    = Cp - Cp*L'*iC*L*Cp;
    %----------------------------------------------------------------------
    Cp    = sparse(Ns,1);
    for j = 1:Np
        try
            Cp = Cp + hp(j)*sum(Qp{j}.q.^2,2);
        catch
            Cp = Cp + hp(j)*diag(Qp{j});
        end
    end
    Cq{i} = Cp - sum(LCp.*M{i}')';

    % re-scale and evaluate conditional expectation (of the sum over trials)
    %----------------------------------------------------------------------  
    SSR   = 0;
    SST   = 0;
    for j = 1:Nt

        Y{i,j}  = R{i}*Y{i,j}*S;
        J{i,j}  = M{i}*Y{i,j};

        % sum of squares
        %------------------------------------------------------------------
        SSR   = SSR + sum(var((Y{i,j} - L{i}*J{i,j}),0,2));
        SST   = SST + sum(var(Y{i,j},0,2));
 
    end

    % Assess accuracy; signal to noise (over sources)
    %======================================================================
    R2(i)  = 100*(SST - SSR)/SST;
    fprintf('Variance explained %.2f (percent)\n',R2)

    % Save results
    %======================================================================
    inverse.type   = type                  % inverse model
    inverse.smooth = s;                    % smoothness (0 - 1)
    inverse.xyz    = xyz;                  % VOI (XYZ)
    inverse.rad    = rad;                  % VOI (rad)

    inverse.M      = M{i};                    % MAP projector (reduced)
    inverse.J      = J(i,:);               % Conditional expectation
    inverse.Y      = Y(i,:);               % ERP data (reduced)
    inverse.L      = L{i};                 % Lead-field  (reduced)
    inverse.R      = R{i};                 % Re-referencing matrix
    inverse.qC     = Cq{i};                % spatial  covariance
    inverse.qV     = Vq;                   % temporal correlations
    inverse.T      = T;                    % temporal subspace
    inverse.U      = speye(Nc(i),Nc(i));   % spatial  subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = It;                   % Indices of time bins
    inverse.Ic     = Ic;                   % Indices of good channels
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.pst    = pst;                  % pers-stimulus time
    inverse.dct    = dct;                  % frequency range
    inverse.F      = F;                    % log-evidence
    inverse.R2     = R2(i);                % variance accounted for (%)

    % save in struct
    %--------------------------------------------------------------------------
    D{i}.inv{D{i}.val}.inverse = inverse;
    D{i}.inv{D{i}.val}.method  = 'Imaging';

    % and delete old contrasts
    %--------------------------------------------------------------------------
    try
        D{i}.inv{D{i}.val} = rmfield(D{i}.inv{D{i}.val},'contrast');
    end


    % display
    %==========================================================================
    spm_eeg_invert_display(D{i});
    drawnow

end

if length(D) == 1, D = D{:}; end
return





% NOTES
%==========================================================================

% informed spatial basis functions using a generlised eigemode solution
%--------------------------------------------------------------------------
Cpos  = Y*(T'*diag(pst > 0)*T)*Y' + speye(Nc,Nc)*exp(-8);
Cpre  = Y*(T'*diag(pst < 0)*T)*Y' + speye(Nc,Nc)*exp(-8);
[U S] = spm_svd(inv(Cpos)*(Cpre),0);
U     = U(:,end - Nm + 1:end);

% display a selected basis function
%--------------------------------------------------------------------------
ql    = QG(:,Ip(64));
il    = find(ql > 1/64);
subplot(2,1,1)
spm_mip(ql(il),vert(il,:)',6);
axis image
