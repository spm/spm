function [D] = spm_eeg_invert(D)
% ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the 
% spatio-temporal hierarchy implicit in EEG data
% Requires:
% D.inv{val}.inverse:
%
%     inverse.trials - indices of D.events.types to invert
%     inverse.smooth - smoothness of source priors (0 to 1)
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.type   - 'MSP' multiple sparse priors
%                        'LOR' LORETA-like model
%                        'IID' LORETA and minimum norm
%     inverse.xyz    - (n x 3) locations of spherical VOIs
%     inverse.rad    - radius (mm) of VOIs
%
% Evaluates:
% 
%     inverse.M      - MAP projector
%     inverse.J      - Conditional expectation
%     inverse.L      - Lead field (reduced)
%     inverse.R      - Rre-referencing matrix
%     inverse.qC     - spatial  covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal subspace
%     inverse.U      - spatial  subspace
%     inverse.Is     - Indices of active dipoles
%     inverse.Nd     - number of dipoles
%     inverse.Nt     - numner of trials
%     inverse.pst    - pers-stimulus time
%     inverse.F      - log-evidence
%     inverse.R2     - variance accounted for (%)
%__________________________________________________________________________


% D - SPM data structure
%==========================================================================
model = D.inv{D.val};

% defaults
%--------------------------------------------------------------------------
try, trial = model.inverse.trials; catch, trial = 1;           end
try, type  = model.inverse.type;   catch, type  = 'MSP';       end
try, s     = model.inverse.smooth; catch, s     = 0.6;         end
try, Np    = model.inverse.Np;     catch, Np    = 128;         end
try, xyz   = model.inverse.xyz;    catch, xyz   = [0 0 0];     end
try, rad   = model.inverse.rad;    catch, rad   = 128;         end


% Load Gain or Lead field matrix
%--------------------------------------------------------------------------
L     = load(model.forward.gainmat);
name  = fieldnames(L);
L     = sparse(getfield(L, name{1}));

% parameters
%==========================================================================
Nb    = D.Nsamples;                              % number of time bins
Nc    = size(L,1);                               % number of channels
Na    = 1024;                                    % number of active sources
Nd    = size(L,2);                               % number of dipoles
Nm    = 64;                                      % number of spatial modes
Nv    = size(xyz,1);                              % number of VOI

% assume radii are the same for all VBOI
%--------------------------------------------------------------------------
if length(rad) ~= Nv
    rad = rad(1)*ones(Nv,1);
else
    rad = rad(:);
end

% Peri-stimulus time
%--------------------------------------------------------------------------
It    = 1:Nb;                                    % bins
pst   = (It - D.events.start)/D.Radc*1000;       % peristimulus time (ms}
dur   = (pst(end) - pst(1))/1000;                % duration (s)
dct   = (It - 1)/2/dur;                          % DCT frequenices (Hz)

% Confounds and temporal subspace
%--------------------------------------------------------------------------
T     = spm_dctmtx(Nb,Nb);                       % pst subspace
T     = T(:,find((dct >= 1) & (dct <= 64)));     % 2 - 64 Hz
Nr    = size(T,2);                               % number of temporal modes
Nm    = min(Nm,Nr);

% Serial correlations
%--------------------------------------------------------------------------
K     = exp(-(pst - pst(1)).^2/(2*16));
K     = toeplitz(K);
qV    = sparse(K*K');
qV    = qV/mean(diag(qV));
iV    = inv(T'*qV*T);

% get data
%==========================================================================

% single-trial analysis (or single ERP)
%--------------------------------------------------------------------------
Ic    = setdiff(D.channels.eeg, D.channels.Bad);
YY    = sparse(0);
for i = trial
    if isfield(D.events,'reject')
          c = find(D.events.code == D.events.types(i) & ~D.events.reject);
    else
 	  c = find(D.events.code == D.events.types(i));
    end
    Nt(i) = length(c);
    Y{i}  = sparse(0); 
    for j = 1:Nt(i)
        Yi   = squeeze(D.data(Ic,It,c(j)))*T;
        Y{i} = Y{i} + Yi;
        YY   = YY + Yi*iV*Yi';
    end
end


% Re-reference matrix (R)
%--------------------------------------------------------------------------
[i j] = min(diag(YY));                           % minimum variance channel
R     = speye(Nc,Nc) - sparse(1:Nc,j,1,Nc,Nc);   % re-referencing matrix
YY    = R*YY*R';
L     = R*L;

% Project from channel space to modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(YY,0);
U     = U(:,1:Nm);
YY    = U'*YY*U;
L     = U'*L;

% Restrict source space
%==========================================================================
vert  = model.mesh.tess_mni.vert;
face  = model.mesh.tess_mni.face;
Is    = sparse(Nd,1);
for i = 1:Nv
        Iv    = sum([vert(:,1) - xyz(i,1), ...
                     vert(:,2) - xyz(i,2), ...
                     vert(:,3) - xyz(i,3)].^2,2) < rad(i)^2;
        Is    = Is | Iv;
end
Is    = find(Is);
vert  = vert(Is,:);
L     = L(:,Is);
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
    QG    = QG(Is,Is);
    QG    = QG*QG;
    QG    = QG.*(QG > exp(-8));
    QG    = spm_cov2corr(QG);
    fprintf(' - done\n')
end


% covariance components
%==========================================================================

% sensor noise (accommodating re-reference)
%--------------------------------------------------------------------------
Qe{1} = U'*R*R'*U;

% create source compeonts
%--------------------------------------------------------------------------
switch(type)

    case {'MSP'}

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
            LQpL{end + 1}.q = L*q;
            
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                             vert(:,2) - vert(Ip(i),2), ...
                             vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
            % bilateral
            %--------------------------------------------------------------
            q               = QG(:,Ip(i)) + QG(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
        end

    case {'LOR'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = L*L';

        % add smoothness component in source space
        %------------------------------------------------------------------
        Qp{2}   = QG;
        LQpL{2} = L*Qp{2}*L';


    case {'IID'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = L*L';

end

% Inverse solution
%==========================================================================

% ReML - ARD
%--------------------------------------------------------------------------
Q           = {Qe{:} LQpL{:}};
[Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,Nr*sum(Nt));


% Covariances: sensor space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
Ne    = length(Qe);
Np    = length(Qp);
Ce    = sparse(Nm,Nm);
LCp   = sparse(Nm,Ns);
he    = h([1:Ne]);
hp    = h([1:Np] + Ne);
for i = 1:Ne
    try
        Ce = Ce + he(i)*Qe{i}.q*Qe{i}.q';
    catch
        Ce = Ce + he(i)*Qe{i};
    end
end
for i = 1:Np
    try
        LCp = LCp + hp(i)*LQpL{i}.q*Qp{i}.q';
    catch
        LCp = LCp + hp(i)*L*Qp{i};
    end
end

% MAP estimates of instantaneous sources
%==========================================================================
iC    = inv(LCp*L' + Ce);
M     = LCp'*iC;

% select the most 'energetic' dipoles and store their indices in j
%--------------------------------------------------------------------------
[i j] = sort(-sum(M.^2,2));
j     = j(1:Na);
M     = M(j,:);
L     = L(:,j);
Is    = Is(j);

% conditional covariance (leading diagonal)
fprintf('Computing conditional covariances:')
% Cq    = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
Cp    = sparse(Na,1);
for i = 1:Np
    try
        Cp = Cp + hp(i)*sum(Qp{i}.q(j,:).^2,2);
    catch
        Cp = Cp + hp(i)*diag(Qp{i}(j,j));
    end
end
Cq    = Cp - sum(LCp(:,j).*(iC*LCp(:,j)))';


fprintf(' - done\n')

% re-scale and evaluate conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
SSR   = 0;
SST   = 0;
for i = 1:length(Y)
    
    Y{i}  = U'*R*Y{i};
    J{i}  = M*Y{i};

    % sum of squares
    %----------------------------------------------------------------------
    SSR   = SSR + sum(var((Y{i} - L*J{i}),0,2));
    SST   = SST + sum(var(Y{i},0,2));

end

% Assess accuracy; signal to noise (over sources)
%==========================================================================
R2    = 100*(SST - SSR)/SST;
fprintf('Variance explained %.2f (percent)\n',R2)

% Save results
%==========================================================================
model.inverse.type   = type                  % inverse model
model.inverse.smooth = s;                    % smoothness (0 - 1)
model.inverse.xyz    = xyz;                  % VOI (XYZ)
model.inverse.rad    = rad;                  % VOI (rad)

model.inverse.M      = M;                    % MAP projector
model.inverse.J      = J;                    % Conditional expectation
model.inverse.Y      = Y;                    % Sensor data (reduced)
model.inverse.L      = L;                    % Lead field  (reduced)
model.inverse.R      = R;                    % Re-referencing matrix
model.inverse.qC     = Cq;                   % spatial  covariance
model.inverse.qV     = qV;                   % temporal correlations
model.inverse.T      = T;                    % temporal subspace
model.inverse.U      = U;                    % spatial  subspace
model.inverse.Is     = Is;                   % Indices of active dipoles
model.inverse.It     = It;                   % Indices of time bins
model.inverse.Ic     = Ic;                   % Indices of good channels
model.inverse.Nd     = Nd;                   % number of dipoles
model.inverse.Nt     = Nt;                   % number of trials
model.inverse.pst    = pst;                  % pers-stimulus time
model.inverse.F      = F;                    % log-evidence
model.inverse.R2     = R2;                   % variance accounted for (%)


% save in struct
%--------------------------------------------------------------------------
D.inv{D.val}         = model;
D.inv{D.val}.method  = 'Imaging';

% and delete old contrasts
%--------------------------------------------------------------------------
try
    D.inv{D.val}= rmfield(D.inv{D.val},'contrast');
end


% display
%==========================================================================
spm_eeg_invert_display(D);

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
