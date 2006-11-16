function [D] = spm_eeg_invert(D)
% ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the 
% spatio-temporal hierarchy implicit in EEG data
% Requires:
%
%     D.inverse.con    - condition or trial type
%     D.inverse.smooth - smoothness of source priors (0 to 1)
%     D.inverse.Np     - number of sparse priors
%     D.inverse.type   - 'MSP' multiple sparse priors
%                        'LOR' LORETA-like model
%                        'IID' LORETA and WMN
%__________________________________________________________________________


% D - SPM data structure
%==========================================================================
try
    model = D.inv{D.val};
catch
    model = D.inv{end};
end

% Load Gain or Lead field matrix
%--------------------------------------------------------------------------
j     = setdiff(D.channels.eeg, D.channels.Bad);
L     = load(model.forward.gainmat);
name  = fieldnames(L);
L     = sparse(getfield(L, name{1}));
L     = L(j,:);

% parameters
%==========================================================================
Nb    = D.Nsamples;                              % number of time bins
Nc    = size(L,1);                               % number of channels
Ns    = size(L,2);                               % number of sources
Na    = 1024;                                    % number of active sources
Nd    = Ns;                                      % number of dipoles
Nm    = 64;                                      % number of spatial modes

% defaults
%--------------------------------------------------------------------------
try, type = model.inverse.type;   catch, type = 'MSP'; end
try, con  = model.inverse.con;    catch, con  = 1;     end
try, s    = model.inverse.smooth; catch, s    = 1/2;   end
try, Np   = model.inverse.Np    ; catch, Np   = 128;   end

% Peri-stimulus time
%--------------------------------------------------------------------------
k     = 1:Nb;                                    % bins
pst   = (k - D.events.start)/D.Radc*1000;        % peristimulus time (ms}
dur   = (pst(end) - pst(1))/1000;                % duration (s)
dct   = (k - 1)/2/dur;                           % DCT frequenices (Hz)

% Confounds and temporal subspace
%--------------------------------------------------------------------------
X     = ones(Nc,1);                              % confounds
R     = speye(Nc,Nc) - X*pinv(X);                % re-referencing matrix
T     = spm_dctmtx(Nb,Nb);                       % pst subspace
T     = T(:,find((dct > 4) & (dct < 64)));       % 4 - 64 Hz
Nr    = size(T,2);                               % number of temporal modes


% Get serial correlations from [known] filtering
%--------------------------------------------------------------------------
K     = sparse(1,1,1,Nb,1);
K     = filtfilt(D.filter.para{1},D.filter.para{2},full(K));
K     = toeplitz(K);
Vq    = sparse(K*K');
Vq    = Vq/mean(diag(Vq));
iV    = inv(T'*Vq*T);


% get data
%==========================================================================

% single-trial analysis (or single ERP)
%--------------------------------------------------------------------------
Y     = sparse(0);
YY    = sparse(0);
c     = find(D.events.code == D.events.types(con));
Nt    = length(c);
for i = 1:Nt
    Yi  = R*(squeeze(D.data(j,k,c(i)))*T);
    Y   = Y  + Yi;
    YY  = YY + Yi*iV*Yi';
end


% scale data and sources
%--------------------------------------------------------------------------
Yscal = sqrt(norm(YY,1));
Lscal = sqrt(norm(L*L',1));
YY    = YY/(Yscal*Yscal);
Y     = Y/Yscal;
L     = L/Lscal;

% Project from channel space to modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(YY,0);
U     = U(:,1:Nm);
YY    = U'*YY*U;
Y     = U'*Y;
L     = U'*L;
X     = U'*X;


% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
fprintf('Computing Green''s function from graph Laplacian:')
%--------------------------------------------------------------------------
load(model.mesh.tess_ctx)
A     = spm_eeg_inv_meshdist(vert,face,0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;
Qi    = speye(Nd,Nd);
QL    = sparse(Nd,Nd);
for i = 1:8
    QL = QL + Qi;
    Qi = Qi*GL/i;
end
QL    = QL*QL;
QL    = spm_cov2corr(QL);

fprintf(' - done\n')

switch(type)

    case {'MSP'}

        % create overMSP spatial basis set in source space
        %------------------------------------------------------------------
        Ip    = fix([1:Np]*Ns/Np);
        for i = 1:Np
            Qp{i}.q   = QL(:,Ip(i));
            LQpL{i}.q = L*Qp{i}.q ;
        end

        % create overMSP spatial basis set in source space
        %------------------------------------------------------------------
        Qe{1}.q = speye(Nm,Nm);

    case {'LOR','IID'}

        % create smoothnes component in source space
        %------------------------------------------------------------------
        Qp{1}   = QL;
        LQpL{1} = L*Qp{1}*L';

        % create overMSP spatial basis set in source space
        %------------------------------------------------------------------
        Qe{1} = speye(Nm,Nm);
end

if strcmp(type,'IID')

    % add minimum norm prior
    %----------------------------------------------------------------------
    Qp{2}     = speye(Ns,Ns);
    LQpL{2}   = L*L';

end

% Inverse solution
%==========================================================================
Q     = {Qe{:} LQpL{:}};

% hyperpriors
%--------------------------------------------------------------------------
Ne    = length(Qe);
Np    = length(Qp);
hE    = [sparse(Ne,1); (sparse(Np,1) - 4)];

% ReML
%--------------------------------------------------------------------------
[Cy,h,Ph,F] = spm_sp_reml(YY,X,Q,Nr*Nt,hE);


% Covariances: sensor space - Ce and source space - L*Cp
%--------------------------------------------------------------------------
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

% select the most active dipoles
%--------------------------------------------------------------------------
[i j] = sort(-sum((M*Y).^2,2));
Is    = j(1:Na);
M     = M(Is,:);
L     = L(:,Is);

% conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
J     = M*Y*T';

% conditional covariance
fprintf('Computing conditional covariances:')
%--------------------------------------------------------------------------
Cp    = sparse(Na,Na);
for i = 1:Np
    try
        Cp = Cp + hp(i)*Qp{i}.q(Is,:)*Qp{i}.q(Is,:)';
    catch
        Cp = Cp + hp(i)*Qp{i}(Is,Is);
    end
end
Cq    = Cp - Cp*L'*iC*L*Cp;
Cq    = diag(Cq);

fprintf(' - done\n')

% Assess accuracy; signal to noise (over sources), SSE and log-evidence
%==========================================================================
SSR    = sum(var((Y - L*J*T),0,2));
SST    = sum(var(Y,0,2));
R2     = 100*(SST - SSR)/SST;
fprintf('Variance explained %.2f (percent)\n',R2)

% Save results
%==========================================================================
model.inverse.con    = con;                  % trial or condition
model.inverse.type   = type                  % inverse model
model.inverse.smooth = s;                    % smoothnes (sd; mm)

model.inverse.MAP    = M*U'/Lscal;           % MAP projector
model.inverse.J      = J*(Yscal/Lscal);      % Conditional expectation
model.inverse.qC     = Cq*(Yscal/Lscal)^2;   % Conditional covariance
model.inverse.qV     = Vq;                   % temporal correlations
model.inverse.T      = T;                    % temproal subspace
model.inverse.Is     = Is;                   % Indices of active dipoles
model.inverse.Nd     = Nd;                   % number of dipoles
model.inverse.Nt     = Nt;                   % numner of trials
model.inverse.pst    = pst;                  % pers-stimulus time
model.inverse.F      = F;                    % log-evidence
model.inverse.R2     = R2;                   % variance accounted for (%)

D.inv{D.val}         = model;

if str2num(version('-release')) >= 14
    save(fullfile(D.path,D.fname), '-V6','D');
else
    save(fullfile(D.path,D.fname),'D');
end

% display
%==========================================================================
spm_eeg_invert_display(D);

return

% display a selected basis function
%--------------------------------------------------------------------------
ql    = QL(:,Ip(64));
il    = find(ql > 1/64);
subplot(2,1,1)
spm_mip(ql(il),vert(il,:)',6);
