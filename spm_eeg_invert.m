function [D] = spm_eeg_invert(D)
% ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert(D)
% ReML estimation of regularisation hyperparameters using the 
% spatio-temporal hierarchy implicit in EEG data
% Requires:
%
%     D.inverse.con    - condition or trial type
%     D.inverse.smooth - smoothness of source priors (mm)
%     D.inverse.type   - 'MSP' multiple sparse priors
%                        'LOR' LORETA-like model
%                        'IID' LORETA and WMN
%__________________________________________________________________________


% D - SPM data structure
%==========================================================================
model = D.inv{D.val};

% Load Gain or Lead field matrix
%--------------------------------------------------------------------------
L     = load(model.forward.gainmat);
name  = fieldnames(L);
L     = sparse(getfield(L, name{1}));

% parameters
%==========================================================================
Nb    = D.Nsamples;                              % number of time bins
Nc    = size(L,1);                               % number of channels
Ns    = size(L,2);                               % number of sources
Nd    = Ns;                                      % number of dipoles
Np    = 128;                                     % number of patches
Nm    = 64;                                      % number of spatial modes

% defaults
%--------------------------------------------------------------------------
try, con  = model.inverse.con;    catch, con  = 1;     end
try, type = model.inverse.type;   catch, type = 'MSP'; end
try, s    = model.inverse.smooth; catch, s    = 32;    end

% Peri-stimulus time
%--------------------------------------------------------------------------
k     = 1:Nb;                                    % bins
pst   = (k - D.events.start)/D.Radc*1000;        % peristimulus time (ms}
dct   = k*(pst(end) - pst(1))/1000/2;            % DCT frequenices (Hz)

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
j     = D.channels.eeg;
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


% get geodesic distance matrix (s = smoothness in mm)
%--------------------------------------------------------------------------
Mdist = load(model.mesh.CtxGeoDist);
name  = fieldnames(Mdist);
Mdist = sparse(getfield(Mdist,name{1}));

switch(type)

    case {'MSP'}

        % create overMSP spatial basis set in source space
        %------------------------------------------------------------------
        Ip    = fix([1:Np]*Ns/Np);
        for i = 1:Np
            j         = Ip(i);
            k         = find(Mdist(:,j));
            Qi        = sparse(k,1,exp(-Mdist(k,j).^2/(2*s^2)),Ns,1);
            Qi(j)     = 1;
            Qp{i}.q   = Qi;
            LQpL{i}.q = L*Qi;
        end

        % create overMSP spatial basis set in source space
        %------------------------------------------------------------------
        Qe{1}.q = speye(Nm,Nm);

    case {'LOR','IID'}

        % create smoothnes component in source space
        %------------------------------------------------------------------
        Qp{1} = speye(Ns,Ns);
        [I,J] = find(Mdist);
        for k = 1:length(I)
            if I(k) ~= J(k)
                Qp{1}(I(k),J(k)) = exp(-Mdist(I(k),J(k)).^2/(2*s^2));
            end
        end
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
M     = LCp'*inv(LCp*L' + Ce);

% select the most active dipoles
%--------------------------------------------------------------------------
Ns    = 512;
[i j] = sort(-sum((M*Y).^2,2));
Is    = j(1:Ns);
M     = M(Is,:);
L     = L(:,Is);

% conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
J     = M*Y*T';

% conditional covariance
fprintf('Computing conditional covariances:')
%--------------------------------------------------------------------------
Cp    = sparse(Ns,Ns);
for i = 1:Np
    try
        Cp = Cp + hp(i)*Qp{i}.q(Is,:)*Qp{i}.q(Is,:)';
    catch
        Cp = Cp + hp(i)*Qp{i}(Is,Is);
    end
end
Cq    = inv(L'*inv(Ce)*L + inv(Cp + speye(Ns,Ns)*exp(-16)));
Cq    = diag(Cq);

fprintf(' - done\n')

% display
%==========================================================================
Fgraph   = spm_figure('GetWin','Graphics');
clf(Fgraph)
figure(Fgraph)
tess_ctx = model.mesh.tess_ctx;
load(tess_ctx)


% maximum response and confidence intervals
%--------------------------------------------------------------------------
[i j] = max(max(abs(J),[],2));
Jt    = J(j,:);                      % over time
ci    = Nt*sqrt(Cq(j))*1.64;
[i j] = max(max(abs(J),[],1));
Js    = J(:,j);                      % over sources
Jmax  = abs(sparse(Is,1,Js,Nd,1));

% maximum response - space
%--------------------------------------------------------------------------
subplot(3,2,1)
spm_eeg_inv_render(Jmax,tess_ctx)
view([180 -90])
title('esimated response')

% maximum response - time
%--------------------------------------------------------------------------
subplot(3,2,2)
plot(pst,Jt,pst,Jt + ci,':',pst,Jt - ci,':')
xlabel('time  ms')
title('esimated response')
axis square

% observed and predicted repeonse in sensor sapce
%--------------------------------------------------------------------------
subplot(4,2,3)
imagesc(U*L*J)
title('predicted channel data')
axis image off

subplot(4,2,4)
imagesc(U*Y*T')
title('observed channel data')
axis image off

% PPM
%==========================================================================
subplot(2,1,2)
Z     = abs(Js)./(Nt*sqrt(Cq));
i     = find(Z > spm_invNcdf(.95));

spm_mip(Jmax(Is(i)),vert(Is(i),:)',6);
axis image
title('esimated response')
drawnow

% Assess accuracy; signal to noise (over sources), SSE and log-evidence
%==========================================================================
SSR    = sum(var(Y - L*J*T,0,2));
SST    = sum(var(Y,0,2));
R2     = 100*(SST - SSR)/SST;

% Save results
%==========================================================================
model.inverse.con    = con;
model.inverse.type   = type
model.inverse.smooth = s;


model.inverse.MAP    = M*U'/Lscal;
model.inverse.qC     = Cq*(Yscal/Lscal)^2;
model.inverse.qV     = Vq;
model.inverse.T      = T;
model.inverse.Is     = Is;
model.inverse.Nd     = Nd;
model.inverse.pst    = pst;
model.inverse.F      = F;
model.inverse.R2     = R2;

D.inv{D.val}         = model;

if str2num(version('-release')) >= 14
    save(fullfile(D.path,D.fname), '-V6','D');
else
    save(fullfile(D.path,D.fname),'D');
end
