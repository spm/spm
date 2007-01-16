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
%     D.inverse.xyzr   - x,y,z and radius of [bilateral] spherical VOI
%     D.inverse.type   - 'MSP' multiple sparse priors
%                        'LOR' LORETA-like model
%                        'IID' LORETA and minimum norm
%
% Evaluates:
% 
% model.inverse.M      - MAP projector
% model.inverse.J      - Conditional expectation
% model.inverse.L      - Lead field (reduced)
% model.inverse.R      - Rre-referencing matrix
% model.inverse.qC     - spatial  covariance
% model.inverse.qV     - temporal correlations
% model.inverse.T      - temporal subspace
% model.inverse.U      - spatial  subspace
% model.inverse.Is     - Indices of active dipoles
% model.inverse.Nd     - number of dipoles
% model.inverse.Nt     - numner of trials
% model.inverse.pst    - pers-stimulus time
% model.inverse.F      - log-evidence
% model.inverse.R2     - variance accounted for (%)
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
Na    = 1024;                                    % number of active sources
Nd    = Ns;                                      % number of dipoles
Nm    = 64;                                      % number of spatial modes

% defaults
%--------------------------------------------------------------------------
try, type = model.inverse.type;   catch, type = 'MSP';       end
try, con  = model.inverse.con;    catch, con  = 1;           end
try, s    = model.inverse.smooth; catch, s    = 0.4;         end
try, Np   = model.inverse.Np;     catch, Np   = 64;          end
try, xyzr = model.inverse.xyzr;   catch, xyzr = [0 0 0 128]; end

% Peri-stimulus time
%--------------------------------------------------------------------------
k     = 1:Nb;                                    % bins
pst   = (k - D.events.start)/D.Radc*1000;        % peristimulus time (ms}
dur   = (pst(end) - pst(1))/1000;                % duration (s)
dct   = (k - 1)/2/dur;                           % DCT frequenices (Hz)

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

% get data
%==========================================================================

% single-trial analysis (or single ERP)
%--------------------------------------------------------------------------
j     = setdiff(D.channels.eeg, D.channels.Bad);
c     = find(D.events.code == D.events.types(con));
Nt    = length(c);
iV    = inv(T'*qV*T);
Y     = sparse(0);
YY    = sparse(0);
for i = 1:Nt
    Yi  = squeeze(D.data(j,k,c(i)))*T;
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

% Re-reference matrix (R)
%--------------------------------------------------------------------------
[i j] = min(var(Y,0,2));                         % minimum variance channel
R     = speye(Nc,Nc) - sparse(1:Nc,j,1,Nc,Nc);   % re-referencing matrix
YY    = R*YY*R';
Y     = R*Y;
L     = R*L;

% Project from channel space to modes (U)
%--------------------------------------------------------------------------
U     = spm_svd(Y,0);
U     = U(:,1:Nm);
YY    = U'*YY*U;
Y     = U'*Y;
L     = U'*L;


% covariance components
%==========================================================================

% sensor noise (accommodating re-reference)
%--------------------------------------------------------------------------
Qe{1} = U'*R*R'*U;

% Compute spatial coherence: Diffusion on a normalised graph Laplacian GL
%==========================================================================

if ~strcmp(type,'IID')
    
    fprintf('Computing Green''s function from graph Laplacian:')
    %----------------------------------------------------------------------
    vert  = model.mesh.tess_mni.vert;
    face  = model.mesh.tess_mni.face;
    A     = spm_eeg_inv_meshdist(vert,face,0);
    GL    = A - spdiags(sum(A,2),0,Nd,Nd);
    GL    = GL*s/2;
    Qi    = speye(Nd,Nd);
    QL    = sparse(Nd,Nd);
    for i = 1:8
        QL = QL + Qi;
        Qi = Qi*GL/i;
    end
    clear Qi
    QL    = QL*QL;
    QL    = QL.*(QL > exp(-8));
    QL    = spm_cov2corr(QL);
    fprintf(' - done\n')
end

switch(type)

    case {'MSP'}

        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Np    = fix(Np/2);
        Qp    = {};
        LQpL  = {};
        Ip    = sum([vert(:,1) + abs(xyzr(1)), ...
                     vert(:,2) - xyzr(2), ...
                     vert(:,3) - xyzr(3)].^2,2) < xyzr(4)^2;
        Ip    = find(vert(:,1) < 0 & Ip);
        Ip    = Ip(ceil([1:Np]*length(Ip)/Np));
        for i = 1:Np
            
            % left hemisphere
            %--------------------------------------------------------------
            q               = QL(:,Ip(i));
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
            % right hemisphere
            %--------------------------------------------------------------
            [d j] = min(sum([vert(:,1) + vert(Ip(i),1), ...
                             vert(:,2) - vert(Ip(i),2), ...
                             vert(:,3) - vert(Ip(i),3)].^2,2));
            q               = QL(:,j);
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = L*q;
            
            % bilateral
            %--------------------------------------------------------------
            q               = QL(:,Ip(i)) + QL(:,j);
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
        Qp{2}   = QL;
        LQpL{2} = L*Qp{2}*L';


    case {'IID'}

        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = L*L';

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
[Cy,h,Ph,F] = spm_sp_reml(YY,[],Q,Nr*Nt,hE);


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

% select the most 'energetic' dipoles and store their indices in Is
%--------------------------------------------------------------------------
[i j] = sort(-sum((M*Y).^2,2));
Is    = j(1:Na);
M     = M(Is,:);
L     = L(:,Is);

% conditional expectation (of the sum over trials)
%--------------------------------------------------------------------------
J     = M*Y;

% conditional covariance (leading diagonal)
fprintf('Computing conditional covariances:')
% Cq    = Cp - Cp*L'*iC*L*Cp;
%--------------------------------------------------------------------------
Cp    = sparse(Na,1);
for i = 1:Np
    try
        Cp = Cp + hp(i)*sum(Qp{i}.q(Is,:).^2,2);
    catch
        Cp = Cp + hp(i)*diag(Qp{i}(Is,Is));
    end
end
Cq    = Cp - sum(LCp(:,Is).*(iC*LCp(:,Is)))';


fprintf(' - done\n')

% Assess accuracy; signal to noise (over sources), SSE and log-evidence
%==========================================================================
SSR   = sum(var((Y - L*J),0,2));
SST   = sum(var(Y,0,2));
R2    = 100*(SST - SSR)/SST;
fprintf('Variance explained %.2f (percent)\n',R2)

% Save results
%==========================================================================
model.inverse.con    = con;                  % trial or condition
model.inverse.type   = type                  % inverse model
model.inverse.smooth = s;                    % smoothness (0 - 1)
model.inverse.xyzr   = xyzr;                 % VOI (XYZ, Radius)

model.inverse.M      = M/Lscal;              % MAP projector
model.inverse.J      = J*(Yscal/Lscal);      % Conditional expectation
model.inverse.L      = L*Lscal;              % Lead field (reduced)
model.inverse.R      = R;                    % Rre-referencing matrix
model.inverse.qC     = Cq*(Yscal/Lscal)^2;   % spatial  covariance
model.inverse.qV     = qV;                   % temporal correlations
model.inverse.T      = T;                    % temporal subspace
model.inverse.U      = U;                    % spatial  subspace
model.inverse.Is     = Is;                   % Indices of active dipoles
model.inverse.Nd     = Nd;                   % number of dipoles
model.inverse.Nt     = Nt;                   % numner of trials
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
ql    = QL(:,Ip(64));
il    = find(ql > 1/64);
subplot(2,1,1)
spm_mip(ql(il),vert(il,:)',6);
