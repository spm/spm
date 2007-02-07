function DCM = spm_dcm_erp(DCM)   
% Estimate parameters of a DCM model
% FORMAT DCM = spm_dcm_erp(DCM)   
%
% DCM     
%    name: name string
%       M:  Forward model
%              M.dipfit - leadfield specification
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%       A: {[nr x nr double]  [nr x nr double]  [nr x nr double]}
%       B: {[nr x nr double], ...}   Connection constraints
%       C: [nr x 1 double]
%
%   options.Tdcm         - [start end] time window in ms
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.h            - number of DCT drift terms (usually 1 or 2)
%   options.Nmodes       - number of spatial models to invert
%   options.type         - 1 - 'ECD (EEG)'
%                          2 - 'ECD (MEG)'
%                          3 - 'Imaging'
%                          4 - 'LFP' 
%                          (see spm_erp_L)
%__________________________________________________________________________
% %W% Karl Friston %E%


% check options 
%==========================================================================

% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                   catch, DCM.name = 'DCM_ERP'; end
try, h     = DCM.options.h;      catch, h        = 1;         end
try, nm    = DCM.options.Nmodes; catch, nm       = 8;         end
try, onset = DCM.options.onset;  catch, onset    = 60;        end



% Data and spatial model
%==========================================================================
DCM    = spm_dcm_erp_data(DCM);
DCM    = spm_dcm_erp_dipfit(DCM);
xY     = DCM.xY;
try
    xU = DCM.U;
catch
    xU = DCM.xU;
end

% dimensions
%--------------------------------------------------------------------------
Nt     = length(xY.xy);                 % number of trials
Nr     = length(DCM.A{1});              % number of sources
Nc     = size(xY.xy{1},2);              % number of channels
Ns     = size(xY.xy{1},1);              % number of time bins
nu     = size(xU.X,2);                  % number of inputs
nx     = Nr*9 + 1;                      % number of states

% check the number of modes is greater or equal to the number of sources
%--------------------------------------------------------------------------
nm     = max(nm,Nr);

% Re-reference matrix (R)
%--------------------------------------------------------------------------
[i j]  = min(var(xY.y));                        % minimum variance channel
R      = speye(Nc,Nc) - sparse(1:Nc,j,1,Nc,Nc); % re-referencing matrix

% confounds - DCT: T - an idemptoent matrix spans the interesting temporal
% space
%--------------------------------------------------------------------------
if h == 0
    X0 = sparse(Ns,1);
else
    X0 = spm_dctmtx(Ns,h);
end
warning off
X0     = kron(speye(Nt,Nt),X0);
T      = speye(Ns*Nt) - X0*inv(X0'*X0)*X0';  % null space of confounds
xY.X0  = X0;
warning on

% Feature selection using principal components (U) of channel space
%--------------------------------------------------------------------------
try
    U  = M.S(:,[1:nm]);  % use previous basis (S) if specified
catch
    y  = T*xY.y*R';
    U  = spm_svd(y'*y,0);
    U  = U(:,[1:nm]);
end

% assume noise variance is the same over modes
%--------------------------------------------------------------------------
xY.Q   = {kron(speye(nm),kron(speye(Nt),speye(Ns)))};


% Inputs
%==========================================================================

% trial effects
%--------------------------------------------------------------------------
try
    if size(xU.X,2) - length(DCM.B)
        warndlg({'please ensure number of trial specific effects', ...
            'encoded by DCM.xU.X & DCM.B are the same'})
    end
catch
    DCM.B = {};
end
try
    xU.u  = kron(xU.X,ones(Ns,1));
catch
    xU.u  = sparse(Nt*Ns,0);
end

% stimulus parameters
%--------------------------------------------------------------------------
xU.dt   = xY.dt;
xU.dur  = xU.dt*(Ns - 1);

% model specification and nonlinear system identification
%==========================================================================

% make model specification a global variable
%--------------------------------------------------------------------------
global M; M = DCM.M;

% adjust onset relative to pst
%--------------------------------------------------------------------------
dur     = xU.dur;
ons     = onset - xY.Time(1);

% prior moments
%--------------------------------------------------------------------------
[pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,M.dipfit,dur,ons);

% likelihood model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_erp';
M.g   = 'spm_gx_erp';
M.IS  = 'spm_int_U';
M.FS  = 'y*E';
M.x0  = 'spm_x_erp';
M.x   = sparse(nx,1);
M.pE  = pE;
M.pC  = pC;
M.m   = nu;
M.n   = nx;
M.l   = Nc;
M.ns  = Ns*Nt;
M.E   = U;

% EM: inversion
%--------------------------------------------------------------------------
[Qp,Cp,Ce,F] = spm_nlsi_GN(M,xU,xY);

% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on

% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
type          = M.dipfit.type;
M.dipfit.type = 'LFP';
x             = feval(M.IS,Qp,M,xU);        % prediction (source space)
M.dipfit.type = type;                       % reset type
y             = feval(M.IS,Qp,M,xU);        % prediction (sensor space)
r             = T*(xY.y - y);               % prediction error
y             = y*M.E*M.E';                 % remove spaital confounds
r             = r*M.E*M.E';                 % remove spaital confounds

% trial specific respsonses (in mode, channel and source space)
%--------------------------------------------------------------------------
for i = 1:Nt
    j     = [1:Ns] + (i - 1)*Ns;
    H{i}  = y(j,:)*M.E;
    E{i}  = r(j,:)*M.E;
    Hc{i} = y(j,:);
    Ec{i} = r(j,:);
    K{i}  = x(j,:);
end

% store estimates in DCM
%--------------------------------------------------------------------------
DCM.M    = M;                    % model specification
DCM.xY   = xY;                   % data structure
DCM.xU   = xU;                   % input structure
DCM.Ep   = Qp;                   % conditional expectation
DCM.Cp   = Cp;                   % conditional covariances
DCM.Pp   = Pp;                   % conditional probability
DCM.H    = H;                    % conditional responses (y), projected space
DCM.Hc   = Hc;                   % conditional responses (y), channel space
DCM.K    = K;                    % conditional responses (x)
DCM.R    = E;                    % conditional residuals (y)
DCM.Rc   = Ec;                   % conditional residuals (y), channel space
DCM.Ce   = Ce;                   % ReML error covariance
DCM.F    = F;                    % Laplace log evidence

% store estimates in D
%--------------------------------------------------------------------------
if strcmp(M.dipfit.type,'Imaging')
    
    
    % Assess accuracy; signal to noise (over sources), SSE and log-evidence
    %---------------------------------------------------------------------
    SSR   = sum(var((T*xY.y*U - T*y*U)));
    SST   = sum(var(T*xY.y*U));
    R2    = 100*(SST - SSR)/SST;
    
    
    % reconsttuct sources in dipole space
    %----------------------------------------------------------------------
    Nd    = M.dipfit.Nd;
    G     = sparse(Nd,Nr);
    for i = 1:Nr
        G(M.dipfit.Ip{i},i) = M.dipfit.U{i}*Qp.L(:,i);
    end
    Is    = find(any(G,2));
    G     = G(Is,:);
    for i = 1:Nt
        J{i} = G*K{i}';
    end

    % get dipole space lead field
    %----------------------------------------------------------------------
    L     = load(M.dipfit.gainmat);
    name  = fieldnames(L);
    L     = sparse(getfield(L, name{1}));
    L     = U'*L(:,Is);
    
    % reduced data (for each trial
    %----------------------------------------------------------------------
    for i = 1:Nt
        Y{i} = U'*xY.xy{i}'*T;
    end

    inverse.trials = DCM.options.trials;   % trial or condition
    inverse.type   = 'DCM';                % inverse model

    inverse.J      = J;                    % Conditional expectation
    inverse.L      = L;                    % Lead field (reduced)
    inverse.R      = R;                    % Re-referencing matrix
    inverse.T      = T;                    % temporal subspace
    inverse.U      = U;                    % spatial  subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = DCM.xY.It;            % Indices of time bins
    inverse.Ic     = DCM.xY.Ic;            % Indices of good channels
    inverse.Y      = Y;                    % reduced data
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.Nt     = Nt;                   % numner of trials
    inverse.pst    = xY.Time;              % pers-stimulus time
    inverse.F      = DCM.F;                % log-evidence
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.dipfit = M.dipfit;             % forward model for DCM

    % save in struct
    %----------------------------------------------------------------------
    try, val = DCM.val;  catch, val = 1; end
    D        = spm_eeg_ldata(DCM.xY.Dfile);
    D.inv{val}.inverse = inverse;

    if spm_matlab_version_chk('7.1') >= 0
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
        save(fullfile(D.path, D.fname), 'D');
    end

end

% and save
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7.1') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end

return



% % Code for generalised SVD feature selection:
% % Project modelled data onto the principal components of channel space
% %--------------------------------------------------------------------------
% y     = R0*xY.y;
% S     = spm_svd(y'*y);
% 
% % EM: This is a loop that will perfom cnaonical feature selection
% %--------------------------------------------------------------------------
% ML    = exp(-32);
% for i = 1:8
%     
%     % inversion
%     %----------------------------------------------------------------------
%     [Qp,Cp,Ce,F] = spm_nlsi_GN(M,xU,xY);
%     r            = y - R0*feval(M.IS,Qp,M,xU);
%     
%     if F < ML(end), break, end
% 
%     % Canonical feature selection
%     %----------------------------------------------------------------------
%     M.E   = spm_svd(inv(S'*r'*r*S)*(S'*y'*y*S),0);
%     M.E   = orth(full(S*M.E(:,1:nm)));
%     
%     % reset stating estimates
%     %----------------------------------------------------------------------
%     ML(i) = F;
%     M.P   = Qp;
%     
% end


