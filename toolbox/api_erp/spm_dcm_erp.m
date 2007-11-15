function DCM = spm_dcm_erp(DCM)   
% Estimate parameters of a DCM model (Newton's methods)
% FORMAT DCM = spm_dcm_erp(DCM)   
%
% DCM     
%    name: name string
%       Lpos:  Source locations
%       xY:    data   [1x1 struct]
%       xU:    design [1x1 struct]
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
clear spm_erp_L

% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                   catch, DCM.name  = 'DCM_ERP'; end
try, h     = DCM.options.h;      catch, h         = 1;         end
try, nm    = DCM.options.Nmodes; catch, nm        = 8;         end
try, onset = DCM.options.onset;  catch, onset     = 80;        end
try, model = DCM.options.model;  catch, model     = 'ERP';     end


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

% check the number of modes is greater or equal to the number of sources
%--------------------------------------------------------------------------
nm     = max(nm,Nr);

% Re-reference matrix (R)
%--------------------------------------------------------------------------
[i j]  = min(var(xY.y));                        % minimum variance channel
R      = speye(Nc,Nc) - sparse(1:Nc,j,1,Nc,Nc); % re-referencing matrix

% confounds - DCT: T - an idempotent matrix spanning temporal subspace
%--------------------------------------------------------------------------
if h == 0
    X0 = sparse(Ns,1);
else
    X0 = spm_dctmtx(Ns,h);
end
hT = [-12 0];
if hT(2)
    i       = find(DCM.xY.Time > hT(1) & DCM.xY.Time < hT(2));
    Ni      = length(i);
    Nh      = min(Ni,8);
    H0      = sparse(Ns,Nh);
    H0(i,:) = diag(hanning(Ni))*spm_dctmtx(Ni,Nh);
    X0      = [X0 H0];
    
end

warning off
Ti     = speye(Ns) - X0*inv(X0'*X0)*X0';         % null space of confounds
T      = kron(speye(Nt,Nt),Ti);
xY.X0  = kron(speye(Nt,Nt),X0);
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
M       = DCM.M;
try, M  = rmfield(M,'g'); end

% adjust onset relative to pst
%--------------------------------------------------------------------------
M.dur   = xU.dur;
M.ons   = onset - xY.Time(1);

switch lower(model)
    
    % linear David et al model (linear in states)
    %======================================================================
    case{'erp'}

        % prior moments on parameters
        %--------------------------------------------------------------------------
        [pE,gE,pC,gC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,M.dipfit,length(M.ons));

        % inital states and equations of motion
        %--------------------------------------------------------------------------
        M.x  =  spm_x_erp(pE);
        M.f  = 'spm_fx_erp';
        M.G  = 'spm_lx_erp';
        
    % linear David et al model (linear in states) - fast version for SEPs
    %======================================================================
    case{'sep'}

        % prior moments on parameters
        %--------------------------------------------------------------------------
        [pE,gE,pC,gC] = spm_sep_priors(DCM.A,DCM.B,DCM.C,M.dipfit,length(M.ons));

        % inital states
        %--------------------------------------------------------------------------
        M.x  = spm_x_erp(pE);
        M.f  = 'spm_fx_erp';
        M.G  = 'spm_lx_sep';
        
    % Neural mass model (nonlinear in states)
    %======================================================================
    case{'nmm'}

        % prior moments on parameters
        %--------------------------------------------------------------------------
        [pE,gE,pC,gC] = spm_nmm_priors(DCM.A,DCM.B,DCM.C,M.dipfit,length(M.ons));

        % inital states
        %--------------------------------------------------------------------------
        M.x  = spm_x_nmm(pE);
        M.f  = 'spm_fx_nmm';
        M.G  = 'spm_lx_nmm';      
        
    otherwise
        warndlg('Unknown model')
end

% likelihood model
%--------------------------------------------------------------------------
M.FS  = 'spm_fy_erp';
M.IS  = 'spm_int_U';
M.fu  = 'spm_erp_u';
M.pE  = pE;
M.pC  = pC;
M.gE  = gE;
M.gC  = gC;
M.m   = nu;
M.n   = length(spm_vec(M.x));
M.l   = Nc;
M.ns  = Ns*Nt;
M.E   = U;


% EM: inversion
%--------------------------------------------------------------------------
[Qp,Qg,Cp,Cg,Ce,F] = spm_nlsi_N(M,xU,xY);

% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on

% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
L   = feval(M.G, Qg,M);           % get gain matrix
x   = feval(M.IS,Qp,M,xU);        % prediction (source space)
y   = x*L';                       % prediction (sensor space)
r   = T*(xY.y - y);               % prediction error
y   = y*M.E*M.E';                 % remove spaital confounds
r   = r*M.E*M.E';                 % remove spaital confounds
x   = x(:,find(any(L)));

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
DCM.M  = M;                    % model specification
DCM.xY = xY;                   % data structure
DCM.xU = xU;                   % input structure
DCM.Ep = Qp;                   % conditional expectation f(x,u,p)
DCM.Cp = Cp;                   % conditional covariances G(g)
DCM.Eg = Qg;                   % conditional expectation
DCM.Cg = Cg;                   % conditional covariances
DCM.Pp = Pp;                   % conditional probability
DCM.H  = H;                    % conditional responses (y), projected space
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.K  = K;                    % conditional responses (x)
DCM.R  = E;                    % conditional residuals (y)
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Ce = Ce;                   % ReML error covariance
DCM.F  = F;                    % Laplace log evidence

DCM.options.h      = h;
DCM.options.Nmodes = nm;
DCM.options.onset  = onset;
DCM.options.model  = model;

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
        G(M.dipfit.Ip{i},i) = M.dipfit.U{i}*Qg.L(:,i);
        try
            G(M.dipfit.Ip{i},i + Nr) = M.dipfit.U{i}*Qg.L(:,i + Nr);
        end
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
        Y{i} = U'*xY.xy{i}'*Ti;
    end

    inverse.trials = DCM.options.trials;   % trial or condition
    inverse.type   = 'DCM';                % inverse model

    inverse.J      = J;                    % Conditional expectation
    inverse.L      = L;                    % Lead field (reduced)
    inverse.R      = R;                    % Re-referencing matrix
    inverse.T      = Ti;                   % temporal subspace
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

    % append DCM results and save in struct
    %----------------------------------------------------------------------
    try, val = DCM.val;  catch, val = 1; end
    D        = spm_eeg_ldata(DCM.xY.Dfile);
    
    D.inv{end + 1}      = D.inv{val};
    D.inv{end}.date     = date;
    D.inv{end}.comment  = {'DCM'};
    D.inv{end}.inverse  = inverse;
    D.val               = length(D.inv);
    try
        D.inv{end}      = rmfield(D.inv{end},'contrast');
    end
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
assignin('base','DCM',DCM)
return
