function DCM = spm_dcm_erp(DCM)   
% Estimate parameters of a DCM model
% FORMAT DCM = spm_dcm_erp(DCM)   
%
% DCM     
%    name: name string
%       M: Forward model
%              M.dipfit - leadfield specification
%       Y: data   [1x1 struct]
%       U: design [1x1 struct]
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
%   options.Spatial_type - EEG-MEG-LFP swtich (see spm_erp_L)
%__________________________________________________________________________
% %W% Karl Friston %E%


% check options 
%==========================================================================

% Filename and options
%--------------------------------------------------------------------------
try,   DCM.name;                 catch, DCM.name = 'DCM_ERP';  end

try,   D  = DCM.options.D;       catch, D = 1;                 end
try,   h  = DCM.options.h;       catch, h = 1;                 end
try,   nm = DCM.options.Nmodes;  catch, nm = 4;                end
try,   T1 = DCM.options.Tdcm(1); catch, T1 = DCM.Y.Time(1);    end
try,   T2 = DCM.options.Tdcm(2); catch, T2 = DCM.Y.Time(end);  end

% data type
%--------------------------------------------------------------------------
try,   DCM.M.Spatial_type = DCM.options.Spatial_type;
catch, DCM.M.Spatial_type = 1;                                 end

% Data
%==========================================================================
nt      = length(DCM.Y.xy);                 % number of trials
nr      = length(DCM.A{1});                 % number of sources
nc      = size(DCM.Y.xy{1},2);              % number of channels
nu      = size(DCM.U.X,2);                  % number of inputs
nx      = nr*9 + 1;                         % number of states
nm      = DCM.options.Nmodes;               % number of modes

% time window and bins for modelling
%--------------------------------------------------------------------------
[m, T1] = min(abs(DCM.Y.Time - T1));
[m, T2] = min(abs(DCM.Y.Time - T2));
j       = [T1:D:T2]';                       % time bins
xY.Time = DCM.Y.Time(j);                    % Time [ms] of downsampled data
for i = 1:nt
    xY.xy{i} = DCM.Y.xy{i}(j,:);
end
xY.y  = spm_cat(xY.xy(:));                  % concatenated response
xY.dt = DCM.Y.dt*D/1000;                    % sampling in seconds
ns    = length(j);                          % number of time samples

% confounds - DCT and Gamma functions
%--------------------------------------------------------------------------
if h == 0
    X0 = sparse(ns,1);
else
    X0 = spm_dctmtx(ns,h);
end
warning off
X0     = kron(speye(nt,nt),X0);
R0     = speye(ns*nt) - X0*inv(X0'*X0)*X0';  % null space of confounds
xY.X0  = X0;
warning on

% assume noise variance is the same over modes
%--------------------------------------------------------------------------
xY.Q   = {speye(ns*nm*nt,ns*nm*nt)};

% Inputs
%==========================================================================

% trial effects
%--------------------------------------------------------------------------
try
    xU.u = kron(DCM.U.X,ones(ns,1));
catch
    xU.u = sparse(nt*ns,0);
end

% stimulus parameters
%--------------------------------------------------------------------------
xU.dt   = xY.dt;
xU.dur  = xU.dt*(ns - 1);
xU.name = DCM.U.name;


% model specification and nonlinear system identification
%==========================================================================

% make model specification a global variable
%--------------------------------------------------------------------------
global M; M = DCM.M;

% prior moments
%--------------------------------------------------------------------------
[pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,DCM.M.dipfit.L,xU.dur);

% likelihood model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_erp';
M.g   = 'spm_gx_erp';
M.IS  = 'spm_int_U';
M.FS  = 'y*E';
M.x   = sparse(nx,1);
M.pE  = pE;
M.pC  = pC;
M.m   = nu;
M.n   = nx;
M.l   = nc;

% Project data onto the principal components of channel space
%--------------------------------------------------------------------------
y     = R0*xY.y;
S     = spm_svd(y'*y);
M.E   = S;
M.E   = M.E(:,[1:nm]);

% EM: This is a hidden loop that will perfom cnaonical feature selection
%--------------------------------------------------------------------------
ML    = exp(32);                          % switch off feature selection
for i = 1:8
    
    % inversion
    %----------------------------------------------------------------------
    [Qp,Cp,Ce,F] = spm_nlsi_GN(M,xU,xY);
    r            = y - R0*feval(M.IS,Qp,M,xU);
    
    if F < ML(end), break, end

    % Canonical feature selection
    %----------------------------------------------------------------------
    M.E   = spm_svd(inv(S'*r'*r*S)*(S'*y'*y*S),0);
    M.E   = orth(full(S*M.E(:,1:nm)));
    
    % reset stating estimates
    %----------------------------------------------------------------------
    ML(i) = F;
    M.P   = Qp;
    
end

% Bayesian inference {threshold = 0}
%--------------------------------------------------------------------------
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on

% neuronal responses (x)
%--------------------------------------------------------------------------
M.Spatial_type = 4;
x              = feval(M.IS,Qp,M,xU);
M.Spatial_type = DCM.M.Spatial_type;

% trial specific respsonses (in mode, channel and source space)
%--------------------------------------------------------------------------
for  i = 1:nt
    j     = [1:ns] + (i - 1)*ns;
    H{i}  = x(j,:)*M.L'*M.E;
    E{i}  = r(j,:)*M.E;
    Hc{i} = x(j,:)*M.L';
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

% and save
%--------------------------------------------------------------------------
save(DCM.name,'DCM');

