function DCM = spm_dcm_erp(DCM)   
% Estimate parameters of a DCM model
% FORMAT DCM = spm_dcm_erp(DCM)   
%
% DCM     
%    name: name string
%       Y: data   [1x1 struct]
%       U: design [1x1 struct]
%       L: [nc x nr double] lead feild
%   Sname: cell of source name strings
%       A: {[nr x nr double]  [nr x nr double]  [nr x nr double]}
%       B: {[nr x nr double], ...}   Connection constaints
%       C: [nr x 1 double]
%     swd: working directory
%___________________________________________________________________
% %W% Karl Friston %E%

% swd
%-------------------------------------------------------------------
try
    swd     = DCM.swd;
    cd(swd)
catch
    swd     = pwd;
    DCM.swd = swd;
end

% Data
%===================================================================
nt    = length(DCM.Y.xy);                   % number of trials
nr    = size(DCM.L,2);                      % number of sources
nc    = size(DCM.L,1);                      % number of channels
ns    = size(DCM.Y.xy{1},1);                % number of samples
nu    = size(DCM.U.X,2);                    % number of inputs
nx    = nr*9 + 1;                           % number of states

% decimate (to about 8ms) and window
%-------------------------------------------------------------------
R     = fix(8/DCM.Y.dt);                    % decimation factor
j     = [R:R:ns]';                          % time bins
for i = 1:nt
    xy{i} = DCM.Y.xy{i}(j,:);
end
xY.y  = spm_cat(xy(:));                     % concatenated response
xY.dt = DCM.Y.dt*R/1000;                    % sampling in seconds
ns    = length(j);                          % number of time samples

% confounds - DCT and Gamma functions
%-------------------------------------------------------------------
try
    h = DCM.Y.h;
catch
    h = 3;
end
X0    = spm_dctmtx(ns,h);
X0    = kron(speye(nt,nt),X0);
R0    = speye(ns*nt) - X0*inv(X0'*X0)*X0';  % null space of confounds
xY.X0 = X0;

% orthogonalise response y (and reduce dimension if necessary)
%-------------------------------------------------------------------
nc    = min(nc,8);
xY.S  = spm_svd(xY.y'*R0*xY.y);
xY.S  = xY.S(:,[1:nc]);
xY.y  = xY.y*xY.S;
L     = xY.S'*DCM.L;

% Inputs
%===================================================================

% trial effects
%-------------------------------------------------------------------
try
    xU.u = kron(DCM.U.X,ones(ns,1));
catch
    xU.u = sparse(nt*ns,0);
end

% stimulus parameters
%-------------------------------------------------------------------
xU.dt   = xY.dt;
xU.dur  = xU.dt*(ns - 1);
xU.name = DCM.U.name;


% prior moments
%-------------------------------------------------------------------
[pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,L,xU.dur);

% model specification and nonlinear system identification
%-------------------------------------------------------------------
M.f   = 'spm_fx_erp';
M.g   = 'spm_gx_erp';
M.x   = sparse(nx,1);
M.pE  = pE;
M.pC  = pC;
M.m   = nu;
M.n   = nx;
M.l   = nc;
M.IS  = 'spm_int_U';


% EM estimation
%-------------------------------------------------------------------
[Ep,Cp,Ce,F] = spm_nlsi_GN(M,xU,xY);

% Bayesian inference {threshold = 0}
%-------------------------------------------------------------------
warning off
Pp    = 1 - spm_Ncdf(0,abs(Ep - pE),diag(Cp));
Pp    = spm_erp_pack(Pp,M.m,nr);
Qp    = spm_erp_pack(Ep,M.m,nr);
warning on

% predicted responses (y) and residuals (r) (in channel space)
%-------------------------------------------------------------------
y     = feval(M.IS,Ep,M,xU);
r     = R0*(xY.y - y);
y     = y*xY.S';
r     = r*xY.S';

% neuronal responses (x)
%-------------------------------------------------------------------
Qp.L  = [];
x     = feval(M.IS,spm_vec(Qp),M,xU);

% trial specific respsonses
%-------------------------------------------------------------------
for  i = 1:nt
    j    = [1:ns] + (i - 1)*ns;
    H{i} = y(j,:);
    E{i} = r(j,:);
    K{i} = x(j,:);
end


% store estimates in DCM
%-------------------------------------------------------------------
DCM.M    = M;                             % model specification
DCM.xY   = xY;                            % data structure
DCM.xU   = xU;                            % input structure
DCM.Ep   = Ep;                            % conditional expectation
DCM.Cp   = Cp;                            % conditional covariances
DCM.Qp   = Qp;                            % conditional expectation
DCM.Pp   = Pp;                            % conditional probability
DCM.H    = H;                             % conditional responses (y)
DCM.K    = K;                             % conditional responses (x)
DCM.R    = E;                             % conditional residuals (y)
DCM.Ce   = Ce;                            % ReML error covariance
DCM.F    = F;                             % Laplace log evidence

% and save
%-------------------------------------------------------------------
save(['DCM', DCM.name],'DCM');

