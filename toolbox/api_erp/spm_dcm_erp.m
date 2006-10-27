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
%   options.Spatial_type - EEG-MEG-LFP swtich (see spm_erp_L)
%__________________________________________________________________________
% %W% Karl Friston %E%


% check options 
%==========================================================================

% Filename and options
%--------------------------------------------------------------------------
try,   DCM.name;                 catch, DCM.name = 'DCM_ERP';  end
try,   h  = DCM.options.h;       catch, h  = 1;                end
try,   nm = DCM.options.Nmodes;  catch, nm = 4;                end

% data type
%--------------------------------------------------------------------------
try,   DCM.M.Spatial_type = DCM.options.Spatial_type;
catch, DCM.M.Spatial_type = 1;                                 end


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
nt     = length(xY.xy);                 % number of trials
nr     = length(DCM.A{1});                 % number of sources
nc     = size(xY.xy{1},2);              % number of channels
ns     = size(xY.xy{1},1);              % number of time bins
nu     = size(xU.X,2);                  % number of inputs
nx     = nr*9 + 1;                         % number of states


% confounds - DCT
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
xY.Q   = {speye(ns*nt*nm,ns*nt*nm)};

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
    xU.u  = kron(xU.X,ones(ns,1));
catch
    xU.u  = sparse(nt*ns,0);
end

% stimulus parameters
%--------------------------------------------------------------------------
xU.dt   = xY.dt;
xU.dur  = xU.dt*(ns - 1);

% model specification and nonlinear system identification
%==========================================================================

% make model specification a global variable
%--------------------------------------------------------------------------
global M; M = DCM.M;

% prior moments
%--------------------------------------------------------------------------
Lpos    = M.dipfit.L.pos;
[pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,Lpos,xU.dur);

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
M.l   = nc;
M.ns  = ns;

% Feature selection using principal components of channel space
%--------------------------------------------------------------------------
try
    M.S = M.S(:,nm);
catch
    y   = R0*xY.y;
    S   = spm_svd(y'*y);
    M.E = S(:,[1:nm]);
end

% EM: inversion
%--------------------------------------------------------------------------
[Qp,Cp,Ce,F] = spm_nlsi_GN(M,xU,xY);

% Bayesian inference {threshold = 0}
%==========================================================================
warning off
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning on

% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
M.Spatial_type = 4;
x              = feval(M.IS,Qp,M,xU);        % prediction (source space)
M.Spatial_type = DCM.M.Spatial_type;
y              = feval(M.IS,Qp,M,xU);        % prediction (sensor space)
r              = R0*(xY.y - y);              % prediction error


% trial specific respsonses (in mode, channel and source space)
%--------------------------------------------------------------------------
L     = spm_erp_L(Qp);
for i = 1:nt
    j     = [1:ns] + (i - 1)*ns;
    H{i}  = x(j,:)*L'*M.E;
    E{i}  = r(j,:)*M.E;
    Hc{i} = x(j,:)*L';
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


