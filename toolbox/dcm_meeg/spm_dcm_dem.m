function DCM = spm_dcm_dem(DCM)
% Estimate parameters of a DCM-DEM model
% FORMAT DCM = spm_dcm_dem(DCM)
%
% DCM
%    name: name string
%       Lpos:  Source locations
%       xY:    data   [1x1 struct]
%       xU:    design [1x1 struct]
%
%   Sname: cell of source name strings
%
%   options.trials       - indices of trials
%   options.Lpos         - source location priors
%   options.Tdcm         - [start end] time window in ms
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.h            - number of DCT drift terms (usually 1 or 2)
%   options.Nmodes       - number of spatial models to invert
%   options.spatial      - 'ERP', 'LFP' or 'IMG'
%   options.onset        - stimulus onset (ms)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_dem.m 4281 2011-03-31 19:49:57Z karl $

% check options
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);

% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                   catch, DCM.name  = name;      end
try, h     = DCM.options.h;      catch, h         = 1;         end
try, Nm    = DCM.options.Nmodes; catch, Nm        = 8;         end
try, onset = DCM.options.onset;  catch, onset     = 60;        end
try, model = DCM.options.model;  catch, model     = 'DEM';     end


% Data and spatial model (use h only for de-trending data)
%==========================================================================
DCM    = spm_dcm_erp_data(DCM,h);
DCM    = spm_dcm_erp_dipfit(DCM);
xY     = DCM.xY;
xU     = DCM.xU;
M      = DCM.M;

% dimensions
%--------------------------------------------------------------------------
Nt     = length(xY.xy);                 % number of trials
Nr     = size(DCM.Lpos,2);              % number of sources
Ns     = size(xY.xy{1},1);              % number of time bins
Nc     = size(xY.xy{1},2);              % number of channels
Nx     = size(xU.X,2);                  % number of trial-specific effects

% check the number of modes is greater or equal to the number of sources
%--------------------------------------------------------------------------
Nm     = max(Nm,Nr);

% confounds - DCT: (force a parameter per channel = activity under x = 0)
%--------------------------------------------------------------------------
if h == 0
    X0 = zeros(Ns, h);
else
    X0 = spm_dctmtx(Ns, h);
end
T0     = speye(Ns) - X0*inv(X0'*X0)*X0';
xY.X0  = X0;

% Serial correlations (precision components) AR(1/4) model
%--------------------------------------------------------------------------
xY.Q   = {spm_Q(1/4,Ns,1)};


%-Inputs
%==========================================================================

% within-trial effects: adjust onset relative to PST
%--------------------------------------------------------------------------
M.ons  = onset - xY.pst(1);
xU.dt  = xY.dt;


%-Model specification and nonlinear system identification
%==========================================================================

% prior moments on parameters
%--------------------------------------------------------------------------
pE.R   = [0 3];
pE.f   = sparse(2,2);
pE.g   = 0;
pC     = speye(length(spm_vec(pE)))/32;

% priors on spatial model
%--------------------------------------------------------------------------
M.dipfit.model = model;
[gE,gC] = spm_L_priors(M.dipfit);

% augment with J (mapping) from hidden states to dipoles
%--------------------------------------------------------------------------
j = find(any(DCM.A{1}));
J = DCM.A{1}(:,1:max(j));

% allow the contribution of second (or further) states to be optimised
%--------------------------------------------------------------------------
cJ    = J;
eJ    = J*0;
for i = 1:Nr
    j       = find(J(i,:),1);
    cJ(i,j) = 0;
    eJ(i,j) = 1;
end
gE.J  = eJ;
gC    = spm_cat(spm_diag({gC, diag(spm_vec(cJ))}));
n     = size(J,2);

% likelihood model
%--------------------------------------------------------------------------
M.IS  = 'DEM_demo_MMN_gen';
M.FS  = 'spm_fy_erp';
M.G   = 'spm_lx_dem';
M.x   = sparse(n,1);
M.pE  = pE;
M.pC  = pC;
M.gE  = gE;
M.gC  = gC;
M.ns  = Ns;

%-Feature selection using principal components (U) of lead-field
%==========================================================================

% Spatial modes
%--------------------------------------------------------------------------
if Nc < Nm
    U     = speye(Nc);
    M.E   = U;
else
    dGdg  = spm_diff(M.G,gE,M,1);
    L     = spm_cat(dGdg);
    U     = spm_svd(L*L',exp(-8));
    try
        U = U(:,1:Nm);
    end
    M.E   = U;
end
Nm    = size(U,2);


% EM: inversion
%==========================================================================
[Qp,Qg,Cp,Cg,Ce,F] = spm_nlsi_N(M,xU,xY);

% Data ID
%==========================================================================
if isfield(M,'FS')
    try
        ID  = spm_data_id(feval(M.FS,xY.y,M));
    catch
        ID  = spm_data_id(feval(M.FS,xY.y));
    end
else
    ID  = spm_data_id(xY.y);
end


% Bayesian inference
%--------------------------------------------------------------------------
sw  = warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning(sw);


% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
L       = feval(M.G, Qg,M);                 % get gain matrix
[x DEM] = feval(M.IS,Qp,M,xU);              % prediction (source space)


% trial-specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
for i = 1:Nt
    y{i} = x{i}*L'*M.E;                 % prediction (sensor space)
    r{i} = T0*(xY.y{i}*M.E - y{i});     % residuals  (sensor space)
end


% store estimates in DCM
%--------------------------------------------------------------------------
DCM.M   = M;                    % model specification
DCM.xY  = xY;                   % data structure
DCM.xU  = xU;                   % input structure
DCM.Ep  = Qp;                   % conditional expectation f(x,u,p)
DCM.Cp  = Cp;                   % conditional covariances G(g)
DCM.Eg  = Qg;                   % conditional expectation
DCM.Cg  = Cg;                   % conditional covariances
DCM.Ce  = Ce;                   % conditional error covariance
DCM.Pp  = Pp;                   % conditional probability
DCM.H   = y;                    % conditional responses (y), projected space
DCM.K   = x;                    % conditional responses (x)
DCM.R   = r;                    % conditional residuals (y)
DCM.F   = F;                    % Laplace log evidence
DCM.ID  = ID;                   % data ID
DCM.DEM = DEM;                  % array of DEM structures

DCM.options.h      = h;
DCM.options.Nmodes = Nm;
DCM.options.onset  = onset;
DCM.options.model  = model;

% store estimates in D
%--------------------------------------------------------------------------
if strcmp(M.dipfit.type,'IMG')

    % Assess accuracy; signal to noise (over sources), SSE and log-evidence
    %----------------------------------------------------------------------
    for i = 1:Nt
        SSR(i) = sum(var(r{i}));
        SST(i) = sum(var(y{i} + r{i}));
    end
    R2    = 100*(sum(SST - SSR))/sum(SST);


    % reconstruct sources in dipole space
    %----------------------------------------------------------------------
    Nd    = M.dipfit.Nd;
    G     = sparse(Nd,0);

    % one dipole per source (i)
    %----------------------------------------------------------------------
    for i = 1:Nr
        G(M.dipfit.Ip{i},end + 1) = M.dipfit.U{i}*Qg.L(:,i);
    end

    clear J
    Is    = find(any(G,2));
    G     = G(Is,:);
    for i = 1:Nt
        J{i} = G*Qg.J*x{i}';
    end

    % get D and dipole space lead field
    %----------------------------------------------------------------------
    try, val = DCM.val;  catch, val = 1; end
    D     = spm_eeg_load(DCM.xY.Dfile);
    L     = spm_eeg_lgainmat(D,Is);
    L     = U'*L;

    % reduced data (for each trial)
    %----------------------------------------------------------------------
    for i = 1:Nt
        Y{i} = U'*xY.y{i}'*T0;
    end

    % fill in fields of inverse structure
    %----------------------------------------------------------------------
    inverse.trials = DCM.options.trials;   % trial or condition
    inverse.type   = 'DCM';                % inverse model
    inverse.J      = J;                    % Conditional expectation
    inverse.L      = L;                    % Lead field (reduced)
    inverse.R      = speye(Nc,Nc);         % Re-referencing matrix
    inverse.T      = T0;                   % temporal subspace
    inverse.U      = U;                    % spatial subspace
    inverse.Is     = Is;                   % Indices of active dipoles
    inverse.It     = DCM.xY.It;            % Indices of time bins
    inverse.Ic     = DCM.xY.Ic;            % Indices of good channels
    inverse.Y      = Y;                    % reduced data
    inverse.Nd     = Nd;                   % number of dipoles
    inverse.Nt     = Nt;                   % number of trials
    inverse.pst    = xY.pst;               % peri-stimulus time
    inverse.F      = DCM.F;                % log-evidence
    inverse.R2     = R2;                   % variance accounted for (%)
    inverse.dipfit = M.dipfit;             % forward model for DCM

    % append DCM results and save in structure
    %----------------------------------------------------------------------
    D.inv{end + 1}      = D.inv{val};
    D.inv{end}.date     = date;
    [pathstr,fname]     = fileparts(DCM.name);
    D.inv{end}.comment  = {fname};
    D.inv{end}.DCMfile  = DCM.name;
    D.inv{end}.inverse  = inverse;
    D.inv{end}.method   = 'Imaging';
    D.val               = length(D.inv);
    try
        D.inv{end}      = rmfield(D.inv{end},'contrast');
    end
    save(D);
end

% and save
%--------------------------------------------------------------------------
if spm_check_version('matlab','7') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end
assignin('base','DCM',DCM)
return
