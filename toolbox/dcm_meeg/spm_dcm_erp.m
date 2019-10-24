function [DCM,dipfit] = spm_dcm_erp(DCM)
% Estimate parameters of a DCM model (Variational Lapalce)
% FORMAT [DCM,dipfit] = spm_dcm_erp(DCM)
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
%   options.trials       - indices of trials
%   options.Tdcm         - [start end] time window in ms
%   options.D            - time bin decimation       (usually 1 or 2)
%   options.h            - number of DCT drift terms (usually 1 or 2)
%   options.Nmodes       - number of spatial models to invert
%   options.analysis     - 'ERP', 'SSR' or 'IND'
%   options.model        - 'ERP', 'SEP', 'CMC', 'CMM', 'NMM' or 'MFM'
%   options.spatial      - 'ECD', 'LFP' or 'IMG'
%   options.onset        - stimulus onset (ms)
%   options.dur          - and dispersion (sd)
%   options.CVA          - use CVA for spatial modes [default = 0]
%   options.Nmax         - maxiumum number of iterations [default = 64]
%
% dipfit - Dipole structure (for electromagnetic forward model)
%        See spm_dcm_erp_dipfit:  this field is removed from DCM.M to save
%        memory - and is offered as an output argument if needed
%
% The scheme can be initialised with parameters for the neuronal model
% and spatial (observer) model by specifying the fields DCM.P and DCM.Q, 
% respectively. If previous priors (DCM.M.pE and pC or DCM.M.gE and gC or 
% DCM.M.hE and hC) are specified, they will be used. Explicit priors can be
% useful for Bayesian parameter averaging - but would not normally be
% called upon - because prior constraints are specified by DCM.A, DCM.B,...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp.m 7679 2019-10-24 15:54:07Z spm $

% check options (and clear persistent variables)
%==========================================================================
drawnow
clear spm_erp_L
name = sprintf('DCM_%s',date);


% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                        catch, DCM.name  = name;        end
try, DCM.xU;                          catch, DCM.xU.X  = sparse(1,0); end
try, Nm       = DCM.options.Nmodes;   catch, Nm        = 8;           end
try, onset    = DCM.options.onset;    catch, onset     = 60;          end
try, dur      = DCM.options.dur;      catch, dur       = 16;          end
try, model    = DCM.options.model;    catch, model     = 'CMC';       end
try, lock     = DCM.options.lock;     catch, lock      = 0;           end
try, multC    = DCM.options.multiC;   catch, multC     = 0;           end
try, symm     = DCM.options.symmetry; catch, symm      = 0;           end
try, CVA      = DCM.options.CVA;      catch, CVA       = 0;           end
try, Nmax     = DCM.options.Nmax;     catch, Nmax      = 64;          end
try, DATA     = DCM.options.DATA;     catch, DATA      = 1;           end
try, Nmax     = DCM.M.Nmax;           catch, Nmax      = Nmax;        end


% symmetry contraints for ECD models only
%--------------------------------------------------------------------------
if ~strcmp(DCM.options.spatial,'ECD'), symm = 0; end

% disallow IMG solutions for generic DCMs
%--------------------------------------------------------------------------
if isstruct(model) && strcmp(DCM.options.spatial,'IMG')
    DCM.options.spatial = 'ECD';
end


% Data and spatial model
%==========================================================================
if DATA
    DCM = spm_dcm_erp_data(DCM);
    DCM = spm_dcm_erp_dipfit(DCM,1);
end
xY      = DCM.xY;
xU      = DCM.xU;
M       = DCM.M;

if ~isfield(xY,'X0'),    xY.X0 = sparse(size(xY.y{1},1),0); end
if ~isfield(xU,'X'),     xU.X  = sparse(1,0);               end
if ~isfield(xY,'scale'), xY.scale  = 1;                     end

% dimensions
%--------------------------------------------------------------------------
Nt      = length(xY.y);                  % number of trials
Nr      = size(DCM.C,1);                 % number of sources
Nu      = size(DCM.C,2);                 % number of exogenous inputs
Ns      = size(xY.y{1},1);               % number of time bins
Nc      = size(xY.y{1},2);               % number of channels
Nx      = size(xU.X,2);                  % number of trial-specific effects

% check the number of modes is greater or equal to the number of sources
%--------------------------------------------------------------------------
Nm      = max(Nm,Nr);

% confounds - residual forming matrix
%--------------------------------------------------------------------------
if isfield(xY,'R')
    M.R = xY.R;
else
    M.R = speye(Ns) - xY.X0*((xY.X0'*xY.X0)\xY.X0');
end


% Serial correlations (precision components) AR model
%--------------------------------------------------------------------------
xY.Q   = {spm_Q(1/2,Ns,1)};


%-Inputs
%==========================================================================

% between-trial effects
%--------------------------------------------------------------------------
try
    if length(DCM.B) < Nx
        for i = 1:Nx
            DCM.B{i} = sparse(Nr,Nr);
        end
    end
catch
    xU.X  = sparse(1,0);
    DCM.B = {};
end

% within-trial effects: adjust onset relative to PST
%--------------------------------------------------------------------------
M.ons  = onset - xY.pst(1);
M.dur  = dur;
xU.dt  = xY.dt;


%-Model specification and nonlinear system identification
%==========================================================================
try, M = rmfield(M,'g'); end

% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);

% check for trial specific inputs
%--------------------------------------------------------------------------
if multC
    pE.C(:,:,2) = pE.C(:,:,1);
    pC.C(:,:,2) = pC.C(:,:,1);
end

% priors on spatial model
%--------------------------------------------------------------------------
M.dipfit.model = model;
[gE,gC]        = spm_L_priors(M.dipfit);

% Set prior correlations (locking trial effects and dipole orientations
%--------------------------------------------------------------------------
if lock, pC = spm_dcm_lock(pC);      end
if symm, gC = spm_dcm_symm(gC,gE);   end


% hyperpriors (assuming a high signal to noise)
%--------------------------------------------------------------------------
hE      = 6;
hC      = 1/128;

% check for previous priors
%--------------------------------------------------------------------------
try
    pE  = M.pE;
    pC  = M.pC;
    fprintf('Using specified priors (for neural model)\n')
end
try
    gE  = M.gE;
    gC  = M.gC;
    fprintf('Using specified priors (for spatial model)\n')
end
try
    hE  = M.hE;
    hC  = M.hC;
    fprintf('Using specified priors (for noise precision)\n')
end


%-Feature selection using (canonical) eigenmodes of lead-field
%==========================================================================
if CVA
    M.U  = spm_dcm_eeg_channelmodes(M.dipfit,Nm,xY);
else
    M.U  = spm_dcm_eeg_channelmodes(M.dipfit,Nm);
end

% scale data features
%--------------------------------------------------------------------------
scale    = std(spm_vec(spm_fy_erp(xY.y,M)));
xY.y     = spm_unvec(spm_vec(xY.y)/scale,xY.y);
xY.scale = xY.scale/scale;


% likelihood model
%==========================================================================

% Use TFM integration scheme (with plasticity) if indicated
%--------------------------------------------------------------------------
if isfield(M,'TFM')
     IS = 'spm_csd_int';
else
     IS = 'spm_gen_erp';
end

% intial states and equations of motion
%--------------------------------------------------------------------------
[x,f,h] = spm_dcm_x_neural(pE,model);

M.FS   = 'spm_fy_erp';
M.G    = 'spm_lx_erp';
M.IS   = IS;
M.f    = f;
M.h    = h;
M.x    = x;
M.pE   = pE;
M.pC   = pC;
M.gE   = gE;
M.gC   = gC;
M.hE   = hE;
M.hC   = hC;
M.m    = Nu;
M.n    = length(spm_vec(M.x));
M.l    = Nc;
M.ns   = Ns;
M.Nmax = Nmax;

% re-intialise states
%--------------------------------------------------------------------------
M.x    = spm_dcm_neural_x(pE,M);
    

% EM: inversion
%==========================================================================
[Qp,Qg,Cp,Cg,Ce,F,LE] = spm_nlsi_N(M,xU,xY);


% Data ID
%==========================================================================
if isfield(M,'FS')
    try
        ID = spm_data_id(feval(M.FS,xY.y,M));
    catch
        ID = spm_data_id(feval(M.FS,xY.y));
    end
else
    ID = spm_data_id(xY.y);
end


% Bayesian inference
%--------------------------------------------------------------------------
sw  = warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning(sw);


% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
L   = feval(M.G, Qg,M);                 % get gain matrix
x   = feval(M.IS,Qp,M,xU);              % prediction (source space)


% trial-specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
try
    j = find(kron(Qg.J,ones(1,Nr)));    % Indices of contributing states
catch
    j = find(spm_cat(Qg.J));
end
x0  = ones(Ns,1)*spm_vec(M.x)';         % expansion point for states
for i = 1:Nt
    K{i} = x{i} - x0;                   % centre on expansion point
    y{i} = M.R*K{i}*L'*M.U;             % prediction (sensor space)
    r{i} = M.R*xY.y{i}*M.U - y{i};      % residuals  (sensor space)
    K{i} = K{i}(:,j);                   % Depolarization in sources
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
DCM.Ce = Ce;                   % conditional error
DCM.Pp = Pp;                   % conditional probability
DCM.H  = y;                    % conditional responses (y), projected space
DCM.K  = K;                    % conditional responses (x) (contributing)
DCM.x  = x;                    % conditional responses (x) (all states)
DCM.R  = r;                    % conditional residuals (y)
DCM.F  = F;                    % Laplace log evidence
DCM.L  = LE;                   % Laplace log evidence components
DCM.ID = ID;                   % data ID

DCM.options.Nmodes   = size(M.U,2);
DCM.options.onset    = onset;
DCM.options.dur      = dur;
DCM.options.model    = model;
DCM.options.lock     = lock;
DCM.options.symm     = symm;
DCM.options.analysis = 'ERP';


% store estimates in D
%--------------------------------------------------------------------------
if strcmp(M.dipfit.type,'IMG') && DATA

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

    % one dipole per subpopulation (p)
    %----------------------------------------------------------------------
    if iscell(Qg.L)
        for p = 1:length(Qg.L)
            for i = 1:Nr
                G(M.dipfit.Ip{i},end + 1) = M.dipfit.U{i}*Qg.L{p}(:,i);
            end
        end

        % one dipole per source (i)
        %----------------------------------------------------------------------
    else
        for i = 1:Nr
            G(M.dipfit.Ip{i},end + 1) = M.dipfit.U{i}*Qg.L(:,i);
        end
        G = kron(Qg.J,G);
    end
    Is    = find(any(G,2));
    Ix    = find(any(G,1));
    G     = G(Is,Ix);
    for i = 1:Nt
        J{i} = G*K{i}';
    end

    % get D and dipole space lead field
    %----------------------------------------------------------------------
    try, val = DCM.val; catch, val = 1; end
    D     = spm_eeg_load(DCM.xY.Dfile);
    L     = spm_eeg_lgainmat(D, Is, DCM.xY.name);
    L     = M.U'*L;

    % reduced data (for each trial
    %----------------------------------------------------------------------
    for i = 1:Nt
        Y{i} = M.U'*xY.y{i}'*M.R;
    end

    % fill in fields of inverse structure
    %----------------------------------------------------------------------
    inverse.trials   = DCM.options.trials;   % trial or condition
    inverse.modality = {DCM.xY.modality};    % modality
    inverse.type     = 'DCM';                % inverse model
    inverse.J        = J;                    % Conditional expectation
    inverse.L        = L;                    % Lead field (reduced)
    inverse.R        = speye(Nc,Nc);         % Re-referencing matrix
    inverse.T        = M.R;                  % temporal subspace
    inverse.U        = M.U;                  % spatial subspace
    inverse.Is       = Is;                   % Indices of active dipoles
    inverse.It       = DCM.xY.It;            % Indices of time bins
    inverse.Ic       = DCM.xY.Ic;            % Indices of good channels
    inverse.Y        = Y;                    % reduced data
    inverse.Nd       = Nd;                   % number of dipoles
    inverse.Nt       = Nt;                   % number of trials
    inverse.pst      = xY.pst;               % peri-stimulus time
    inverse.F        = DCM.F;                % log-evidence
    inverse.R2       = R2;                   % variance accounted for (%)
    inverse.dipfit   = M.dipfit;             % forward model for DCM

    % append DCM results and save in structure
    %----------------------------------------------------------------------
    D.inv{end + 1}      = D.inv{val};
    D.inv{end}.date     = date;
    [dummy,fname]       = fileparts(DCM.name);
    D.inv{end}.comment  = {fname};
    D.inv{end}.DCMfile  = DCM.name;
    D.inv{end}.inverse  = inverse;
    D.val               = length(D.inv);
    try
        D.inv{end}      = rmfield(D.inv{end},'contrast');
    end
    save(D);
end

% remove dipfit stucture to save memory
%--------------------------------------------------------------------------
dipfit = DCM.M.dipfit;
DCM.M  = rmfield(DCM.M,'dipfit');

% and save
%--------------------------------------------------------------------------
if DATA
    save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
end
