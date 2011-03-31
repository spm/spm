function DCM = spm_dcm_phase(DCM)   
% Estimate parameters of a DCM model of phase-coupled responses
% FORMAT DCM = spm_dcm_phase(DCM)   
%
% DCM     
%    name: name string
%       M:  Forward model
%              M.dipfit - leadfield specification
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%
%   Sname: cell of source name strings
%
%   Connection constraints:
%
%       A: {[nr x nr double] }
%       B: {[nr x nr double]}   for GUI specification
%                               (Nfourier=1 & only sine terms)
%   or
%
%       As: [nr x nr x Nfourier]
%       Ac: [nr x nr x Nfourier]
%       Bs: [nr x nr x Nfourier]   
%       Bc: [nr x nr x Nfourier]   for script specification
%
%   options.type         - 'ECD' 
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging
%
% Will Penny
% $Id: spm_dcm_phase.m 4281 2011-03-31 19:49:57Z karl $


% check options 
%==========================================================================
clear spm_erp_L

% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                  catch, DCM.name           = 'DCM_PHASE'; end
try, DCM.options.Nmodes;        catch, DCM.options.Nmodes = 4;         end
try, h     = DCM.options.h;     catch, h                  = 1;         end
try, onset = DCM.options.onset; catch, onset              = 80;        end

% Data and spatial model
%==========================================================================
DCM    = spm_dcm_erp_dipfit(DCM, 1);

% Get data if not gotten already
%if ~isfield(DCM.xY,'source')   
    DCM  = spm_dcm_phase_data(DCM);
%end

disp('Estimating Model ...');

xY     = DCM.xY;
xU     = DCM.xU;
xU.dt  = xY.dt;

if isfield(DCM,'A')
    % If using UI, copy model structures into As and Bs matrices
    DCM.As=DCM.A{1};
    DCM.Bs{1}=DCM.B{1};
end

% dimensions
%--------------------------------------------------------------------------
Nt     = length(xY.y);                  % number of trials
Nr     = size(DCM.As,1);                % number of sources
nx     = Nr;                            % number of states

Ns     = size(xY.y{1},1);               % number of samples
xU.u   = zeros(Ns,1);

Nu     = size(xU.X,2);                  % number of trial-specific effects

% Set error covariance components - one per region
%--------------------------------------------------------------------------
xY.Q   = spm_Ce(Ns*Nt*ones(1,Nr));

% Offsets
xY.X0  = sparse(Ns,0);

% Inputs
%==========================================================================

% trial-specific effects
%--------------------------------------------------------------------------
try
    if length(DCM.Bs) ~= Nu;
        warndlg({'Please ensure number of trial specific effects', ...
                 'encoded by DCM.xU.X & DCM.B are the same'})
    end
catch
    DCM.Bs = {};
end

% model specification and nonlinear system identification
%==========================================================================
M      = DCM.M;

M.freq = mean(DCM.options.Fdcm);
M.fb = 0.5*(DCM.options.Fdcm(2)-DCM.options.Fdcm(1));

try, M = rmfield(M,'g');  end
try, M = rmfield(M,'FS'); end

if ~isfield(M,'dipfit')
    M.dipfit=[];
end

% prior moments
%--------------------------------------------------------------------------
[pE,gE,pC,gC] = spm_phase_priors(DCM,M.fb,M.dipfit);

% likelihood model
%--------------------------------------------------------------------------
M.IS  = 'spm_gen_phase';
M.f   = 'spm_fx_phase';
M.G   = 'spm_lx_phase';

M.pE  = pE;
M.pC  = pC;
M.gE  = gE;
M.gC  = gC;
M.n   = nx;
M.l   = Nr;
M.ns  = Ns;

% Don't plot progress
M.nograph=0;

% Initial state
try
    M.x;
catch
    M.x=zeros(Nr,1);
end

% Set up initial phases 
for n=1:length(DCM.xY.y)
    xx=DCM.xY.y{n}(1,:);
    M.trial{n}.x = double(xx(:)');
end


% keyboard
% myx=spm_gen_phase(pE,M,xU);

% EM: inversion
%--------------------------------------------------------------------------
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


% Bayesian inference {threshold = prior} NB Prior on A,B  and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on','SPM:negativeVariance');


% neuronal and sensor responses (x and y)
%--------------------------------------------------------------------------
x   = feval(M.IS,Qp,M,xU);        % prediction (source space)
L   = feval(M.G, Qg,M);           % get gain matrix


% trial-specific responses (in mode, channel and source space)
%--------------------------------------------------------------------------
for i = 1:Nt
    s    = x{i};                  % prediction (source space)
    y{i} = s*L';                  % prediction (sensor space)
    %r  = xY.y{i} - y;            % residuals  (sensor space)
    
end

% Change back design matrix to user specified
xU.X=DCM.xU.oldX;

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
DCM.Ce = Ce;                   % conditional error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
DCM.y  = y;                    % Model predictions

% and save
%--------------------------------------------------------------------------
if spm_check_version('matlab','7') >= 0
    save(DCM.name, '-V6', 'DCM');
else
    save(DCM.name, 'DCM');
end
assignin('base','DCM',DCM)
return
