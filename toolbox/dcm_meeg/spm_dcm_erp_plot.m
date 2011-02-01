function x = spm_dcm_erp_plot(DCM)
% plots predicted source activity
% FORMAT x = spm_dcm_erp_plot(DCM)
%
% DCM - DCM structure:
% store estimates in DCM
%--------------------------------------------------------------------------
% DCM.M  - model specification
% DCM.xY - data structure
% DCM.xU - input structure
% DCM.Ep - conditional expectation f(x,u,p)
% DCM.Cp - conditional covariances G(g)
% DCM.Eg - conditional expectation
% DCM.Cg - conditional covariances
% DCM.Pp - conditional probability
% DCM.H  - conditional responses (y), projected space
% DCM.K  - conditional responses (x)
% DCM.R  - conditional residuals (y)
% DCM.F  - Laplace log evidence
% DCM.L  - Laplace log evidence components
% DCM.ID - data ID
% 
% 
% DCM.options.h
% DCM.options.Nmodes
% DCM.options.onset
% DCM.options.model
% DCM.options.lock
% DCM.options.symm
%
% x{i}   - source activity contributing sources {trial i}
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_erp_plot.m 4186 2011-02-01 20:11:32Z karl $

% check options
%==========================================================================

% unpack
%--------------------------------------------------------------------------
xY     = DCM.xY;
xU     = DCM.xU;
M      = DCM.M;

% dimensions
%--------------------------------------------------------------------------
Nr     = size(DCM.C,1);                 % number of sources
Nu     = size(DCM.C,2);                 % number of exogenous inputs
Ns     = size(xY.xy{1},1);              % number of time bins

% confounds and parameters
%--------------------------------------------------------------------------
T0     = speye(Ns) - xY.X0*((xY.X0'*xY.X0)\xY.X0');
Qp     = DCM.Ep;

% neuronal responses
%--------------------------------------------------------------------------
x     = feval(M.IS,Qp,M,xU);

% trial-specific responses (in source space)
%--------------------------------------------------------------------------
x0    = ones(Ns,1)*spm_vec(M.x)';       % expansion point for states
j     = kron(DCM.Eg.J,ones(1,Nr));      % Indices of contributing states
j     = logical(j);
for i = 1:numel(x)
    x{i} = x{i} - x0;                   % centre on expansion point
    x{i} = x{i}(:,j);                   % Depolarization in sources
end


% store estimates in DCM
%--------------------------------------------------------------------------
DCM.K  = x;                             % conditional sensitivity

% display
%--------------------------------------------------------------------------
spm_dcm_erp_results(DCM,'ERPs (sources)');
