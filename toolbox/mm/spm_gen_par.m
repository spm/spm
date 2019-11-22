function Q = spm_gen_par(P,M,U)
% Generates condition specific parameters using DCM for M/EEG
% FORMAT Q = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.X  - between-trial effects (encodes the number of trials)
%   U.dt - time bins for within-trial effects
%
% Q   - Condition specific parameters 
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging
% Amirhossein Jafarian
% $Id $
%--------------------------------------------------------------------------
if nargin < 3, U.X = sparse(1,0); end

% check input u = f(t,P,M) and switch off full delay operator
%--------------------------------------------------------------------------
try M.fu; catch, M.fu  = 'spm_erp_u'; end
try M.ns; catch, M.ns  = 128;         end
try M.N;  catch, M.N   = 0;           end
try U.dt; catch, U.dt  = 0.004;       end

% peristimulus time
%--------------------------------------------------------------------------
if nargout > 1
    pst = (1:M.ns)*U.dt - M.ons/1000;
end

% within-trial (exogenous) inputs
%==========================================================================
if ~isfield(U,'u')
    
    % peri-stimulus time inputs
    %----------------------------------------------------------------------
    U.u = feval(M.fu,(1:M.ns)*U.dt,P,M);
    
end

if isfield(M,'u')
    
    % remove M.u to preclude endogenous input
    %----------------------------------------------------------------------
    M = rmfield(M,'u');
    
end

% between-trial (experimental) inputs
%==========================================================================
if isfield(U,'X')
    X = U.X;
else
    X = sparse(1,0);
end

if ~size(X,1)
    X = sparse(1,0);
end

% trials
%==========================================================================
y      = cell(size(X,1),1);
c = P.xc ;
% condition-specific parameters
%----------------------------------------------------------------------
Q   = spm_gen_Q(P,X(c,:));

