function [y,pst] = spm_gen_erp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y,pst] = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.X  - between-trial effects (encodes the number of trials)
%   U.dt - time bins for within-trial effects
%
% y   - {[ns,nx];...} - predictions for nx states {trials}
%                     - for ns samples
% pst - peristimulus time (seconds)
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_erp.m 5407 2013-04-12 19:03:29Z karl $

% default inputs - one trial (no between-trial effects)
%--------------------------------------------------------------------------
if nargin < 3, U.X = sparse(1,0); end

% check input u = f(t,P,M)
%--------------------------------------------------------------------------
try, M.fu; catch, M.fu  = 'spm_erp_u'; end
try, M.ns; catch, M.ns  = 128;         end
try, U.dt; catch, U.dt  = 0.004;       end

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
try
    X = U.X;
    if ~size(X,1)
        X = sparse(1,0);
    end
catch
    X = sparse(1,0);
end


% cycle over trials
%--------------------------------------------------------------------------
y      = cell(size(X,1),1);
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
    
    % trial-specific effects on C (first effect only)
    %----------------------------------------------------------------------
    try
        Q.C = Q.C(:,:,1) + X(c,1)*Q.C(:,:,2);
    end
    
    % trial-specific effects on A (connections)
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        % extrinsic (driving) connections
        %------------------------------------------------------------------
        for j = 1:length(Q.A)
            Q.A{j} = Q.A{j} + X(c,i)*P.B{i};
        end
        
        % modulatory connections
        %------------------------------------------------------------------
        try
            Q.M  = Q.M + X(c,i)*P.N{i};
        end
        
        % intrinsic connections
        %------------------------------------------------------------------
        try
            Q.G(:,1) = Q.G(:,1) + X(c,i)*diag(P.B{i});
        end
        
    end
    
    % solve for steady-state - for each condition
    %----------------------------------------------------------------------
    M.x  = spm_dcm_neural_x(Q,M);
    
    % integrate DCM for this trial
    %----------------------------------------------------------------------
    y{c} = spm_int_L(Q,M,U);
    
end

