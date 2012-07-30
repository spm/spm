function [y] = spm_gen_erp(P,M,U)
% Generates a prediction of trial-specific source activity
% FORMAT [y] = spm_gen_erp(P,M,U)
%
% P - parameters
% M - neural-mass model structure
% U - trial-effects
%   U.X  - between-trial effects (encodes the number of trials)
%   U.dt - time bins for within-trial effects
%
% y - {[ns,nx];...} - predictions for nx states {trials}
%                   - for ns samples
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_gen_erp.m 4814 2012-07-30 19:56:05Z karl $

% default inputs - one trial (no between-tria effects)
%--------------------------------------------------------------------------
if nargin < 3, U.X = sparse(1,0); end

% within-trial (exogenous) inputs
%==========================================================================
if ~isfield(U,'u')
    
    % check input u = f(t,P,M)
    %----------------------------------------------------------------------
    try, fu  = M.fu; catch,  fu  = 'spm_erp_u'; end
    try, ns  = M.ns; catch,  ns  = 128;         end
    try, dt  = U.dt; catch,  dt  = 0.004;       end
    
    % peri-stimulus time inputs
    %----------------------------------------------------------------------
    U.u = feval(fu,(1:ns)*U.dt,P,M);
    
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
    
    % trial-specific effects
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        % extrinsic (forward and backwards) connections
        %------------------------------------------------------------------
        for j = 1:length(Q.A)
            Q.A{j} = Q.A{j} + X(c,i)*P.B{i};
        end
        
        % intrinsic connections
        %------------------------------------------------------------------
        try
            Q.H(:,1) = Q.H(:,1) + X(c,i)*diag(P.B{i});
        catch
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

