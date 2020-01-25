function [dfdx,dfdu,dgdx] = spm_dcm2ssm(P,M)
% linearises a dynamic causal model about an expansion point
% FORMAT [dfdx,dfdu,dgdx] = spm_dcm2ssm(P,M)
%
% P - model parameters
% M - model (with flow M.f and expansion point M.x and M.u)
%    M.f     - dx/dt = f(x,u,P,M)  {function string or m-file}
%    M.g     - y     = g(x,u,P,M)  {function string or m-file}
%    M.x [default: sparse(M.n,1)]
%    M.u [default: sparse(M.m,1)]
% 
% dfdx - Jacobian
% dfdu - input  matrix
% dgdx - output matrix
%
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm2ssm.m 7774 2020-01-25 18:07:03Z karl $


% get local linear approximation
%==========================================================================

% check expansion points
%--------------------------------------------------------------------------
try, M.x; catch,      M.x = sparse(M.n,1); end
try, M.u; catch, try, M.u = sparse(M.m,1); catch, M.u = M.x; end, end

% delay operator - if not specified already
%--------------------------------------------------------------------------
if isfield(M,'D')
    
    dfdx = spm_diff(M.f,M.x,M.u,P,M,1);
    dfdu = spm_diff(M.f,M.x,M.u,P,M,2);
    D    = M.D;
    
else
    
    if nargout(M.f) == 4
        [f,dfdx,D,dfdu] = feval(M.f,M.x,M.u,P,M);
        
    elseif nargout(M.f) == 3
        [f,dfdx,D]      = feval(M.f,M.x,M.u,P,M);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        
    elseif nargout(M.f) == 2
        [f,dfdx]        = feval(M.f,M.x,M.u,P,M);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        D               = 1;
    else
        dfdx            = spm_diff(M.f,M.x,M.u,P,M,1);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        D               = 1;
    end
end

% state-space matrix operators
%==========================================================================
if isfield(M,'g')
    if nargout(M.g) == 2
        [g,dgdx] = feval(M.g,M.x,M.u,P,M);
    else
        dgdx     = spm_diff(M.g,M.x,M.u,P,M,1);
    end
else
    dgdx     = speye(numel(M.u),numel(M.x));
end
dfdx  = D*dfdx;
dfdu  = D*dfdu;
