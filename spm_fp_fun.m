function [y] = spm_fp_fun(P,M,U)
% returns the predicted diffusion for Fokker Planck optimisation
% FORMAT [y] = spm_fp_fun(P,M,U)
%
% P = spm_vec(P)
% P.a - 0th order coefficients of force
% P.b - 1st order coefficients of force
% P.c - 2nd order coefficients of force
%
% M   - model specifying flow(M(1).f; density M(1).fq and support M(1).X
% U   - inputs
%
% y   - prediction
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp_fun.m 2030 2008-09-02 18:28:40Z karl $
 
% default: first level of hierarchical model
%--------------------------------------------------------------------------
M   = M(1);
 
% default expansion point for inputs
%--------------------------------------------------------------------------
try
    u = U.u;
catch
    u = sparse(M.m,1);
end
 
% predicted dispersion
%--------------------------------------------------------------------------
N     = length(M.X);
for i = 1:N
    f      = feval(M.f,M.X(i,:)',u,P);
    p      = feval(M.fq,M.X(i,:));
    dfdx   = spm_diff(M.f,M.X(i,:)',u,P,1);
    dpdx   = spm_diff(M.fq,M.X(i,:),1);
    y(i,1) = trace(p*dfdx) + dpdx*f;
end


    
