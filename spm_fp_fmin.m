function [P] = spm_fp_fmin(M)
% optimises the parameters with respect to an equilibrium density
% FORMAT [P] = spm_fp_fmin(M)
%
% M   - model structure with desired density specified by M(1).fq and
%       support specified by M(1).X = spm_ndgrid(x)
%
% P   - optimised parameters
%
%--------------------------------------------------------------------------
% This routine uses EM (spm_nlsi_NG) and the Fokker Planck formulation to
% minimise the difference between the flow and dispersion terms induced by
% the free parameters of the flow (M(1),f).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fp_fmin.m 4136 2010-12-09 22:22:28Z guillaume $
 
 
% specify function returning the flow-dependent part of dp/dt 'spm_fp_fun'
%--------------------------------------------------------------------------
M      = M(1);
try, M = rmfield(M,'hE'); end
try, M = rmfield(M,'hC'); end
U      = [];
M.IS   = 'spm_fp_fun';
 
% Dispersion
%--------------------------------------------------------------------------
N     = size(M.X,1);
D     = inv(M.W)/2;
for i = 1:N
    Y.y(i,1) = trace(D*spm_cat(spm_diff(M.fq,M.X(i,:),[1 1])')');
end

% Optimise
%--------------------------------------------------------------------------
P     = spm_nlsi_GN(M,U,Y);


