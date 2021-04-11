function J = spm_ness_J(P,M,X)
% returns the Jacobian given a polynomial parameterisation
% FORMAT J = spm_ness_J(P,M,X)
%--------------------------------------------------------------------------
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness_hd.m 8089 2021-04-04 12:12:06Z karl $

% numerical evaluation of Jacobian
%--------------------------------------------------------------------------
J = full(spm_diff('spm_NESS_gen',P,M,X,3));

return
