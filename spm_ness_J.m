function J = spm_ness_J(P,M,X)
% Return the Jacobian given a polynomial parameterisation
% FORMAT J = spm_ness_J(P,M,X)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% numerical evaluation of Jacobian
%--------------------------------------------------------------------------
J = full(spm_diff('spm_NESS_gen',P,M,X,3));
