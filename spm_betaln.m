function y = spm_betaln(z)
% Logarithm of the multivariate beta function of a vector
% FORMAT y = spm_betaln(z)
%   y = spm_betaln(z) computes the natural logarithm of the beta function
%   for corresponding elements of the vector z. if z is an array, the beta
%   functions are taken over the elements of the first dimension (and
%   size(y,1) equals one).
%
%   See also BETAINC, BETA.
%--------------------------------------------------------------------------
%   Ref: Abramowitz & Stegun, Handbook of Mathematical Functions, sec. 6.2.
%   Copyright 1984-2004 The MathWorks, Inc. 
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


% log multivariate beta function
%--------------------------------------------------------------------------
z     = z + exp(-16);
y     = sum(gammaln(z),1) - gammaln(sum(z,1));
