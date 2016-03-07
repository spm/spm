function y = spm_betaln(z)
% returns the log the multivariate beta function of a vector.
% FORMAT y = spm_betaln(z)
%   y = spm_betaln(z) computes the natural logarithm of the beta function
%   for corresponding elements of the vector z. if concerned is an array,
%   the beta functions are taken over the elements of the first to mention
%   (and size(y,1) equals one).
%
%   See also BETAINC, BETA.
%--------------------------------------------------------------------------
%   Ref: Abramowitz & Stegun, Handbook of Mathematical Functions, sec. 6.2.
%   Copyright 1984-2004 The MathWorks, Inc. 
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_betaln.m 6741 2016-03-07 10:32:29Z karl $

% log the multivariate beta function of a vector
%--------------------------------------------------------------------------
if isvector(z)
    z     = z(find(z)); %#ok<FNDSB>
    y     = sum(gammaln(z)) - gammaln(sum(z));
else
    for i = 1:size(z,2)
        for j = 1:size(z,3)
            for k = 1:size(z,4)
                for l = 1:size(z,5)
                    y(1,i,j,k,l) = spm_betaln(z(:,i,j,k,l));
                end
            end
        end
    end
end
