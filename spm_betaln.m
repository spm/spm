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
% $Id: spm_betaln.m 7508 2018-12-21 09:49:44Z thomas $

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
                    for m = 1:size(z,6)
                        y(1,i,j,k,l,m) = spm_betaln(z(:,i,j,k,l,m));
                    end
                end
            end
        end
    end
end
