function D = spm_mldivide(A, B)
% Regularised variant of mldivide(A, B) or A \ B, similar to spm_inv(A) * B
% FORMAT D = spm_mldivide(A, B)
%
% D = inv(A) * B, or if A is near singular D = inv(A + TOL*eye(size(A)) * B
% 
% where TOL is adaptively increased if necessary.
%
% This function should be preferable to spm_inv(A) * B if A is large and
% sparse or if B has few columns, since the inverse need not be explicitly
% computed (the linear system can be solved with the backslash operator).
%
% See also: spm_mrdivide
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Ged Ridgway
% $Id: spm_mldivide.m 5219 2013-01-29 17:07:07Z spm $

% A problem with this (and with original spm_inv) is that negative definite
% matrices (or those with mixed positive and negative eigenvalues) will not
% be sensibly regularised by adding a scaled identity. However, it is both
% expensive to compute the eigenvalues and also difficult to find a general
% strategy to handle a mixed-sign set of eigenvalues.

% Stop warning being displayed, but clear lastwarn so can observe occurence
sw = warning('off', 'MATLAB:nearlySingularMatrix');
lastwarn('');

D = A \ B;

if ~isempty(lastwarn)
    % Warning would have occurred, but was turned off, try regularising,
    % starting with low TOL and increasing if required...
    [m,n] = size(A);
    TOL = max(eps(norm(A, 'inf')) * max(m, n), exp(-32));
    while ~isempty(lastwarn)
        lastwarn('');
        D = (A + TOL * speye(m, n)) \ B;
        % ...increase regularisation in case warning obtained again...
        TOL = TOL * 100;
    end
end

warning(sw);
