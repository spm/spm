function [D] = spm_KL_cat(Q,P)
% KL divergence between two categorical distributions
% FORMAT [D] = spm_KL_cat(Q,P)
%
% Calculates KL(Q||P) = <log Q/P> wrt Q between two categorical 
% distributions Q and P
%
% if supplied with arrays, the KL divergence will be summed over the first
% dimension. The arrays can be normalised (c.f., Dirichlet parameters).
%
% see also: spm_kl_dirichlet.m (for row vectors)
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% Karl Friston 
% $Id: spm_KL_dir.m 8183 2021-11-04 15:25:19Z guillaume $

% ensure sum to one constraint over first mention
%--------------------------------------------------------------------------
Q = spm_norm(Q(:,:));
P = spm_norm(P(:,:));

%  KL divergence
%--------------------------------------------------------------------------
D = sum( sum(Q.*(spm_log(Q) - spm_log(P))) );


return

function A  = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
A           = bsxfun(@rdivide,A,sum(A,1));
A(isnan(A)) = 1/size(A,1);


function A  = spm_log(A)
% log of numeric array plus a small constant
%--------------------------------------------------------------------------
A           = log(A);
A(isinf(A)) = -32;