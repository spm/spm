function [d] = spm_kl_gamma(b_q,c_q,b_p,c_p)
% KL divergence between two Gamma densities
% FORMAT [d] = spm_kl_gamma(b_q,c_q,b_p,c_p)
%
% KL (Q||P) = <log Q/P> where avg is wrt Q
%
% b_q, c_q    Parameters of first Gamma density
% b_p, c_p    Parameters of second Gamma density
%__________________________________________________________________________

% Will Penny 
% Copyright (C) 2004-2022 Wellcome Centre for Human Neuroimaging


digamma_c_q = psi(c_q);
d = (c_q-1)*digamma_c_q - log(b_q) - c_q - gammaln(c_q);
d = d + gammaln(c_p) + c_p*log(b_p) - (c_p-1)*(digamma_c_q+log(b_q));
d = d + b_q * c_q / b_p;
