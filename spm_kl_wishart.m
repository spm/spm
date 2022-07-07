function [kl] = spm_kl_wishart(q,Q,p,P)
% KL divergence between two Wishart densities
% FORMAT [kl] = spm_kl_wishart(q,Q,p,P)
%
% Calculate KL (Q||P) = <log Q/P> where avg is wrt Q
% between two Wishart densities Q and P
%
% q,Q      Parameters of first density
% p,P      Parameters of first density
%__________________________________________________________________________

% Will Penny 
% Copyright (C) 2007-2022 Wellcome Centre for Human Neuroimaging


logdetQ = log(det(Q));
logdetP = log(det(P));

% Get dimension 
d = size(Q,1);

LqQ = 0;
for i=1:d
  LqQ = LqQ+psi((q+1-i)/2);
end  
LqQ = LqQ+d*log(2)-logdetQ;

LpP = 0;
for i=1:d
  LpP = LpP+psi((p+1-i)/2);
end  
LpP = LpP+d*log(2)-logdetP;

logZqQ = 0.5*q*d*log(2)-0.5*q*logdetQ+spm_lg_gamma(d,q/2);
logZpP = 0.5*p*d*log(2)-0.5*p*logdetP+spm_lg_gamma(d,p/2);

Qinv = inv(Q);
kl = 0.5*(q-d-1)*LqQ-0.5*(p-d-1)*LpP-0.5*q*d+0.5*q*trace(P*Qinv);
kl = kl + logZpP - logZqQ;
