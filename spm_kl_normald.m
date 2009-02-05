function [d] = spm_kl_normald (m_q,c_q,m_p,c_p)
% KL divergence between two Gaussians with, possibly, diagonal covariances
% FORMAT [d] = spm_kl_normald (m_q,c_q,m_p,c_p)
%
% Calculate the KL distance KL (Q||P) = <log Q/P> where avg is wrt Q
% between two Normal densities Q and P
%
% m_q, c_q    Mean and covariance of first Normal density
% m_p, c_p    Mean and covariance of second Normal density
% 
% If c_q and c_p are diagonal, pass them as vectors, and KL will
% be computed more efficiently. Both must be full or both must be diagonal
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kl_normald.m 2696 2009-02-05 20:29:48Z guillaume $

dd=length(m_q);

m_q=m_q(:);
m_p=m_p(:);

[d1,d2]=size(c_q);

if (d1==d2)
    % Full covariances
    Term1=0.5*spm_logdet(c_p)-0.5*spm_logdet(c_q);
    inv_c_p=inv(c_p);
    Term2=0.5*trace(inv_c_p*c_q)+0.5*(m_q-m_p)'*inv_c_p*(m_q-m_p);
else
    % Diagonal covariances
    c_q=c_q(:);
    c_p=c_p(:);
    Term1=0.5*sum(log(c_p))-0.5*sum(log(c_q));
    Term2=0.5*sum(c_q./c_p);
    e=sqrt(1./c_p).*(m_q-m_p);
    Term2=Term2+0.5*e'*e;
end
d=Term1+Term2-0.5*dd;
