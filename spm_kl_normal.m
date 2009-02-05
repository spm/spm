function [d] = spm_kl_normal (m_q,c_q,m_p,c_p)
% KL divergence between two multivariate normal densities
% FORMAT [d] = spm_kl_normal (m_q,c_q,m_p,c_p)
%
% KL (Q||P) = <log Q/P> where avg is wrt Q
%
% between two Normal densities Q and P
%
% m_q, c_q    Mean and covariance of first Normal density
% m_p, c_p    Mean and covariance of second Normal density
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kl_normal.m 2696 2009-02-05 20:29:48Z guillaume $

d=length(m_q);
m_q=m_q(:);
m_p=m_p(:);

Term1=0.5*spm_logdet(c_p)-0.5*spm_logdet(c_q);

inv_c_p=inv(c_p);
Term2=0.5*trace(inv_c_p*c_q)+0.5*(m_q-m_p)'*inv_c_p*(m_q-m_p);

d=Term1+Term2-0.5*d;



