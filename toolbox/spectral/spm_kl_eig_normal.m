function [d] = spm_kl_eig_normal (m_q,c_q,c_p)
% KL divergence between normal densities using eigendecomposition 
% function [d] = spm_kl_eig_normal (m_q,c_q,c_p)
%
% Calculate the KL distance 
%
% KL (Q||P) = <log Q/P> where avg is wrt Q
%
% between two Normal densities Q and P where P is
% zero mean and has a diagonal covariance.
%
% m_q, c_q    Mean and covariance of first Normal density
% c_p         Covariance of second (zero-mean) Normal density
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_kl_eig_normal.m 1143 2008-02-07 19:33:33Z spm $

d=length(m_q);
m_q=m_q(:);

ldcp=0.5*sum(log(diag(c_p)));
[v,eigvals]=eig(c_q);
ldcq=-0.5*sum(log(diag(eigvals)));
Term1=ldcp+ldcq;

inv_c_p=diag(1./diag(c_p));
Term2=0.5*trace(inv_c_p*c_q)+0.5*m_q'*inv_c_p*m_q;
d=Term1+Term2-0.5*d;


