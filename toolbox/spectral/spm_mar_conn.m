function [psig,chi2] = spm_mar_conn (mar,conn)
% Test for significance of connections
% FORMAT [psig,chi2] = spm_mar_conn (mar,conn)
%
% mar            MAR data structure (see spm_mar.m)
% conn           conn(i,j)=1 if we are testing significance
%                of connection from time series i to time
%                series j - zero otherwise
% 
% psig           significance of connection
% chi2           associated Chi^2 value
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar_conn.m 1143 2008-02-07 19:33:33Z spm $

pmat=eye(mar.p);
for i=1:mar.p,
   C(i).mat=kron(pmat(:,i),conn);
end

C_con=[];
for i=1:mar.p,
   C_con=[C_con, C(i).mat(:)];
end

w_hat=mar.wmean(:);

% Pick off all AR coefficients for that connection
% (ie. at all lags)
c_tilde=C_con'*w_hat;

% Get corresponding part of AR coeff covariance matrix
Sigma_tilde=C_con'*mar.w_cov*C_con;

chi2=c_tilde'*inv(Sigma_tilde)*c_tilde;

psig=1-chi2cdf(chi2,mar.p);
