function [Omega] = spm_get_omega (p,d,w_cov,xtx)
% Get expected error of MAR model
% FORMAT [Omega] = spm_get_omega (p,d,w_cov,xtx)
%
% p         Number of time lags
% d         Dimension of time series
% w_cov     Uncertainty in MAR coefficients
% xtx       X'X where X is design matrix (ie. from lagged data)
%
% Omega     Expected error
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


Omega=zeros(d,d);
% Submatrix size
s=p*d;
for di=1:d
    for dj=di:d
        istart=1+(di-1)*s;
        istop=istart+s-1;
        jstart=1+(dj-1)*s;
        jstop=jstart+s-1;
        w_cov_i_j=w_cov(istart:istop,jstart:jstop);
        Omega(di,dj)=trace(w_cov_i_j*xtx);
    end
end

Omega = Omega+Omega'-diag(diag(Omega));
