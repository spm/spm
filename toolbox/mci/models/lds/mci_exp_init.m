function [w0,a] = mci_exp_init (Y,M,doplot)
% Exponentially interpolate to t=0
% FORMAT [w0,a] = mci_exp_init (Y,M,doplot)
%
% Y         Cell of data from multiple subjects
%           Y{n}.y, Y{n}.ind for n=1..N
% M         Model structure
% doplot    plot fits
%
% w0        [d x N] matrix of initial states
%           where d is number of states
% a         [d x N] matrix of exponential coefficients
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_exp_init.m 6548 2015-09-11 12:39:47Z will $

if nargin < 3
    doplot=0;
end

N=length(Y);
d=size(Y{1}.y,2);

for n=1:N,
    for j=1:d,
        % Fit
        keep=find(Y{n}.y(:,j)>0);
        ydat=Y{n}.y(keep,j);
        ind=Y{n}.ind(keep);
        Nt=size(ydat,1);
        xd=[ind(:)*M.dt,ones(Nt,1)];
        beta=pinv(xd)*log(ydat);
        
        % Extrapolate
        y0=exp(beta(2));
        yhat=y0*exp(beta(1)*M.t);
        
        if doplot
            figure;plot(ind*M.dt,ydat,'x');
            hold on; plot(M.t,yhat,'r');
        end
        w0(j,n)=y0;
        a(j,n)=beta(1);
    end
end