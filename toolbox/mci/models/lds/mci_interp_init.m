function [w0] = mci_interp_init (Y,M)
% Linear interpolate to t=0
% FORMAT [w0] = mci_interp_init (Y,M)
%
% Y     Cell of data from multiple subjects
%       Y{n}.y, Y{n}.ind for n=1..N
% M     Model structure
%
% w0    [d x N] matrix of initial states
%       where d is number of states
%__________________________________________________________________________

% Will Penny
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

N=length(Y);
d=size(Y{1}.y,2);
doplot=0;

for n=1:N,
    for j=1:d,
        
        % Fit
        Nt=size(Y{n}.y,1);
        xd=[Y{n}.ind(:),ones(Nt,1)];
        yd=Y{n}.y(:,j);
        beta=pinv(xd)*yd;
        
        % Extrapolate
        xt=[[1:M.N]',ones(M.N,1)];
        yhat=xt*beta;
        
        if doplot
            figure;plot(Y{n}.ind,yd,'x');
            hold on; plot([1:M.N]',yhat,'r');
        end
        w0(j,n)=yhat(1);
    end
end