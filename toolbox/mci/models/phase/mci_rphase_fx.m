function [f] = mci_rphase_fx (x,u,P,M)
% Flow function for phase model
% FORMAT [f] = mci_rphase_fx (x,u,P,M)
%
% x      state vector
% u      inputs
% P      parameter vector
% M      model structure
%
% f      dx/dt
%__________________________________________________________________________

% Will Penny and Biswa Sengupta
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

params = spm_unvec (P,M.pE);

aconn = params.aconn;
bconn = params.bconn;

Nc=length(M.conn);
a=zeros(M.n,M.n);
b=zeros(M.n,M.n);
for c=1:Nc,
    i=M.conn(c).regions(1);
    j=M.conn(c).regions(2);
    a(i,j)=aconn(c);
    b(i,j)=bconn(c);
end

for i = 1:M.n,
    f(i)=M.freq;
    for j = 1:M.n,
        if ~(i==j)
            f(i) = f(i) + a(i,j)*sin(x(i)-x(j)) + b(i,j)*cos(x(i)-x(j));
        end
    end
end
f=2*pi*f(:);

