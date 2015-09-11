function [dfdx] = mci_rphase_dfdx (x,u,P,M)
% State sensitivity for phase model (reduced connectivity)
% FORMAT [dfdx] = mci_rphase_dfdx (x,u,P,M)
%
% x      state vector
% M      model structure
% P      parameter vector
%
% dfdx   Jacobian wrt states
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_rphase_dfdx.m 6548 2015-09-11 12:39:47Z will $

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

D = M.n;
dfdx = zeros(D,D);

ind=1:D;

for i = 1:D,
    ind_sum = [setdiff(ind,i) setdiff(i,ind)];
    for j = 1:D,
        if i == j 
            % diagonal elements
            tmp=0;
            for k = 1:length(ind_sum)
                jj=ind_sum(k);
                tmp = tmp + a(i,jj).*cos(x(i) - x(jj)) - b(i,jj).*sin(x(i) - x(jj)) ;
            end
            dfdx(i,j) = tmp;
        else
            % off-diagonal elements
            dfdx(i,j) = -a(i,j).*cos(x(i) - x(j)) + b(i,j).*sin(x(i) - x(j));
        end
    end
end

dfdx=2*pi*dfdx;

%dfdx_chk = spm_diff(M.f,x,u,P,M,1)
