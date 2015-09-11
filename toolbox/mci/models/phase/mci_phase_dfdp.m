function [dfdp] = mci_phase_dfdp (x,u,P,M)
% Parameter sensitivity for phase model
% FORMAT [dfdp] = mci_phase_dfdp (x,u,P,M)
%
% x      State vector
% u      inputs
% P      parameter vector
% M      model structure
%
% dfdp   Jacobian wrt. parameters, df/dp
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_phase_dfdp.m 6548 2015-09-11 12:39:47Z will $

D = M.n;
df_dp=zeros(D,length(P));

for k=1:D,
    % kth state variable
    for i=1:D,
        for j=1:D,
            % i,jth parameters in a,b matrices
            if i==j
                da(i,j)=0;
                db(i,j)=0;
            else
                % derivative only non-zero for k=i
                da(i,j)=(k==i)*sin(x(i)-x(j));
                db(i,j)=(k==i)*cos(x(i)-x(j));
            end
        end
        dw(i)=(k==i);
    end
    % Concatenate in order produced by spm_vec 
    dfdp(k,:)=[db(:)',da(:)',dw(:)'];
end



