function [f] = mci_phase_fx (x,u,P,M)
% Flow function for phase model
% FORMAT [f] = mci_phase_fx (x,u,P,M)
%
% x      state vector
% u      inputs
% P      parameter vector
% M      model structure
%
% f      dx/dt
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_phase_fx.m 6548 2015-09-11 12:39:47Z will $

params = spm_unvec (P,M.pE);

a = params.sinCoeff;
b = params.cosCoeff;
w = params.intrinPhase;

D = M.n;
% initialise state matrix
f = zeros(D,1);

for i = 1:D,
    f(i)=w(i);
    for j = 1:D,
        if ~(i==j)
            f(i) = f(i) + a(i,j)*sin(x(i)-x(j)) + b(i,j)*cos(x(i)-x(j));
        end
    end
end


