function [F] = mci_ramsay_fx (x,U,P,M)
% State equation for Ramsay model
% FORMAT [F] = mci_ramsay_fx (x,U,P,M)
%
% x     State vector
%       x(1) Voltage variable
%       x(2) Recovery variable
% U     inputs
% P     vector of model parameters - 2 params only
% M     model
%
% F     dx/dt
%
% J Ramsay et al (2007) Parameter estimation for differential equations:
% a generalised smoothing approach. J Roy Stat Soc B, 69(5):741-796.
%
% See also section 10 (page 26) and contribution by W.Penny on page 75 of: 
%
% Girolami and Calderhead (2011) Riemann manifold Langevin and Hamiltonian
% Monte Carlo methods. J Roy Stat Soc B,73(2):123-214.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_ramsay_fx.m 6548 2015-09-11 12:39:47Z will $

c=3;

P(1)=exp(P(1));
P(2)=exp(P(2));

F=zeros(2,1);
F(1)=c*(x(1)-(1/3)*x(1)^3+x(2));
F(2)=-(1/c)*(x(1)-P(1)+P(2)*x(2));

