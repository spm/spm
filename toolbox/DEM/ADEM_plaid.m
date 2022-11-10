function [y,n] = ADEM_plaid(x,n)
% creates a Gaussian modulated n x n visual plaid stimulus
% FORMAT [y,n] = ADEM_plaid(x,[n])
% x(1) - horizontal displacement
% x(2) - vertical displacement
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% default stimulus
%--------------------------------------------------------------------------
try
    n;
catch
    n = 4;
end
 
% stimulus
%--------------------------------------------------------------------------
sx  = [1:n]';
sx  = sx - n/2 - n*x(1)/8;
sy  = [1:n] ;
sy  = sy - n/2 - n*x(2)/8;
sx  = exp(-sx.^2/(2*(n/6)^2)).*cos(2*pi*sx*2/n);
sy  = exp(-sy.^2/(2*(n/6)^2)).*cos(2*pi*sy*2/n);

% vectorise under defaults
%--------------------------------------------------------------------------
if nargin > 1
    y = sx*sy;
else
    y = spm_vec(sx*sy);
end
