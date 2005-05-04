function [W] = spm_Volt_W(u)
% returns basis functions used for Volterra expansion
% FORMAT [W] = spm_Volt_W(u);
% u  - times {seconds}
% W  - basis functions (mixture of Gammas)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_Volt_W.m 112 2005-05-04 18:20:52Z john $


u     = u(:);
W     = [];
for i = 2:4
	m   = (2^i);
	s   = sqrt(m);
	W   = [W spm_Gpdf(u,(m/s)^2,m/s^2)];
end
