function [W] = spm_Volt_W(u)
% returns basis functions used for Volterra expansion
% FORMAT [W] = spm_Volt_W(u);
% u  - times {seconds}
% W  - basis functions (mixture of Gammas)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Volt_W.m 1143 2008-02-07 19:33:33Z spm $


u     = u(:);
W     = [];
for i = 2:4
    m   = (2^i);
    s   = sqrt(m);
    W   = [W spm_Gpdf(u,(m/s)^2,m/s^2)];
end
