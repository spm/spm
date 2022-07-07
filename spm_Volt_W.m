function [W] = spm_Volt_W(u)
% Return basis functions used for Volterra expansion
% FORMAT [W] = spm_Volt_W(u)
% u  - times {seconds}
% W  - basis functions (mixture of Gammas)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1997-2022 Wellcome Centre for Human Neuroimaging


u     = u(:);
W     = [];
for i = 2:4
    m = (2^i);
    s = sqrt(m);
    W = [W spm_Gpdf(u,(m/s)^2,m/s^2)];
end
