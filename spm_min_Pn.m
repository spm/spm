function x = min_Pn(n, p,s,u,V)
% A JB special
% FORMAT
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

x = (spm_Pn(n,s,u,V) - p).^2;

