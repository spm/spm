function x = min_Pz(z, p,s,V)
% A JB special
% FORMAT
%_______________________________________________________________________
% %W% Jean-Baptiste Poline %E%

x = (spm_Pz(s,z,V) - p).^2;

