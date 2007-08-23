function [A,B,C,H] = spm_dcm_reshape(P,m,n,r)
% converts free parameter vector to matrices
% FORMAT [A B C H] = spm_dcm_reshape(P,m,n,[r]);
% P     - parameter vector
% m     - number of inputs
% n     - number of regions
% [r]   - returns relative connections {without scaling by P(1)}
%
% A...  - intrinsic connections
% B     - modulatory connections
% C     - direct connections
% H     - hemodynamic parameters
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_dcm_reshape.m 895 2007-08-23 21:55:21Z klaas $
 
 
% scale intrinsic connections {A}
%--------------------------------------------------------------------------
P     = full(P);
if nargin == 4
    q = 1;
else
    q = exp(P(1));
end
P(1)  = [];
 
% fill in intrinsic connections {A}
%--------------------------------------------------------------------------
j     = 1:n*n;
A     = reshape(P(j),n,n)*q;
P(j)  = [];
 
% fill in modulatory connections {B}
%--------------------------------------------------------------------------
j     = 1:n*n*m;
B     = reshape(P(j),n,n,m)*q;
P(j)  = [];
 
% fill in direct connections {C}
%--------------------------------------------------------------------------
j     = 1:n*m;
C     = reshape(P(j),n,m);
P(j)  = [];

% fill in hemodynamic parameters {H}
%--------------------------------------------------------------------------
% determine how many free hemodynamic parameters were used per region:
% 5 (classical hemodynamic model) or 6 (revised hemodynamic model)
hp    = length(P)/n;
j     = 1:n*hp;
H     = reshape(P(j),n,hp);
