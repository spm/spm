function [A,B,C,H] = spm_dcm_reshape(P,m,n,r)
% converts free parameter vector to matrices
% FORMAT [A B C H] = spm_dcm_reshape(P,m,n,[r]);
% P     - parameter vector
% m     - number of inputs
% n     - number of states
% [r]   - returns relative connections {without scaling by P(1)}
%
% A...  - intrinsic connections
% B     - modulatory connections
% C     - direct connections
% H     - hemodynamic parameters
%___________________________________________________________________________
% @(#)spm_dcm_reshape.m	2.2 Karl Friston 02/02/22

% scale intrinsic connections {A}
%---------------------------------------------------------------------------
P     = full(P);
if nargin == 4
	q  = 1;
else
	q  = P(1);
end
P(1)  = [];

% fill in intrinsic connections {A}
%---------------------------------------------------------------------------
j     = 1:n*n;
A     = reshape(P(j),n,n)*q;
P(j)  = [];

% fill in modulatory connections {B}
%---------------------------------------------------------------------------
j     = 1:n*n*m;
B     = reshape(P(j),n,n,m)*q;
P(j)  = [];

% fill in direct connections {C}
%---------------------------------------------------------------------------
j     = 1:n*m;
C     = reshape(P(j),n,m);
P(j)  = [];

% fill in hemodynamic parameters {H}
%---------------------------------------------------------------------------
j     = 1:n*5;
H     = reshape(P(j),n,5);
