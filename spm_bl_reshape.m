function [A,B,C,D,L,O] = spm_bl_reshape(P,m,n,l)
% embeds free parameter vector in Lie operators for a Bilinear model
% FORMAT [A,B,C,D,L,O] = spm_bl_reshape(P,m,n,l);
% P     - parameter vector
% m     - number of inputs
% n     - number of states
% l     - number of outputs
%
% A...  - Lie operators (matices)
% L     - 1st order output matrix
% O     - 2nd order output matrix
%___________________________________________________________________________
% %W% Karl Friston %E%

% get A,...,O {where P = [A(:); B(:); C(:); D(:); L(:); O(:)]}
%---------------------------------------------------------------------------
P     = full(P);
A     = reshape(P([1:n*n]),n,n);
B     = reshape(P([1:n*n*m] + n*n),n,n,m);
C     = reshape(P([1:n*m] + n*n + n*n*m),n,m);
D     = reshape(P([1:n] + n*n + n*n*m + n*m),n,1);
L     = reshape(P([1:n*l] + n*n + n*n*m + n*m + n),n,l);
O     = reshape(P([1:n*n*l] + n*n + n*n*m + n*m + n + n*l),n,n,l);
