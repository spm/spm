function [x] = spm_expm(J,x)
% approximate matrix exponential using a Taylor expansion
% FORMAT [y] = spm_expm(J,x)
% y          = expm(J)*x:
%
% This routine covers and extends expm  functionality  by  using  a
% comoutationally  expedient  approximation  that can handle sparse
% matrices when dealing with the special case of expm(J)*x, where x
% is a vector, in an efficient fashion
%___________________________________________________________________________
% %W% Karl Friston %E%


% compute y  = expm(J)*x
%            = (1 + J + J*J/2! + J*J*J/3!  + ...)*x
%---------------------------------------------------------------------------
J     = sparse(J);
x     = sparse(x);
fx    = J*x;
j     = 1;
while norm(fx,1) > 1e-16

	% acculuate high order terms
	%-------------------------------------------------------------------
	j  = j + 1;
	x  = x + fx;
	fx = J*fx/j;

end
