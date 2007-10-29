function [x] = spm_expm(J,x)
% approximate matrix exponential using a Taylor expansion
% FORMAT [y] = spm_expm(J,x)
% FORMAT [y] = spm_expm(J)
% y          = expm(J)*x:
% y          = expm(J);
%
% This routine covers and extends expm  functionality  by  using  a
% comoutationally  expedient  approximation  that can handle sparse
% matrices when dealing with the special case of expm(J)*x, where x
% is a vector, in an efficient fashion
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_expm.m 984 2007-10-29 19:36:49Z karl $


% expm(J) use Pade approximation
%---------------------------------------------------------------------------
if nargin == 1

	% ensure norm is < 1/2 by scaling by power of 2
	%-------------------------------------------------------------------
    J     = sparse(J);
	[f,e] = log2(norm(J,'inf'));
	s     = max(0,e+1);
	J     = J/2^s;

	X     = J; 
	c     = 1/2;
	E     = speye(size(J)) + c*J;
	D     = speye(size(J)) - c*J;
	q     = 6;
	p     = 1;
	for k = 2:q
		c   = c * (q-k+1) / (k*(2*q-k+1));
 		X   = J*X;
 		cX  = c*X;
		E   = E + cX;
		if p
			D = D + cX;
		else
			D = D - cX;
		end
		p = ~p;
	end
	E = D\E;

	% Undo scaling by repeated squaring
	%-------------------------------------------------------------------
	for k = 1:s
		E = E*E;
	end
	x     = E;
    
else

	% compute y = expm(J)*x = (1 + J + J*J/2! + J*J*J/3!  + ...)*x
	%-------------------------------------------------------------------
	J     = sparse(J);
	x     = sparse(x);
	x0    = x;
	fx    = J*x;
	j     = 1;

	while norm(fx,1) > 1e-16

		% accumulate high-order terms until convergence
		%-----------------------------------------------------------
		j  = j + 1;
		x  = x + fx;
		fx = J*fx/j;

		% revert to Pade approximation if numerical overflow
		%-----------------------------------------------------------
		if norm(x,1) > 1e16
			x = spm_expm(J)*x0;
			return
		end
	end
end
