function [x] = spm_expm(J,x,n)
% approximate matrix exponential using a Taylor expansion
% FORMAT [y] = spm_expm(J,x,n)
% y          = expm(J)*x:
%
% (n - optional parameter to deal with oveflow situations)
%
% This routine covers and extends expm  functionality  by  using  a
% comoutationally  expedient  approximation  that can handle sparse
% matrices when dealing with the special case of expm(J)*x, where x
% is a vector, in an efficient fashion
%___________________________________________________________________________
% %W% Karl Friston %E%


% compute expm(J/n) {expm(J)*x = expm(J/n).....expm(J/n)*x
%---------------------------------------------------------------------------
if nargin == 2
	n = 1;
end
J     = sparse(J/n);
x     = sparse(x);
for i = 1:n

	% compute y  = expm(J)*x
	%            = (1 + J + J*J/2! + J*J*J/3!  + ...)*x
        %-------------------------------------------------------------------
	fx    = J*x;
	for j = 1:1024

		% acculuate high order terms
		%-----------------------------------------------------------
                x  = x + fx;
                fx = J*fx/(j + 1);

		% break if converging
		%-----------------------------------------------------------
		m  = max(abs(fx));
		if m > 1e+8
			warning('results may be inaccurate')
		end
		if m < 1e-8
			break
		end
	end

	% break if converging
	%-----------------------------------------------------------
	if j == 1024
		warning('number of iterations exceeded')
	end
end
