function [J] = spm_diff(varargin)
% matrix differential
% FORMAT [dfdx] = spm_diff(f,x,...,P,n)
% 
% f   - [inline] function f(x,P)
% x   - argument[s]
% P   - parameter[s]
% n   - argument or parameter to differentiate w.r.t.
%
% dfdx - df(x,P)/dx{n}
%___________________________________________________________________________
% %W% Karl Friston %E%


% create inline object
%---------------------------------------------------------------------------
f     = fcnchk(varargin{1});
x     = varargin(2:(end - 1));
n     = varargin{end};
dx    = 1e-6;

if length(n) == 1

	% dfdx
	%------------------------------------------------------------------
	f0    = feval(f,x{:});
	J     = sparse(length(f0(:)),length(x{n}(:)));
	if size(J,1)
		for i = 1:length(x{n}(:))
			xi       = x;
			xi{n}(i) = xi{n}(i) + dx;
			dfdx     = (feval(f,xi{:}) - f0)/dx;
			J(:,i)   = sparse(dfdx(:));
		end
	end
else
	% dfdxdx
	%------------------------------------------------------------------
	f0    = spm_diff(f,x{:},n(1));
	J     = cell(1,length(x{n(2)}(:)));
	for i = 1:length(x{n(2)}(:))
		xi          = x;
		xi{n(2)}(i) = xi{n(2)}(i) + dx;
		dfdx        = (spm_diff(f,xi{:},n(1)) - f0)/dx;
		J{i}        = sparse(dfdx);
	end
end
