function [BASIS] = spm_basis(m,n,h)
% Creates basis functions for stereotactic normalization
% FORMAT [BASIS] = spm_basis(m,n,[h]);
% m   - x dimension
% n   - y dimension
% h   - order
%____________________________________________________________________________
%
% spm_basis creates a series of two dimensional basis functions using
% Fourier-like modes and post hoc orthogonalization.  The order {h}
% specifies the number of modes (of increasing spatial frequency) =
% h^2 + (h + 1)^2
%
% The matrix of basis functions (BASIS) returned by spm_basis are usually
% used in nonlinear spatial normalization problems (see spm_sn.m)
%
%__________________________________________________________________________
% %W% %E%

% create basis functions (BASIS)
%----------------------------------------------------------------------------
if nargin ~= 3; h = 2; end

BASIS = [];
d     = min([m n]);

if d > 1
	x     = pi*[0:(m - 1)]/(m - 1);
	y     = pi*[0:(n - 1)]/(n - 1);
	for i = 1:h
   		for j = 1:h
      			s = sin(x*i)'*sin(y*j);
      			BASIS = [BASIS s(:)];
    		end
	end
	for i = 0:h
   		for j = 0:h
      			s = cos(x*i)'*cos(y*j);
      			BASIS = [BASIS s(:)];
  		 end
	end
else
	m = max([m n]);
	x     = pi*[0:(m - 1)]/(m - 1);
	for i = 0:h
   		for j = 0:h
      			s = cos(x*i);
      			BASIS = [BASIS s(:)];
  		end
	end
end

% orthogonalize and scale to unit maximum
%----------------------------------------------------------------------------
BASIS = orth(BASIS);
BASIS = BASIS/max(BASIS(:));

