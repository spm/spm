function [C] = spm_poly(x,b)
% returns a basis set of polynomials using a specified vector
% FORMAT [C] = spm_poly(x,b)
%
% x   - vector pf parameters
% b   - number of subjects
% C   - basis
%____________________________________________________________________________
%
% returns a basis set of polynomials using a specified vector to order
% 3 and mean corrects. 
%
%__________________________________________________________________________
% %W% %E%

C     = [];
D     = [];
x     = x - mean(x);

for i = 1:3
	D = [D x(:).^i]; 	end

for i = 1:b
	C = [C; D];		end

C     = C - ones(size(C,1),1)*mean(C);
C     = C./(ones(size(C,1),1)*sqrt(sum(C.^2)));
