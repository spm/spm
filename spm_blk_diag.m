function [U,V] = spm_blk_diag(X)
% finds the block-diagonal structure of a matrix
% FORMAT [U,V] = spm_blk_diag(X)
% X - (m x n) matrix
% U - (m x p) matrix of indices for p partitions
% V - (p x m) matrix of indices for p partitions
%___________________________________________________________________________
% %W% Karl Friston %E%

% find block structure
%---------------------------------------------------------------------------
[m n] = size(X);
U     = [];
V     = [];
while 1
	i     = min(find(any(X,1)));
	if length(i)
		u     = any(X(:,i),2);
		while 1
			v  = any( X(u,:) ,1);
			du = any( X(:,v) ,2) | u;
			if any(u - du)
				u  = du;
			else
				U = [U u];
				V = [V;v];
				break
			end
		end
		X(u,v) = 0;
	else
		break
	end
end
