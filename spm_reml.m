function [Ce,h] = spm_reml(Cy,X,Q);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h] = spm_reml(Cy,X,Q);
%
% Cy  - (m x m) data covariance matrix y*y'
% X   - (m x p) design matrix
% Q   - {1 x q} covariance constraints
%
% Ce  - (m x m) estimated errors = h(1)*Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
%___________________________________________________________________________
% %W% John Ashburner, Karl Friston %E%

% order
%---------------------------------------------------------------------------
m     = length(Cy);
q     = length(Q);

% ensure X is not rank deficient
%---------------------------------------------------------------------------
X     = spm_svd(X);

% REML	objective function 2F = -r'*iCe*r - log(det(Ce)) - log(det(XiCeX));
%---------------------------------------------------------------------------
J     = zeros(q,q);
g     = zeros(q,1);
h     = zeros(q,1);
for i = 1:length(Q)
	h(i) = any(diag(Q{i}));
end
for k = 1:32


	% esitmate of Ce
	%------------------------------------------------------------------
	Ce    = sparse(m,m);
	for i = 1:q
		Ce = Ce + h(i)*Q{i};
	end
	iCe   = inv(Ce);

	% E-step {conditional covariance Cby}
	%===================================================================
	iCeX  = iCe*X;
	Cby   = inv(X'*iCeX);

	% M-step {Gauss-Newton ascent on log(p(h|y)}
	%===================================================================

	% 1st derivatives g = 2dF/dhi
	%------------------------------------------------------------------
	R     = iCe - iCeX*Cby*iCeX';
	Cy    = sparse(Cy.*spones(R));
	RCyR  = R*Cy*R';
	for i = 1:q
		QC{i}  = Q{i}*iCe;
		g(i)   = sum(sum(R.*Q{i})) - sum(sum(RCyR.*Q{i}));
	end

	% 2nd derivatives J = 2ddF/didhj
	%------------------------------------------------------------------
	for i = 1:q
	for j = i:q
		J(i,j) = sum(sum(QC{i}.*QC{j}));
		J(j,i) = J(i,j);	
	end
	end
	dh    = -J\g;
	h     = h + dh;

	% Convergence
	%===================================================================
	w     = dh'*dh;
	fprintf('%-30s: %i %30s%e\n','  ReML Iteration',k,'...',full(w));
	if w < 1e-6, break, end

end
