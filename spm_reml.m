function [Ce,h,W] = spm_reml(Cy,X,Q);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h,W] = spm_reml(Cy,X,Q);
%
% Cy  - (m x m) data covariance matrix y*y'
% X   - (m x p) design matrix
% Q   - {1 x q} covariance constraints
%
% Ce  - (m x m) estimated errors = h(1)&Q{1} + h(2)*Q{2} + ...
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

% REML	objective function = r'*iCe*r + log(det(Ce)) + log(det(XiCeX));
%---------------------------------------------------------------------------
iCe   = speye(m,m);
W     = zeros(q,q);
u     = zeros(q,1);
H     = zeros(q,1);
for k = 1:128


	% E-step
	%------------------------------------------------------------------
	iCeX  = iCe*X;
	Cb    = inv(X'*iCeX);
	R     = iCe - iCeX*Cb*iCeX';

	% M-step
	%------------------------------------------------------------------
	RCy   = full(R)*Cy;
	for i = 1:q
		RQ{i}  = R*Q{i};
		u(i)   = sum(sum(RCy.*RQ{i}));		% trace{R*Cy*R*Q{i}}
	end

	for i = 1:q
	for j = 1:q
		W(i,j) = sum(sum(RQ{j}.*RQ{i}));	% trace{R*Q{j}*R*Q{i}}
	end
	end
	h     = pinv(W)*u;

	% convergence
	%------------------------------------------------------------------
	dh    = (H - h)'*(H - h);
	if dh < 1e-6, break, end
        fprintf('%-30s: %i %30s%e\n','REML Iteration',k,'...',full(dh));
	H     = h;

	% esitmate of Ce
	%------------------------------------------------------------------
	Ce    = sparse(m,m);
	for i = 1:q
		Ce = Ce + h(i)*Q{i};
	end
	iCe   = inv(Ce);

end
