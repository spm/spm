function [Ce,h,W] = spm_REML(Cy,X,Q);
% REML estimation of covariance components from Cov{y}
% FORMAT [Ce,h,W] = spm_REML(Cy,X,Q);
%
% Cy  - (m x m) data covariance matrix y*y'
% X   - (m x p) design matrix
% Q   - {1 x q} covariance constraints
%
% Ce  - (m x m) estimated errors = h(1)&Q{1} + h(2)*Q{2} + ...
% h   - (q x 1) hyperparameters
%___________________________________________________________________________
% %W% John Ashburner, Karl Friston %E%

flops(0)

% order
%---------------------------------------------------------------------------
m     = length(Cy);
q     = length(Q);

% REML	objective function = r'*iCe*r + log(det(Ce)) + log(det(XiCeX));
%---------------------------------------------------------------------------
iCe   = speye(m,m);
W     = zeros(q,q);
u     = zeros(q,1);
H     = ones(q,1);
for k = 1:128

	% E-step
	%------------------------------------------------------------------
	iCeX  = iCe*X;
	Cb    = inv(X'*iCeX);
	R     = iCe - iCeX*Cb*iCeX';

	% M-step    NB trace(AB) = sum(sum(A.*B'))
	%------------------------------------------------------------------
	RCy   = R*Cy;
	for i = 1:q
		RQ{i}  = R*Q{i};
		u(i)   = sum(sum(RCy.*RQ{i}'));
	end

	for i = 1:q
	for j = 1:q
		W(i,j) = sum(sum(RQ{i}.*RQ{j}'));
	end
	end
	h     = pinv(W)*u;

	% convergence
	%------------------------------------------------------------------
	dh    = H - h;
	if dh'*dh < 1e-16, break, end
	H     = h;

	% esitmate of Ce
	%------------------------------------------------------------------
	Ce    = sparse(m,m);
	for i = 1:q
		Ce = Ce + h(i)*Q{i};
	end
	iCe   = inv(Ce);

end
f = flops;
disp([f k])
