function [F,df,beta,U,V,W] = spm_mancova(xX,Y,c)
% ManCova for the multivariate linear model
% FORMAT [F,df,beta,U,V,W] = spm_mancova(xX,Y,c);
%
% X     - {m x p} Design matrix (or structure)
% Y     - {m x n} matrix of response {m x 1} variables
% c     - {u x v} contrast
%
% F     - Multivariate F values based on Wilk's Lambda
% df    - {1 x 2} vector of degrees of freedom
% beta  - {p x n} matrix of parameter estimates
% U     - {n x n} matrix of canonical vectors
% V     - {1 x n} vector of canonical values
% W     - {m x n} matrix of canonical variates
%_______________________________________________________________________
%
% spm_ManCova uses a General Least Squares Model of the form
%
%	Y  = X*beta + e
%
% to compute the parameter estimates (beta), make multivariate
% inferences (F) and characterize the responses in terms of Canonical
% Vectors and Variates (U & W).
%
% This analysis is the same as a canonical correlation analysis.  In this
% instance the canonical vectors associated with X are given by 
% beta*U such that W = X*beta*U.
%_______________________________________________________________________
% %W% Karl Friston %E%

% create design matrix structure if necessary
%---------------------------------------------------------------------------
if ~isstruct(xX)
	xX    = spm_sp('Set',xX);
end
if ~isfield(xX,'pX')
	xX.pX = spm_sp('x-',xX);
end

% contrast
%---------------------------------------------------------------------------
if ~isstruct(c)
	xCon  = spm_FcUtil('Set','','F','c',c,xX);
else
	xCon  = c;
end

% pseudoinverse of null space X0
%--------------------------------------------------------------------------
X0        = spm_FcUtil('X0',xCon,xX);
pX0       = pinv(X0);

%-Degrees of freedom
%---------------------------------------------------------------------------
h         = rank(spm_FcUtil('X1o',xCon,xX));
[q p]     = size(Y);
r         = q - xX.rk;
a         = r - (p - h + 1)/2;
if (p + h) == 3;
	b = 1;
else
	b = sqrt((p^2 * h^2 - 4)/(p^2 + h^2 - 5));
end
c         = (p*h - 2)/2;
df        = [p*h (a*b - c)];

		
%-parameter estimates
%---------------------------------------------------------------------------
Y         = Y - X0*(pX0*Y);
beta      = xX.pX*Y;
h         = xX.X*beta;
r         = Y - h;

% canonical vectors and variates
%---------------------------------------------------------------------------
[U V]     = eig(h'*h,r'*r);
V         = diag(real(V));		
[i j]     = sort(-V);
V         = V(j);					% canonical value
U         = real(U(:,j));				% canonical vector
W         = h*U;					% canonical variate

% inference based on the approximation due to Rao (1951)
%---------------------------------------------------------------------------
p         = prod(1./(1 + V))^(1/b);
F         = (1 - p)/p*df(2)/df(1);
