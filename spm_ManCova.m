function [CV,CU,CS,CHI,P,BETA] = spm_ManCova(C,G,X)
% ManCova for the multivariate linear model
% FORMAT [CV,CU,CS,CHI,P,BETA] = spm_ManCova(C,G,X);
%
% C     - {m x h} Design matrix partition of interest
% G     - {m x g} Design matrix partition of confounds
% X     - {m x n} matrix of response {m x 1} variables
%
% CHI   - {1 x h} vector of Chi squared statistics
% P     - {1 x h} vector of P values
% BETA  - {p x n} matrix of parameter estimates
% CV    - {m x p} matrix of canonical variates
% CU    - {n x p} matrix of canonical vectors
% CS    - {1 x p} vector of canonical values
%_______________________________________________________________________
%
% spm_ManCova uses a General Least Squares Model of the form
%
%	X  = [C G]*BETA + e
%
% to compute the parameter estimates (BETA), make multivariate
% inferences (CHI & P) and characterize the responses in terms of Canonical
% Vectors and Variates (CU & CV).  If the size(X,2) > size(X,1)/3 then
% an SVD dimension reduction to a p-variate is used.
%
% This analysis is the same as a canonical correlation analysis.  In this
% instance the canonical vectors associated with C are given by 
% BETA*CU such that CV = C*BETA*CU.
%
% C is orthogonalized with respect to G and the confounds are removed
% prior to SVD.
%
% Orthogonalization of the confounds is required to use the inferences
% about the dimensionality of the response. i.e. P(s) is the P{dimension-
% ality of the response > (s - 1)}.  P(s) is based on CHI(s).  CHI(1)
% corresponds to Wilk's Lambda (after transformation).
%
%_______________________________________________________________________
% %W% Karl Friston %E%

% check design matrix.  If C = [] then assume statistics are required
% for the effect of confounds
%-----------------------------------------------------------------------
q     = size([C G],1);
if size(C,2) == 0; C = G; G = [];  end
if size(G,2) == 0; G = zeros(q,1); end


% orthogonalize C w.r.t. G
%-----------------------------------------------------------------------
C     = C - G*(pinv(G)*C);

% remove confounds
%-----------------------------------------------------------------------
XA    = X - G*(pinv(G)*X);


% MULTIVARIATE ANALYSIS - ManCova
%=======================================================================

% SVD DIMENSION REDUCTION
%-----------------------------------------------------------------------
[e s u] = svd(XA,0);
v       = diag(s).^2;
v       = length(v)*v/sum(v);

% reduce dimensionality (if required)
%-----------------------------------------------------------------------
if size(XA,2) > size(XA,1)/3;

	Q     = find(v > 1);
else
	Q     = find(v > eps);
end
e     = e(:,Q);
u     = u(:,Q);
s     = s(Q,Q);
X     = e*s;


% MANCOVA proper
%=======================================================================

% degrees of freedom
%-----------------------------------------------------------------------
h     = rank(C);					% condition
r     = q - h - rank(G);				% residuals
p     = size(X,2);					% p-variate repsonse

% estimate parameters
%-----------------------------------------------------------------------
BETA  = pinv(C)*X;					% parameter estimates
T     = C*BETA;						% condition
R     = X - T;						% error


% CANONICAL VARIATES ANALYSIS
%=======================================================================

% Using the generalized eigenvalue solution
%-----------------------------------------------------------------------
[E V] = eig(T'*T,R'*R);
V     = diag(real(V));
[d j] = sort(-V);
V     = V(j);
E     = real(E(:,j));

% test for the dimensionality of the alternative hypothesis
%-----------------------------------------------------------------------
P     = [];
CHI   = [];
for t = 0:(min([h p]) - 1)
	d    = (r - ((p - h + 1)/2))*log(prod(1 + V((t + 1):p)));
	P    = [P (1 - spm_Xcdf(d,(p - t)*(h - t)))];
	CHI  = [CHI d];
end

% compute canonical images {vectors} and variates
%-----------------------------------------------------------------------
CV    = T*E;					% canonical variate
CU    = u*E;					% canonical vector
CS    = V';					% canonical values
BETA  = BETA*u';				% parameter estimates

