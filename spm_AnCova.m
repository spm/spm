function [df,F,BETA,T,RES] = spm_AnCova(C,G,sigma,X,CON)
% AnCova for temporally correlated variables
% FORMAT [df,F,BETA,T,RES] = spm_AnCova(C,G,sigma,[X,[CON]]);
%
% C     - Design matrix partition of interest (K*c)
% C     - Design matrix partition of confounds (K*g)
% sigma - SD of Gaussian kernel K
% X     - {m x n} matrix of response {m x 1} variables (K*x)
% CON   - {t x p} matrix of contrasts p = size([C G],2)
%
% df    - {1 x 2} vector of degrees of freedom
% F     - {1 x n} vector of F statistics
% BETA  - {p x n} matrix of parameter estimates
% T     - {t x n} matrix of T values
% RES   - {1 x n} matrix of residual SSQ
%_______________________________________________________________________
%
% spm_AnCova uses a General[ized] Least quares Model of the form
%
%	K*x  = K*[c g]*BETA + K*e
%
% to compute the parameter estimates (BETA) and make inferences (F)
% where K is a convolution matrix corresponding to a Gaussian kernel
% of parameter sigma.  sigma = 0 corresponds to the case of independent
% observations
%
% ref: Seber GAF (1977) Linear regression analysis, Wiley New York
% ref: Worsley KJ & Friston KJ (1995) NeuroImage 2:173-181
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% check design matrix.  If C = [] then assume statistics are required
% for the effect of confounds
%---------------------------------------------------------------------------
q     = size([C G],1);
if size(C,2) == 0; C = G; G = [];  end
H     = [C G];
if size(G,2) == 0; G = zeros(q,1); end


% temporal convolution kernel
%---------------------------------------------------------------------------
if sigma
	k = exp(-[(-6*sigma):(6*sigma)].^2/(2*sigma^2));
	k = k/max(k);
else
	k = 1;
end
v     = conv(k,k);
vv    = conv(v,v);


% cunning piece of code to compute expectations of variances and effective
% degrees of freedom by capitalizing on the sparsity structure of the problem
%---------------------------------------------------------------------------
PH     = pinv(H);
PHV    = PH;
for i  = 1:size(PHV,1)
	d        = conv(PHV(i,:),v);
	PHV(i,:) = d([1:q] + (length(v) - 1)/2);
end
PG     = pinv(G);
PGV    = PG;
for i  = 1:size(PGV,1)
	d        = conv(PGV(i,:),v);
	PGV(i,:) = d([1:q] + (length(v) - 1)/2);
end
PHVV   = PHV;
for i  = 1:size(PHVV,1)
	d         = conv(PHVV(i,:),v);
	PHVV(i,:) = d([1:q] + (length(v) - 1)/2);
end

% traces
%---------------------------------------------------------------------------
d      = conv(ones(1,q),v.^2);
trV    = v((length(v) + 1)/2)*q;
trVV   = sum(d([1:q] + (length(v) - 1)/2));

trQ    = sum(sum(H'.*PHV));
trP    = sum(sum(G'.*PGV));
trQV   = sum(sum(H'.*PHVV));

trQQ   = trace(H*(PHV*H)*PHV);
trPP   = trace(G*(PGV*G)*PGV);
trPQ   = trace(G*(PGV*H)*PHV);

trRV   = trV - trQ;
trRV2  = trVV - 2*trQV + trQQ;
trR0V  = trQ - trP;
trR0V2 = trQQ - 2*trPQ + trPP;


% the [effective] degrees of freedom
%---------------------------------------------------------------------------
df    = [trR0V^2/trR0V2 trRV^2/trRV2];

% if a response variable is specified
%---------------------------------------------------------------------------
T     = [];
F     = [];
BETA  = [];
if nargin > 3

	% estimate parameters and sum of squares due to error
	%-------------------------------------------------------------------
	BETA  = (PH*X);					% parameter estimates
	R     = X - H*BETA;				% residuals under Ha
	N     = X - G*(PG*X);				% residuals under Ho
	RES   = sum(R.^2);				% SSQ under Ha
	NUL   = sum(N.^2);				% SSQ under Ho

	% test for effects with he F ratio of variances
	%-----------------------------------------------------------------
	F     = trRV/trR0V*(NUL - RES)./RES;		% F statistic

	% compute t statistics if contrasts are specified
	%-----------------------------------------------------------------
	if nargin > 4
		T     = zeros(size(CON,1),size(BETA,2));
		for i = 1:size(CON,1)
			c      = CON(i,:);
			RMS    = sqrt((RES/trRV)*(c*PHV*PH'*c'));
			T(i,:) = c*BETA./RMS;		% T statistic
		end
	end
end
