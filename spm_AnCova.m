function varargout = spm_AnCova(C,G,sigma,Y,CON,Xd)
% General Linear Model & inference for serially correlated regression
%
% FORMAT [Fdf,F,BETA,T,RES,BCOV,Xd] = spm_AnCova(C,G,sigma,Y,CON,Xd);
% FORMAT [Fdf,F,BETA,T,RES,BCOV,Xd] = spm_AnCova(C,G,K,Y,CON,Xd);
% FORMAT [Fdf,Xd] = spm_AnCova(C,G,sigma);
%
% C     - Design matrix partition of interest (K*C)
% G     - Design matrix partition of confounds (K*G)
% sigma - SD of Gaussian kernel K (or kernel itself) [default 0]
% Y     - {m x n} matrix of response {m x 1} variables (K*Y)
% CON   - {t x p} matrix of contrasts p = size([C G],2)
% Xd    - Data structure containing working design matrix df and pseudoinverses
%         fieldnames(Xd) = {'C';'G';'sigma';'df';'trRV';'trR0V';'PG';'PH';'PHV'}
%         (See  body of code for definitions.)
%
% Fdf   - {1 x 2} vector of degrees of freedom
% F     - {1 x n} vector of F statistics
% BETA  - {p x n} matrix of parameter estimates
% T     - {t x n} matrix of T values
% RES   - {1 x n} matrix of residual sums of squares (ResSS)
% BCOV  - {p x p} matrix such that cov(BETA) = RES*BCOV
%_______________________________________________________________________
%
% spm_AnCova uses a General Linear Model of the form:
%
%	K*Y  = K*[c g]*BETA + K*e
%
% to compute the parameter estimates (BETA) and make inferences (F)
% where K is a convolution matrix corresponding to a Gaussian kernel
% of parameter sigma.  sigma = 0 corresponds to the standard case of
% independent observations.
%
% spm_AnCova checks to see whether the previously computed df and
% pseudoinverses held in the data structure Xd are appropriate for the
% current computation, and uses them is they are. This facility is
% expedient when spm_AnCova is called repeatedly with the same design
% matrix and sigma.
% 
% ref: Seber GAF (1977) Linear regression analysis, Wiley New York
% ref: Worsley KJ & Friston KJ (1995) NeuroImage 2:173-181
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%

%-Format arguments
%-----------------------------------------------------------------------
if nargin<6, bXd=0; else, bXd=1; end
if nargin<3, sigma=0; end
if nargin<2, error('Insufficient arguments'), end

%-Output arguments
%-----------------------------------------------------------------------
F     = [];
BETA  = [];
T     = [];
RES   = [];
BCOV  = [];

%-Check whether variables in design matrix data structure Xd should be used
%-----------------------------------------------------------------------
if bXd, bXd = all([size(Xd.C)==size(C),size(Xd.G)==size(G)]); end
if bXd, bXd = all(all(Xd.C==C)) & all(all(Xd.G==G)) & (Xd.sigma==sigma); end


%-Check design matrix.  If C = [] then assume statistics are required
% for the effect of confounds
%-----------------------------------------------------------------------
q     = size([C G],1);
if size(C,2) == 0; C = G; G = [];  end
H     = [C G];
if size(G,2) == 0; G = zeros(q,1); end


%-Compute & store / unpack df and psuedoinverses
%-----------------------------------------------------------------------
if ~bXd
	%-Temporal convolution kernel - A sparse Toeplitz matrix
	%---------------------------------------------------------------
	if prod(size(sigma))==1
		K = spm_sptop(sigma,q);
	else
		K = sigma;
	end
	V      = K*K';

	%-Cunning piece of code to compute expectations of variances and 
	% effective degrees of freedom by capitalizing on the sparsity
	% structure of the problem
	%---------------------------------------------------------------
	PH     = pinv(H);
	PHV    = PH*V;
	PG     = pinv(G);
	PGV    = PG*V;
	PHVV   = PHV*V;

	%-Traces
	%---------------------------------------------------------------
	trV    = trace(V);
	trVV   = sum(sum(V.*V));

	trQ    = sum(sum(H'.*PHV));
	trP    = sum(sum(G'.*PGV));
	trQV   = sum(sum(H'.*PHVV));

	trQQ   = sum(sum(H'.*((PHV*H)*PHV)));
	trPP   = sum(sum(G'.*((PGV*G)*PGV)));
	trPQ   = sum(sum(G'.*((PGV*H)*PHV)));

	trRV   = trV - trQ;
	trRV2  = trVV - 2*trQV + trQQ;
	trR0V  = trQ - trP;
	trR0V2 = trQQ - 2*trPQ + trPP;

	%-The [effective] degrees of freedom
	%---------------------------------------------------------------
	df     = [trR0V^2/trR0V2 trRV^2/trRV2];

	%-Save design matrix parameters in data structure Xd
	%---------------------------------------------------------------
	Xd     = struct('C',C,'G',G,'sigma',sigma,'df',df,...
			'trRV',trRV,'trR0V',trR0V,'PG',PG,'PH',PH,'PHV',PHV);
else
	%-Unpack variables from data structure Xd
	%---------------------------------------------------------------
	df     = Xd.df;
	trRV   = Xd.trRV;
	trR0V  = Xd.trR0V;
	PG     = Xd.PG;
	PH     = Xd.PH;
	PHV    = Xd.PHV;
end

Fdf   = df;


%-If a response variable is specified
%-----------------------------------------------------------------------
if nargin>3

	%-Estimate parameters and sum of squares due to error
	%---------------------------------------------------------------
	BETA  = (PH*Y);					% parameter estimates

	%-SS
	%---------------------------------------------------------------
	R     = Y - H*BETA;				% residuals under Ha
	N     = Y - G*(PG*Y);				% residuals under Ho
	RES   = sum(R.^2);				% SSQ under Ha
	NUL   = sum(N.^2);				% SSQ under Ho
	
	%-Test for effects with the F ratio of variances
	%---------------------------------------------------------------
	F     = trRV/trR0V*(NUL - RES)./RES;		% F statistic

	%-Compute t statistics if contrasts are specified
	%---------------------------------------------------------------
	if nargin > 4
		T     = zeros(size(CON,1),size(BETA,2));
		for i = 1:size(CON,1)
			c      = CON(i,:);
			RMS    = sqrt((RES/trRV)*(c*PHV*PH'*c'));
			T(i,:) = c*BETA./RMS;		% T statistic
		end
	end
end

%-Matrix for computing the standard error of the paramter estimates
%-----------------------------------------------------------------------
if nargout>5, BCOV = PHV*PH'/trRV; end


%-Sort output
%-----------------------------------------------------------------------
if nargin<=3
	varargout = {Fdf,Xd};
else
	varargout = {Fdf,F,BETA,T,RES,BCOV,Xd};
end
