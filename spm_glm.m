function varargout = spm_glm(varargin)
% ****
% FORMAT varargout = spm_DesUtil(action,varargin)
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%


%-Parameters
%-----------------------------------------------------------------------
tol = 1e-12;				%-Tolerance for integer rounding

%-Format arguments
%-----------------------------------------------------------------------
if nargin==0, error('do what? no arguments given...')
	else, action = varargin{1}; end


switch lower(action), case 'kcc'                 %-Karl's "cunning" code
%=======================================================================
% xXd = spm_glm('KCC',X,K)

%-Format arguments
if nargin<2, error('insufficient arguments'), end
X = varargin{2};
if nargin<3, K=speye(size(X,1)); else, K=varargin{3}; end

%-Cunning piece of code to compute expectations of variances and 
% effective degrees of freedom by capitalizing on the sparsity
% structure of the problem
%-----------------------------------------------------------------------

%-Pseudoinverses and related matrices
V      = K*K';
KX     = K*X;
PKX    = pinv(K*X);
PKXV   = PKX*V;
PKXVV  = PKXV*V;

%-Traces
%-----------------------------------------------------------------------
trV    = trace(V);
trVV   = sum(sum(V.*V));

trQ    = sum(sum(X'.*PKXV));
trQV   = sum(sum(X'.*PKXVV));

trQQ   = sum(sum(X'.*((PKXV*X)*PKXV)));

trRV   = trV - trQ;
trRV2  = trVV - 2*trQV + trQQ;

%-The [effective] degrees of freedom (round if nearly an integer)
%-----------------------------------------------------------------------
df     = trRV^2/trRV2;
tmp    = round(df);
if abs(df-tmp)<tol, df=tmp; end

%-Matrix for computing the standard error of the paramter estimates
Bcov   = PKXV*PKX'/trRV;

%-Return design matrix parameters in data structure xXd
%-----------------------------------------------------------------------
varargout = {struct(...
		'X',X,'K',K,'KX',KX,'df',df,'Bcov',Bcov,...
		'trRV',trRV,'PKX',PKX,'PKXV',PKXV,...
		'trQ',trQ,'trQQ',trQQ			)};

%PG     = pinv(G);
%PGV    = PG*V;
%trP    = sum(sum(G'.*PGV));
%trPP   = sum(sum(G'.*((PGV*G)*PGV)));
%trPQ   = sum(sum(G'.*((PGV*X)*PKXV)));
%trR0V  = trQ - trP;
%trR0V2 = trQQ - 2*trPQ + trPP;
%
%'trR0V',trR0V,'PG',PG,


otherwise
%=======================================================================
error('Unknown action string')



%=======================================================================
end
