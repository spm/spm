function [Ep,Cp,S] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using a Levenburg-Marquardt/EM algorithm
% FORMAT [Ep,Cp,Ce] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%___________________________________________________________________________
% M.f  - f(x,u,P)
% M.g  - g(x,u,P)
%	x - state  variables
%	u - inputs or causes
%	P - free parameters
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
%
%
% U.u  - inputs
% U.dt - sampling interval
%
% Y.y  - outputs
% Y.X0 - Confounds or null space
% Y.dt - sampling interval for outputs
% Y.Ce - error covariance contraints
%
%
% Parameter estimates
%---------------------------------------------------------------------------
% Ep  - (p x 1)        	conditional expectation  E{P|y}
% Cp  - (p x p)        	conditional covariance   Cov{P|y}
% Ce  - (v x v)        	[Re]ML estimate of       Cov{e}
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% dynamic MIMO input-state-ouput model under Gaussian assumptions
%
%  		 	dx/dt = f(x,u,P)
%   		 	y     = g(x,u,P) + e
%
% A static nonlinear observation model with fixed input or causes u
% obtains when x = []. i.e.
%
%   		 	y     = g(u,P) + e
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expecation and covariance of P
% If the free-energy starts to increase, a Levenburg-Marquart scheme is
% invoked.  The M-Step estimates the covariance components of e, in terms
% of [Re]ML point estimators.
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% data y
%---------------------------------------------------------------------------
y      = Y.y;
[ns m] = size(y);			% number of samples and responses

% covariance components
%---------------------------------------------------------------------------
if isfield(Y,'Ce')
	Q = Y.Ce;
else
	Q = spm_Ce(ns*ones(1,m));
end
nh    = length(Q);                      % number of variance components
ne    = length(Q{1});                   % number of error terms

% initialise hyperparameters and precision iS = inv(S)
%---------------------------------------------------------------------------
h     = sparse(nh,1);
for i = 1:nh
     h(i) = any(diag(Q{i}));
end
S     = sparse(ne,ne);
for i = 1:nh
    S = S + h(i)*Q{i};
end
iS    = inv(S);

% priors
%---------------------------------------------------------------------------
pE     = M.pE;
pC     = M.pC;

% confounds (if specified)
%---------------------------------------------------------------------------
if isfield(Y,'X0')
	Ju = kron(speye(m,m),Y.X0);
else
	Ju = [];
end

% SVD of prior covariance
%---------------------------------------------------------------------------
Vp     = spm_svd(pC,1e-16);
np     = size(Vp,2);		% number of parameters (effective)
nu     = size(Ju,2);		% number of parameters (confounds)
ip     = [1:np];
iu     = [1:nu] + np;

% treat confounds as fixed effects
%---------------------------------------------------------------------------
uE     = inv(Ju'*Ju)*Ju'*y(:);
uC     = speye(nu,nu)*1e6;


% compbine priors on free parameters and confounds
%---------------------------------------------------------------------------
p      = [sparse(np,1); uE];
pC     = blkdiag(Vp'*pC*Vp, uC);
ipC    = inv(pC);

% EM 
%===========================================================================
dv     = 1e-6;
lm     = 0;
F      = -Inf;
dFdh   = sparse(nh,1);
dFdhh  = sparse(nh,nh);
for  k = 1:32

        % integrate
        %-------------------------------------------------------------------
        p0     = pE + Vp*p(ip);
        fp     = spm_int(p0,M,U,ns);

	% compute partial derivatives [Jp] df(p)/dp*Vp {p = parameters}
	%-------------------------------------------------------------------
	for  i = ip
		pi      = p0 + Vp(:,i)*dv;
		dfp     = spm_int(pi,M,U,ns) - fp;
		Jp(:,i) = dfp(:)/dv;
	end

        % e and dedp = J
	%-------------------------------------------------------------------
        e      = y(:) - fp(:) - Ju*p(iu);
        J      = - [Jp Ju];				

	% M-Step: ReML estimator of variance components:  h = max{F(p)}
	%===================================================================
        for i = 1:8

            % precision and conditional covariance
            %-------------------------------------------------------------------
            S     = sparse(ne,ne);
            for i = 1:nh
                S = S + h(i)*Q{i};
            end
            iS    = inv(S);
            Cp    = inv(J'*iS*J  + ipC);

            % precision operators for M-Step
            %-------------------------------------------------------------------
            for i = 1:nh
                PS{i} = -iS*Q{i};
                P{i}  = PS{i}*iS;
                PJ{i} = P{i}*J;
            end

            % derivatives: dLdh = dL/dh,...
            %---------------------------------------------------------------
            for i = 1:nh
                dFdh(i,1)      =  trace(PS{i})/2 - e'*P{i}*e/2 ...
                                 -sum(sum(Cp.*(J'*PJ{i})))/2;
                for j = i:nh
                    dFdhhij    = -sum(sum(PS{i}.*PS{j}))/2 ...
                                 -sum(sum(Cp.*(PJ{i}'*S*PJ{j})));
                    dFdhh(i,j) =  dFdhhij;
                    dFdhh(j,i) =  dFdhhij;
                end
            end
                
            % update ReML estimate
            %---------------------------------------------------------------
            dh    = inv(-dFdhh)*dFdh;
            h     = h  + dh;

            % convergence
	    %---------------------------------------------------------------
            nmh   = norm(dh);
	    if nmh < 1e-5, break, end

        end

    	% objective function
	%===================================================================
        Fp    = - e'*iS*e/2 ...
                - p'*ipC*p/2 ...
                - log(det(S))/2 ...
                + log(det(Cp))/2;

    	% Levenburg Marquardt regulatisation
	%-------------------------------------------------------------------
        if Fp < F
                if ~lm
                    lm = speye(np + nu)/256;
                else
                    lm = lm*2;
                end
	        p     = p - dp;
        else

	        % update gradients, curvatures and F
	        %----------------------------------------------------------
                dFdp  = -J'*iS*e - ipC*p;
                dFdpp = -J'*iS*J - ipC;
                F     = Fp;

        end

	% E-Step: Conditional estimator of new expansion point E{p|y}
	%===================================================================
        dp    = inv(-dFdpp + lm*norm(full(dFdpp)))*dFdp;
        p     = p + dp;

	% convergence
	%-------------------------------------------------------------------
	nmp   = dp(ip)'*dp(ip);
        fprintf('%-6s: %i %6s %e %6s %e\n','E-Step',k,'F:',F,'dp:',nmp)
	if nmp < 1e-5, break, end

	% graphics
	%-------------------------------------------------------------------
	if length(dbstack) < 3
		bar(pE + Vp*p(ip))
		xlabel('parameter')
		ylabel('conditional expectation')
		title(sprintf('%s: %i','E-Step',j))
		grid on
		drawnow
	end
end

% outputs
%---------------------------------------------------------------------------
Ep     = pE + Vp*p(ip);
Cp     = Vp*Cp(ip,ip)*Vp';










