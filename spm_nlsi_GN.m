function [Ep,Cp,S] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using a Fisher-Scoring/EM algorithm
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
% dynamical MIMO input-state-ouput model under Gaussian assumptions
%
%  		 	dx/dt = f(x,u,P)
%   		 	y     = g(x,u,P) + e
%
% A static nonlinear observation model with fixed input or causes u
% obtains x = []. i.e.
%
%   		 	y     = g(u,P) + e
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. The estimation uses a Gauss-Newton method with MAP
% point estimators at each iteration.  The covariance of e is a maximum
% likelihood estimator based on the residuals (ReML estimators).
% This corresponds to a Gauss-Newton ascent on the conditional probabilty
% p{P|y}
%
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

% data
%---------------------------------------------------------------------------
y     = Y.y;
[v m] = size(y);
ns    = v;

% covariance components
%---------------------------------------------------------------------------
if isfield(Y,'Ce')
	Q = Y.Ce;
else
	Q = spm_Ce(v*ones(1,m));
end
nh    = length(Q);
ne    = length(Q{1});

% initialise hyperparameters and precision
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
np     = size(Vp,2);
nu     = size(Ju,2);
ip     = [1:np];
iu     = [1:nu] + np;

% treat confounds as fixed effects
%---------------------------------------------------------------------------
uE     = inv(Ju'*Ju)*Ju'*y;
uC     = speye(nu,nu)*1e8;


% augment priors with confounds
%---------------------------------------------------------------------------
pE     = [Vp'*pE; uE];
pC     = blkdiag(Vp'*pC*Vp, uC);
ipC    = inv(pC);

% EM 
%===========================================================================
p      = pE;
dv     = 1e-6;
lm     = 0;
F      = -Inf;
dFdh   = sparse(nh,1);
dFdhh  = sparse(nh,nh);
for  k = 1:32

        % integrate
        %-------------------------------------------------------------------
        fp     = spm_int(Vp*p(ip),M,U,v);

	% compute partial derivatives [Jp] df(p)/dp*Vp {p = parameters}
	%-------------------------------------------------------------------
        vp     = Vp*p(ip);
	for  i = ip
		pi      = vp + Vp(:,i)*dv;
		dfp     = spm_int(pi,M,U,v) - fp;
		Jp(:,i) = dfp(:)/dv;
	end

        % E and dEdp = J
	%-------------------------------------------------------------------
        Eq     = y(:) - fp(:) - Ju*p(iu);
        Ep     = p - pE;
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
                dFdh(i,1)      =  trace(PS{i})/2 - Eq'*P{i}*Eq/2 ...
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
        Fp    = - Eq'*iS*Eq/2 ...
                - Ep'*ipC*Ep/2 ...
                - log(det(S))/2 ...
                + log(det(Cp))/2;

    	% Levenburg Marquardt regulatisation
	%-------------------------------------------------------------------
        if Fp < F
                if ~lm
                    lm = speye(np + nu)/128;
                else
                    lm = lm*2;
                end
	        p  = p - dp;
        else
	        % update F
	        %----------------------------------------------------------
                F  = Fp;
        end

	% E-Step: Conditional estimator of new expansion point E{p|y}
	%===================================================================
        dFdp  = -J'*iS*Eq - ipC*Ep;
        dFdpp = -J'*iS*J  - ipC;
        dp    = inv(-dFdpp + lm*norm(full(dFdpp)))*dFdp;
        
	% update - ensuring the system is dissipative (and break if not)
	%-------------------------------------------------------------------
	for i = 1:8
		A    = spm_bi_reduce(M,Vp*(p(ip) + dp(ip)));
		s    = max(real(eig(full(A))));
		if s > 0
			dp = dp/2;
		else
			break
		end
	end
        if i == 8
                warndlg('system is unstable')
                break
        end
        p     = p + dp;

	% convergence
	%-------------------------------------------------------------------
	nmp   = dp(ip)'*dp(ip);
        fprintf('%-6s: %i %6s %e %6s %e\n','E-Step',k,'F:',F,'dp:',nmp)
	if nmp < 1e-5, break, end

	% graphics
	%-------------------------------------------------------------------
	if length(dbstack) < 3
		bar(Vp*p(ip))
		xlabel('parameter')
		ylabel('conditional expectation')
		title(sprintf('%s: %i','E-Step',j))
		grid on
		drawnow
	end

end

% outputs
%---------------------------------------------------------------------------
Ep     = Vp*p(ip);
Cp     = Vp*Cp(ip,ip)*Vp';










