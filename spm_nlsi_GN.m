function [Ep,Cp,S,F] = spm_nlsi_GN(M,U,Y)
% Bayesian Parameter estimation using a Levenburg-Marquardt/EM algorithm
% FORMAT [Ep,Cp,Ce,F] = spm_nlsi_GN(M,U,Y)
%
% Dynamical MIMO models
%___________________________________________________________________________
% M.f  - f(x,u,P)
% M.g  - g(x,u,P)
%	x  - state  variables
%	u  - inputs or causes
%	P  - free parameters
% M.pE - prior expectation  - E{P}
% M.pC - prior covariance   - Cov{P}
% M.IS - integration scheme (function name or handle f(P,M,U,varargin))
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
% log evidence
%---------------------------------------------------------------------------
% F       - [-ve] free energy F = log evidence = p(y|f,g,pE,pC) = p(y|m)
%
%___________________________________________________________________________
% Returns the moments of the posterior p.d.f. of the parameters of a 
% dynamic MIMO input-state-ouput model under Gaussian assumptions
%
%              dx/dt = f(x,u,P)
%              y     = g(x,u,P) + e
%
% A static nonlinear observation model with fixed input or causes u
% obtains when x = []. i.e.
%
%              y     = g(u,P) + e
%
% Priors on the free parameters P are specified in terms of expectation pE
% and covariance pC. The E-Step uses a Fisher-Scoring scheme and a Laplace
% approximation to estimate the conditional expecation and covariance of P
% If the free-energy starts to increase, a Levenburg-Marquardt scheme is
% invoked.  The M-Step estimates the covariance components of e, in terms
% of [Re]ML point estimators.
%
%---------------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_nlsi_GN.m 112 2005-05-04 18:20:52Z john $



% integration scheme
%---------------------------------------------------------------------------
try
    f = fcnchk(M.IS);
catch
    f = 'spm_int';
end

% data y
%---------------------------------------------------------------------------
y       = Y.y;
[ns nr] = size(y);			% number of samples and responses

% covariance components
%---------------------------------------------------------------------------
try
    Q = Y.Ce;
catch
    Q = spm_Ce(ns*ones(1,nr));
end
nh    = length(Q);          % number of variance components
ne    = length(Q{1});       % number of error terms

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

% prior moments
%---------------------------------------------------------------------------
pE    = M.pE;
pC    = M.pC;

% confounds (if specified)
%---------------------------------------------------------------------------
try
    Ju = kron(speye(nr,nr),Y.X0);
catch
    Ju = [];
end

% dimension reduction of parameter space
%---------------------------------------------------------------------------
V     = spm_svd(pC,1e-8);
nu    = size(Ju,2);		              % number of parameters (confounds)
np    = size(V, 2);		              % number of parameters (effective)
ip    = [1:np];
iu    = [1:nu] + np;

% second-order moments (in reduced space)
%---------------------------------------------------------------------------
pC    = V'*pC*V;
uC    = speye(nu)*1e+8;
ipC   = inv(blkdiag(pC, uC));

% intialise
%---------------------------------------------------------------------------
p     = sparse(np + nu,1);
Ep    = pE;
Cp    = pC;

C.F   = -Inf;
dv    = 1/128;
lm    = 0;
dFdh  = sparse(nh,1);
dFdhh = sparse(nh,nh);


% EM
%===========================================================================
for k = 1:64
    
    
    % M-Step: ReML estimator of variance components:  h = max{F(p,h)}
    %===================================================================
    
    % prediction and errors
    %-------------------------------------------------------------------
    g     = feval(f,Ep,M,U,ns);
    e     = y(:) - g(:) - Ju*p(iu);
    
    % compute partial derivatives [Jp] df(p)/dp*V {p = parameters}
    %-------------------------------------------------------------------
    for i = ip
        dV      = dv*sqrt(Cp(i,i));
        dg      = feval(f,Ep + V(:,i)*dV,M,U,ns) - g;
        Jp(:,i) = dg(:)/dV;
    end
    J     = -[Jp Ju];

    % iterate a Fisher scoring scheme to find h = max{F(p,h)}
    %-------------------------------------------------------------------
    for m = 1:64
        
        % precision and conditional covariance
        %---------------------------------------------------------------
        S     = sparse(ne,ne);
        for i = 1:nh
            S = S + h(i)*Q{i};
        end
        iS    = inv(S);
        Cp    = inv(J'*iS*J  + ipC);
        
        % precision operators for M-Step
        %---------------------------------------------------------------
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
        nmh   = norm(dh)/norm(h);
        if nmh < 1e-5, break, end

    end
            
    % E-Step with Levenburg-Marquardt regularisation
    %===================================================================
    
    % objective function: F(p) (= log evidence - divergence)
    %-------------------------------------------------------------------
    F     = - e'*iS*e/2 ...
            - p'*ipC*p/2 ...
            - ne*log(8*atan(1))/2 ...
            + spm_logdet(iS )/2 ...
            + spm_logdet(ipC)/2 ...
            + spm_logdet(Cp )/2;
    
    
    % if F has increased, update gradients and curvatures for E-Step
    %-------------------------------------------------------------------
    if F > C.F
        
        % accept current estimates
        %---------------------------------------------------------------
        C.p   = p;
        C.h   = h;
        C.F   = F;

        % E-Step: Conditional update of gradients and curvature
        %---------------------------------------------------------------
        dFdp  = -J'*iS*e - ipC*p;
        dFdpp = -J'*iS*J - ipC;

        % decrease regularization
        %---------------------------------------------------------------
        lm    = lm/2;
        str   = 'E-Step(-)';
        
        % graphics
        %---------------------------------------------------------------
        if length(dbstack) < 4
            
            % subplot parameters
            %-----------------------------------------------------------
            subplot(2,1,1)
            plot([1:ns]*Y.dt,g),                        hold on
            plot([1:ns]*Y.dt,g + reshape(e,ns,nr),':'), hold off
            xlabel('time')
            title(sprintf('%s: %i','E-Step',k))
            grid on
            drawnow
            
            % subplot parameters
            %-----------------------------------------------------------
            subplot(2,1,2)
            bar(V*p(ip))
            xlabel('parameter')
            title('conditional [minus prior] expectation')
            grid on
            drawnow
        end
        
    else
        
        % reset expansion point
        %---------------------------------------------------------------
        p     = C.p;
        h     = C.h;
    
        % and increase regularization
        %---------------------------------------------------------------
        lm    = max(lm*4,1/512);
        str   = 'E-Step(+)';
        
    end
    
    % E-Step: update
    %==================================================================
    l     = lm*norm(full(dFdpp))*speye(np + nu);
    dp    = inv(-dFdpp + l)*dFdp;
    p     = p  + dp;
    Ep    = pE + V*p(ip);
    
    % convergence
    %-------------------------------------------------------------------
    nmp    = dp'*ipC*dp;
    fprintf('%-6s: %i %6s %e %6s %e\n',str,k,'F:',C.F,'dp:',full(nmp))
    if nmp < 1e-5, break, end
    
end

% outputs
%-----------------------------------------------------------------------
Ep     = pE + V*p(ip);
Cp     = V*Cp(ip,ip)*V';
F      = C.F;

return
