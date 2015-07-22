function [DEM] = spm_DEM(DEM)
% Dynamic expectation maxmisation (Variational Laplacian filtering)
% FORMAT DEM   = spm_DEM(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - response variable, output or data
% DEM.U  - explanatory variables, inputs or prior expectation of causes
% DEM.X  - confounds
%__________________________________________________________________________
%
% generative model
%--------------------------------------------------------------------------
%   M(i).g  = y(t)  = g(x,v,P)    {inline function, string or m-file}
%   M(i).f  = dx/dt = f(x,v,P)    {inline function, string or m-file}
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h log-precision (cause noise)
%   M(i).hC = prior covariances of h log-precision (cause noise)
%   M(i).gE = prior expectation of g log-precision (state noise)
%   M(i).gC = prior covariances of g log-precision (state noise)
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%   M(i).xP = precision (states)
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%
% conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.x    = Conditional expectation of hidden states
%   qU.v    = Conditional expectation of causal states
%   qU.w    = Conditional prediction error (states)
%   qU.z    = Conditional prediction error (causes)
%   qU.C    = Conditional covariance: cov(v)
%   qU.S    = Conditional covariance: cov(x)
%
% conditional moments of model-parameters - q(p)
%--------------------------------------------------------------------------
%   qP.P    = Conditional expectation
%   qP.C    = Conditional covariance
%
% conditional moments of hyper-parameters (log-transformed) - q(h)
%--------------------------------------------------------------------------
%   qH.h    = Conditional expectation (cause noise)
%   qH.g    = Conditional expectation (state noise)
%   qH.C    = Conditional covariance
%
% F         = log evidence = log marginal likelihood = negative free energy
%__________________________________________________________________________
%
% spm_DEM implements a variational Bayes (VB) scheme under the Laplace
% approximation to the conditional densities of states (u), parameters (p)
% and hyperparameters (h) of any analytic nonlinear hierarchical dynamic
% model, with additive Gaussian innovations.  It comprises three
% variational steps (D,E and M) that update the conditional moments of u, p
% and h respectively
%
%                D: qu.u = max <L>q(p,h)
%                E: qp.p = max <L>q(u,h)
%                M: qh.h = max <L>q(u,p)
%
% where qu.u corresponds to the conditional expectation of hidden states x
% and causal states v and so on.  L is the ln p(y,u,p,h|M) under the model
% M. The conditional covariances obtain analytically from the curvature of
% L with respect to u, p and h.
%
% The D-step is embedded in the E-step because q(u) changes with each
% sequential observation.  The dynamical model is transformed into a static
% model using temporal derivatives at each time point.  Continuity of the
% conditional trajectories q(u,t) is assured by a continuous ascent of F(t)
% in generalised co-ordinates.  This means DEM can deconvolve online and
% can represents an alternative to Kalman filtering or alternative Bayesian
% update procedures.
%
%
% To accelerate computations one can specify the nature of the model using
% the field:
%
% M(1).E.linear = 0: full        - evaluates 1st and 2nd derivatives
% M(1).E.linear = 1: linear      - equations are linear in x and v
% M(1).E.linear = 2: bilinear    - equations are linear in x, v and x.v
% M(1).E.linear = 3: nonlinear   - equations are linear in x, v, x.v, and x.x
% M(1).E.linear = 4: full linear - evaluates 1st derivatives (for generalised 
%                                  filtering, where parameters change)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM.m 6502 2015-07-22 11:37:13Z karl $


% check model, data, priors and confounds and unpack
%--------------------------------------------------------------------------
[M,Y,U,X] = spm_DEM_set(DEM);
 
% find or create a DEM figure
%--------------------------------------------------------------------------
Fdem = spm_figure('GetWin','DEM');
 
% tolerance for changes in norm
%--------------------------------------------------------------------------
TOL  = exp(-4);
 
% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d + 1;                   % embedding order of q(v)
n    = M(1).E.n + 1;                   % embedding order of q(x) (n >= d)
 
% number of states and parameters
%--------------------------------------------------------------------------
nY   = size(Y,2);                      % number of samples
nl   = size(M,2);                      % number of levels
nv   = sum(spm_vec(M.m));              % number of v (casual states)
nx   = sum(spm_vec(M.n));              % number of x (hidden states)
ny   = M(1).l;                         % number of y (inputs)
nc   = M(end).l;                       % number of c (prior causes)
nu   = nv*d + nx*n;                    % number of generalised states
 
% number of iterations
%--------------------------------------------------------------------------
try, nD = M(1).E.nD; catch, nD = 1; end
try, nE = M(1).E.nE; catch, nE = 8; end
try, nM = M(1).E.nM; catch, nM = 8; end
try, K  = M(1).E.K;  catch, K  = 1; end
 
 
% initialise regularisation parameters
%--------------------------------------------------------------------------
if nx
    td = 1/nD;                         % integration time for D-Step
else
    td = {2};
end
if M(1).E.linear == 1
    te = 4;                            % integration time for E-Step
else
    te = 0;
end
tm     = 4;                             % integration time for M-Step
 
% precision components Q{} requiring [Re]ML estimators (M-Step)
%==========================================================================
Q     = {};
for i = 1:nl
    v0{i,i} = sparse(M(i).l,M(i).l);
    w0{i,i} = sparse(M(i).n,M(i).n);
end
V0    = kron(sparse(n,n),spm_cat(v0));
W0    = kron(sparse(n,n),spm_cat(w0));
Qp    = blkdiag(V0,W0);
for i = 1:nl
    
    % precision (R) and covariance of generalised errors
    %----------------------------------------------------------------------
    iVv   = spm_DEM_R(n,M(i).sv);
    iVw   = spm_DEM_R(n,M(i).sw);
    
    % noise on causal states (Q)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).Q)
        q          = v0;
        q{i,i}     = M(i).Q{j};
        Q{end + 1} = blkdiag(kron(iVv,spm_cat(q)),W0);
    end
    
    % and fixed components (V)
    %----------------------------------------------------------------------
    q      = v0;
    q{i,i} = M(i).V;
    Qp     = Qp + blkdiag(kron(iVv,spm_cat(q)),W0);
    
    % noise on hidden states (R)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).R)
        q          = w0;
        q{i,i}     = M(i).R{j};
        Q{end + 1} = blkdiag(V0,kron(iVw,spm_cat(q)));
    end
    
    % and fixed components (W)
    %----------------------------------------------------------------------
    q      = w0;
    q{i,i} = M(i).W;
    Qp     = Qp + blkdiag(V0,kron(iVw,spm_cat(q)));
    
end
 
 
% number of hyperparameters
%--------------------------------------------------------------------------
nh    = length(Q);
 
% fixed priors on states (u)
%--------------------------------------------------------------------------
xP    = spm_cat(spm_diag({M.xP}));
Px    = kron(spm_DEM_R(n,0),xP);
Pv    = kron(spm_DEM_R(d,0),sparse(nv,nv));
Pu    = spm_cat(spm_diag({Px Pv}));
Pu    = Pu + speye(nu,nu)*nu*eps;
 
% hyperpriors
%--------------------------------------------------------------------------
ph.h  = spm_vec({M.hE M.gE});               % prior expectation of h
ph.c  = spm_cat(spm_diag({M.hC M.gC}));     % prior covariances of h
qh.h  = ph.h;                               % conditional expectation
qh.c  = ph.c;                               % conditional covariance
ph.ic = spm_pinv(ph.c);                     % prior precision
 
% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
    
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC,0);         % basis for parameters
    M(i).p    = size(qp.u{i},2);            % number of qp.p
    qp.p{i}   = sparse(M(i).p,1);           % initial qp.p
    pp.c{i,i} = qp.u{i}'*M(i).pC*qp.u{i};   % prior covariance
    
end
Up    = spm_cat(spm_diag(qp.u));
 
% initialise and augment with confound parameters B; with flat priors
%--------------------------------------------------------------------------
np    = sum(spm_vec(M.p));                  % number of model parameters
nb    = size(X,1);                          % number of confounds
nn    = nb*ny;                              % number of nuisance parameters
nf    = np + nn;                            % number of free parameters
ip    = (1:np);
ib    = (1:nn) + np;
pp.c  = spm_cat(pp.c);
pp.ic = spm_inv(pp.c);
pp.p  = spm_vec(qp.p);
 
% initialise conditional density q(p) := qp.e (for D-Step)
%--------------------------------------------------------------------------
for i = 1:(nl - 1)
    try
        qp.e{i} = qp.p{i} + qp.u{i}'*(spm_vec(M(i).P) - spm_vec(M(i).pE));
    catch
        qp.e{i} = qp.p{i};
    end
end
qp.e  = spm_vec(qp.e);
qp.c  = sparse(nf,nf);
qp.b  = sparse(ny,nb);
 
 
% initialise dedb
%--------------------------------------------------------------------------
for i = 1:nl
    dedbi{i,1} = sparse(M(i).l,nn);
end
for i = 1:nl - 1
    dndbi{i,1} = sparse(M(i).n,nn);
end
for i = 1:n
    dEdb{i,1}  = spm_cat(dedbi);
end
for i = 1:n
    dNdb{i,1}  = spm_cat(dndbi);
end
dEdb  = [dEdb; dNdb];
 
 
% initialise cell arrays for D-Step; e{i + 1} = (d/dt)^i[e] = e[i]
%==========================================================================
qu.x      = cell(n,1);
qu.v      = cell(n,1);
qu.y      = cell(n,1);
qu.u      = cell(n,1);
[qu.x{:}] = deal(sparse(nx,1));
[qu.v{:}] = deal(sparse(nv,1));
[qu.y{:}] = deal(sparse(ny,1));
[qu.u{:}] = deal(sparse(nc,1));
 
% initialise cell arrays for hierarchical structure of x[0] and v[0]
%--------------------------------------------------------------------------
x         = {M(1:end - 1).x};
v         = {M(1 + 1:end).v};
qu.x{1}   = spm_vec(x);
qu.v{1}   = spm_vec(v);
 
% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dx    = kron(spm_speye(n,n,1),spm_speye(nx,nx,0));
Dv    = kron(spm_speye(d,d,1),spm_speye(nv,nv,0));
Dy    = kron(spm_speye(n,n,1),spm_speye(ny,ny,0));
Dc    = kron(spm_speye(d,d,1),spm_speye(nc,nc,0));
D     = spm_cat(spm_diag({Dx,Dv,Dy,Dc}));
 
% and null blocks
%--------------------------------------------------------------------------
dVdy  = sparse(n*ny,1);
dVdc  = sparse(d*nc,1);
dVdyy = sparse(n*ny,n*ny);
dVdcc = sparse(d*nc,d*nc);
 
% gradients and curvatures for conditional uncertainty
%--------------------------------------------------------------------------
dWdu  = sparse(nu,1);
dWdp  = sparse(nf,1);
dWduu = sparse(nu,nu);
dWdpp = sparse(nf,nf);
 
% preclude unnecessary iterations
%--------------------------------------------------------------------------
if ~nh,        nM = 1; end
if ~nf && ~nh, nE = 1; end
 
 
% preclude very precise states from entering free-energy/action
%--------------------------------------------------------------------------
ix     = (1:(nx*n)) + ny*n + nv*n;
iv     = (1:(nv*d)) + ny*n;
je     = diag(Qp) < exp(16);
ju     = [je(ix); je(iv)];
 
 
% E-Step: (with embedded D and M-Steps)
%==========================================================================
Fi     = -Inf;
for iE = 1:nE
    
    % get time and clear persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic;  clear spm_DEM_eval
    
    % [re-]set accumulators for E-Step
    %----------------------------------------------------------------------
    dFdh  = sparse(nh,1);                 % gradient   (hyperparamteres)
    dFdhh = sparse(nh,nh);                % curvatiure (hyperparamteres)
    dFdp  = sparse(nf,1);                 % gradient   (paramteres)
    dFdpp = sparse(nf,nf);                % curvatiure (paramteres)
    qp.ic = sparse(0);                    % conditional precision (p)
    iqu.c = sparse(0);                    % conditional information (p)
    EE    = sparse(0);
    ECE   = sparse(0);
    
    % [re-]set precisions using ReML hyperparameter estimates
    %----------------------------------------------------------------------
    iS    = Qp;
    for i = 1:nh
        iS = iS + Q{i}*exp(qh.h(i));
    end
    
    % [re-]adjust for confounds
    %----------------------------------------------------------------------
    Y     = Y - qp.b*X;
    
    % [re-]set states & their derivatives
    %----------------------------------------------------------------------
    try, qu = qU(1); end
    
    % D-Step: (nD D-Steps for each sample)
    %======================================================================
    for iY = 1:nY
        
        % [re-]set states for static systems
        %------------------------------------------------------------------
        if ~nx, try, qu = qU(iY); end, end
        
        % D-Step: until convergence for static systems
        %==================================================================
        Fd     = -exp(64);
        for iD = 1:nD
            
            % sampling time
            %--------------------------------------------------------------
            ts = iY + (iD - 1)/nD;
            
            % derivatives of responses and inputs
            %--------------------------------------------------------------
            try
                qu.y(1:n) = spm_DEM_embed(Y,n,ts,1,M(1).delays);
                qu.u(1:d) = spm_DEM_embed(U,d,ts);
            catch
                qu.y(1:n) = spm_DEM_embed(Y,n,ts);
                qu.u(1:d) = spm_DEM_embed(U,d,ts);
            end
            
            % compute dEdb (derivatives of confounds)
            %--------------------------------------------------------------
            b     = spm_DEM_embed(X,n,ts);
            for i = 1:n
                dedbi{1}  = -kron(b{i}',speye(ny,ny));
                dEdb{i,1} =  spm_cat(dedbi);
            end
            
            % evaluate functions:
            % E = v - g(x,v) and derivatives dE.dx, ...
            %==============================================================
            [E,dE] = spm_DEM_eval(M,qu,qp);
            
            % conditional covariance [of states {u}]
            %--------------------------------------------------------------
            qu.p   = real(dE.du'*iS*dE.du) + Pu;
            qu.c   = diag(ju)*spm_inv(qu.p)*diag(ju);
            iqu.c  = iqu.c + spm_logdet(qu.c);
            
            % and conditional covariance [of parameters {P}]
            %--------------------------------------------------------------
            dE.dP  = spm_cat({dE.dp dEdb});
            ECEu   = dE.du*qu.c*dE.du';
            ECEp   = dE.dP*qp.c*dE.dP';
            
            if ~nx
                
                % Evaluate objective function L(t) (for static models)
                %----------------------------------------------------------
                L = - trace(real(E'*iS*E))/2 ...      % states (u)
                    - trace(real(iS*ECEp))/2;         % expectation q(p)
                
                % if F is increasing, save expansion point
                %----------------------------------------------------------
                if L > Fd
                    td     = {min(td{1} + 1, 4)};
                    Fd     = L;
                    B.qu   = qu;
                    B.E    = E;
                    B.dE   = dE;
                    B.ECEp = ECEp;
                    
                else
                    
                    % otherwise, return to previous expansion point
                    %------------------------------------------------------
                    qu     = B.qu;
                    E      = B.E;
                    dE     = B.dE;
                    ECEp   = B.ECEp;
                    td     = {min(td{1} - 2,-4)};
                end
            end
            
            % save states at qu(t)
            %--------------------------------------------------------------
            if iD == 1
                qE{iY} = E;
                qU(iY) = qu;
            end
            
            
            % uncertainty about parameters dWdv, ... ; W = ln(|qp.c|)
            %==============================================================
            if np
                for i = 1:nu
                    CJp(:,i)   = spm_vec(qp.c(ip,ip)*dE.dpu{i}'*iS);
                    dEdpu(:,i) = spm_vec(dE.dpu{i}');
                end
                dWdu   = real(CJp'*spm_vec(dE.dp'));
                dWduu  = real(CJp'*dEdpu);
            end
            
            
            % D-step update: of causes v{i}, and hidden states x(i)
            %==============================================================
            
            % conditional modes
            %--------------------------------------------------------------
            q     = {qu.x{1:n} qu.v{1:d} qu.y{1:n} qu.u{1:d}};
            u     = spm_vec(q);
            
            % first-order derivatives
            %--------------------------------------------------------------
            dVdu  = -real(dE.du'*iS*E)     - dWdu/2  - Pu*u(1:nu);
            
            % and second-order derivatives
            %--------------------------------------------------------------
            dVduu = -real(dE.du'*iS*dE.du) - dWduu/2 - Pu;
            dVduy = -real(dE.du'*iS*dE.dy);
            dVduc = -real(dE.du'*iS*dE.dc);
            
            % gradient
            %--------------------------------------------------------------
            dFdu  = spm_vec({dVdu;  dVdy;  dVdc });
            
            % Jacobian (variational flow)
            %--------------------------------------------------------------
            dFduu = spm_cat({dVduu  dVduy  dVduc  ;
                             []     dVdyy  []     ;
                             []     []     dVdcc});
            
            
            % update conditional modes of states
            %==============================================================
            f     = K*dFdu  + D*u;
            dfdu  = K*dFduu + D;
            
            du    = spm_dx(dfdu,f,td);
            q     = spm_unvec(u + du,q);
                        
            % and save them
            %--------------------------------------------------------------
            qu.x(1:n) = q((1:n));
            qu.v(1:d) = q((1:d) + n);
            
            % save Lyapunov exponents (eigenvalues) if requested
            %--------------------------------------------------------------
            if iD == 1 && isfield(DEM,'E')
                DEM.E(:,iY) = eig(full(dfdu));
            end
                        
            % D-Step: break if convergence (for static models)
            %--------------------------------------------------------------
            if ~nx
                qU(iY) = qu;
            end
            if ~nx && ((dFdu'*du < TOL) || (norm(du,1) < TOL))
                break
            end
            
        end % D-Step
        
        % Gradients and curvatures for E-Step: W = tr(C*J'*iS*J)
        %==================================================================
        if np
            for i = ip
                CJu(:,i)   = spm_vec(qu.c*dE.dup{i}'*iS);
                dEdup(:,i) = spm_vec(dE.dup{i}');
            end
            dWdp(ip)       = CJu'*spm_vec(dE.du');
            dWdpp(ip,ip)   = CJu'*dEdup;
        end
        
        
        % Accumulate; dF/dP = <dL/dp>, dF/dpp = ...
        %------------------------------------------------------------------
        dFdp  = dFdp  - dWdp/2  - real(dE.dP'*iS*E);
        dFdpp = dFdpp - dWdpp/2 - real(dE.dP'*iS*dE.dP);
        qp.ic = qp.ic           + real(dE.dP'*iS*dE.dP);
        
        % and quantities for M-Step
        %------------------------------------------------------------------
        EE    = real(E*E') + EE;
        ECE   = ECE + ECEu + ECEp;

        
    end % sequence (nY)
    
    
    % M-step - optimise hyperparameters (mh = total update)
    %======================================================================
    mh     = 0;                      
    for iM = 1:nM
        
        % [re-]set precisions using ReML hyperparameter estimates
        %------------------------------------------------------------------
        iS    = Qp;
        for i = 1:nh
            iS = iS + Q{i}*exp(qh.h(i));
        end
        S     = spm_inv(iS);
        dS    = ECE + EE - S*nY;
        
        % 1st-order derivatives: dFdh = dF/dh
        %------------------------------------------------------------------
        for i = 1:nh
            dPdh{i}   =    Q{i}*exp(qh.h(i));
            dFdh(i,1) = -trace(dPdh{i}*dS)/2;
        end
        
        % 2nd-order derivatives: dFdhh
        %------------------------------------------------------------------
        for i = 1:nh
            for j = 1:nh
                dFdhh(i,j) = -trace(dPdh{i}*S*dPdh{j}*S*nY)/2;
            end
        end
        
        % hyperpriors
        %------------------------------------------------------------------
        qh.e  = qh.h  - ph.h;
        dFdh  = dFdh  - ph.ic*qh.e;
        dFdhh = dFdhh - ph.ic;
        
        % update ReML estimate of parameters
        %------------------------------------------------------------------
        dh    = spm_dx(dFdhh,dFdh,{tm});
        dh    = max(min(dh,2),-2);
        qh.h  = qh.h + dh;
        mh    = mh   + dh;
        
        % conditional covariance of hyperparameters
        %------------------------------------------------------------------
        qh.c  = -spm_inv(dFdhh);
        
        % convergence (M-Step)
        %------------------------------------------------------------------
        if (dFdh'*dh < TOL) || (norm(dh,1) < TOL), break, end
        
    end % M-Step
    
 
    % conditional precision of parameters
    %------------------------------------------------------------------
    qp.ic(ip,ip) = qp.ic(ip,ip) + pp.ic;
    qp.c         = spm_inv(qp.ic);
    
    % evaluate objective function (F)
    %======================================================================
    
    % free-energy and action
    %----------------------------------------------------------------------
    Lu  = - trace(iS(je,je)*EE(je,je))/2 ...          % states (u)
          - n*ny*log(2*pi)*nY/2 ...                   % constant
          + spm_logdet(iS(je,je))*nY/2 ...            % entropy - error
          + iqu.c/(2*nD);                             % entropy q(u)
    
    Lp  = - trace(qp.e'*pp.ic*qp.e)/2  ...            % parameters (p)
          - trace(qh.e'*ph.ic*qh.e)/2  ...            % hyperparameters (h)
          + spm_logdet(qp.c(ip,ip)*pp.ic)/2   ...     % entropy q(p)
          + spm_logdet(qh.c*ph.ic)/2;                 % entropy q(h)
    
    La  = - trace(qp.e'*pp.ic*qp.e)*nY/2   ...        % parameters (p)
          - trace(qh.e'*ph.ic*qh.e)*nY/2   ...        % hyperparameters (h)
          + spm_logdet(qp.c(ip,ip)*pp.ic*nY)*nY/2 ... % entropy q(p)
          + spm_logdet(qh.c*ph.ic*nY)*nY/2;           % entropy q(h)
    
    Li  = Lu + Lp;                                    % free-energy
    Ai  = Lu + La;                                    % free-action
    
    
    % if F is increasing, save expansion point and derivatives
    %------------------------------------------------------------------
    if Li > Fi || iE < 2
        
        
        % Accept free-energy and save current parameter estimates
        %------------------------------------------------------------------
        Fi      = Li;
        te      = min(te + 1/2,4);
        tm      = min(tm + 1/2,4);
        B.qp    = qp;
        B.qh    = qh;
        B.pp    = pp;
        
        % E-step: update expectation (p)
        %==================================================================
        
        % gradients and curvatures
        %------------------------------------------------------------------
        dFdp(ip)     = dFdp(ip)     - pp.ic*(qp.e - pp.p);
        dFdpp(ip,ip) = dFdpp(ip,ip) - pp.ic;
        
        % update conditional expectation
        %------------------------------------------------------------------
        dp      = spm_dx(dFdpp,dFdp,{te});
        qp.e    = qp.e + dp(ip);
        qp.p    = spm_unvec(qp.e,qp.p);
        qp.b    = spm_unvec(dp(ib),qp.b);

    else
        
        % otherwise, return to previous expansion point
        %------------------------------------------------------------------
        nM      = 1;
        qp      = B.qp;
        pp      = B.pp;
        qh      = B.qh;
        te      = min(te - 2, -2);
        tm      = min(tm - 2, -2); 
        
    end
    F(iE)   = Fi;
    A(iE)   = Ai;
        
    
    % save model-states (for each time point)
    %==================================================================
    for t = 1:length(qU)
        v     = spm_unvec(qU(t).v{1},v);
        x     = spm_unvec(qU(t).x{1},x);
        z     = spm_unvec(qE{t}(1:(ny + nv)),{M.v});
        w     = spm_unvec(qE{t}([1:nx] + (ny + nv)*n),{M.x});
        for i = 1:(nl - 1)
            if M(i).m, QU.v{i + 1}(:,t) = spm_vec(v{i}); end
            if M(i).n, QU.x{i}(:,t)     = spm_vec(x{i}); end
            if M(i).n, QU.w{i}(:,t)     = spm_vec(w{i}); end
            if M(i).l, QU.z{i}(:,t)     = spm_vec(z{i}); end
        end
        QU.v{1}(:,t)  = spm_vec(qU(t).y{1}) - spm_vec(z{1});
        if M(nl).l, QU.z{nl}(:,t) = spm_vec(z{nl});      end
        
        % and conditional covariances
        %--------------------------------------------------------------
        i       = (1:nx);
        QU.S{t} = qU(t).c(i,i);
        i       = (1:nv) + nx*n;
        QU.C{t} = qU(t).c(i,i);
    end
    
    % report and break if convergence
    %------------------------------------------------------------------
    spm_figure('Select', Fdem)
    spm_DEM_qU(QU)
    if np
        subplot(2*nl,2,4*nl)
        bar(full(Up*qp.e))
        xlabel({'parameters {minus prior}'})
    end
    if nh
        subplot(2*nl,4,8*nl - 4)
        bar(full(qh.h))
        title({'log-precision'})
    end
    if length(F) > 2
        subplot(2*nl,4,8*nl - 5)
        plot(F - F(1))
        xlabel('updates')
        title('free-energy')
    end
    drawnow
    
    % report (EM-Steps)
    %------------------------------------------------------------------
    str{1} = sprintf('DEM: %i (%i:%i)',iE,iD,iM);
    str{2} = sprintf('F:%.4e',full(F(iE) - F(1)));
    str{3} = sprintf('p:%.2e',full(norm(dp,1)));
    str{4} = sprintf('h:%.2e',full(norm(mh,1)));
    str{5} = sprintf('(%.2e sec)',full(toc));
    fprintf('%-16s%-16s%-14s%-14s%-16s\n',str{:})
    
    
    % Convergence
    %------------------------------------------------------------------
    if (norm(dp,1) < TOL*norm(spm_vec(qp.p),1)) && (norm(mh,1) < TOL), break, end
    if te < -8, break, end
    
end
spm_figure('Focus', Fdem)

% Assemble output arguments
%==========================================================================
 
% conditional moments of model-parameters (rotated into original space)
%--------------------------------------------------------------------------
qP.P     = spm_unvec(Up*qp.e + spm_vec(M.pE),M.pE);
qP.C     = Up*qp.c(ip,ip)*Up';
qP.V     = spm_unvec(diag(qP.C),M.pE);
qP.dFdp  = Up*dFdp(ip);
qP.dFdpp = Up*dFdpp(ip,ip)*Up';
 
% conditional moments of hyper-parameters (log-transformed)
%--------------------------------------------------------------------------
qH.h   = spm_unvec(qh.h,{{M.hE} {M.gE}});
qH.g   = qH.h{2};
qH.h   = qH.h{1};
qH.C   = qh.c;
qH.V   = spm_unvec(diag(qH.C),{{M.hE} {M.gE}});
qH.W   = qH.V{2};
qH.V   = qH.V{1};
 
% assign output variables
%--------------------------------------------------------------------------
DEM.M  = M;
DEM.U  = U;                   % causes
DEM.X  = X;                   % confounds
 
DEM.qU = QU;                  % conditional moments of model-states
DEM.qP = qP;                  % conditional moments of model-parameters
DEM.qH = qH;                  % conditional moments of hyper-parameters
 
DEM.F  = F;                   % [-ve] Free energy
DEM.S  = A;                   % [-ve] Free action
