function [QU,qP,qH,F] = spm_DEM(M,Y,U,X)
% FORMAT [qU,qP,qH,F] = spm_DEM(M,Y,U,X)
% FORMAT DEM          = spm_DEM(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - inputs or data
% DEM.U  - prior expectation of causes
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
%   M(i).hE = prior expectation of h hyper-parameters
%   M(i).hC = prior covariances of h hyper-parameters
%   M(i).Q  = precision components
%   M(i).V  = fixed covariance component
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%
%  [M(i).P  = initial model parameters]  (default: M(i).pE)]
%  [M(i).h  = initial hyper parameters]  (default: M(i).hE)]
%
%
% conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qP.x    = Conditional expectation of hidden states
%   qP.v    = Conditional expectation of causal states
%   qP.e    = Conditional residuals
%   qP.C    = Conditional covariance: cov(v)
%
% conditional moments of model-parameters - q(p)
%--------------------------------------------------------------------------
%   qP.P    = Conditional expectation
%   qP.Pi   = Conditional expectation for each level
%   qP.C    = Conditional covariance
%  
% conditional moments of hyper-parameters (log-transformed) - q(h)
%--------------------------------------------------------------------------
%   qH.h    = Conditional expectation
%   qH.hi   = Conditional expectation for each level
%   qH.C    = Conditional covariance
%   qH.iC   = Component  precision: cov(vec(e[:})) = inv(kron(iC,iV))
%   qH.iV   = Sequential precision
%
% F         = log evidence = marginal likelihood = negative free energy
%__________________________________________________________________________
%
% spm_DEM implements a variational Bayes (VB) scheme under the Laplace
% approximation to the conditional densities of states (u), parameters (p)
% and hyperparameters (h) of any analytic nonlinear hierarchical dynamic
% model, with additive Gaussian innovations.  It comprises three
% variational steps (D,E and M) that update the conditional moments of u, p
% and h respectively
%
%                D: qu.v = max <L>q(p,h)
%                E: qp.p = max <L>q(v,h)
%                M: qh.h = max <L>q(v,p)
%
% where qu.v corresponds to the conditional expectation of states v, and so
% on.  L is the ln p(y,v,p,h) under the model M.  The hidden states x are a
% deterministic function of the causes v and therefore do not enter the
% conditional density.  The conditional covariances obtain analytically
% from the curvature of L with respect to v, p and h.
%
% The D-step is embedded in the E-step because q(v) changes with each
% sequential observation.  The dynamical model is transformed into a static
% model using temporal derivatives at each time point.  Continuity of the
% conditional trajectories q(v,t) is assured by a continuous ascent of F(t)
% evaluated in the future.  This means DEM can deconvolve online and can be
% used as an alternative to Kalman filtering or alternative Bayesian update
% procedures.
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% create a graphics figure
%--------------------------------------------------------------------------
warning off
fig = get(0,'children');
for i = 1:length(fig)
    if strcmp(get(fig(i),'name'),'Dynamic Expectation Maximisation')
        Fdem = gcf;
    end
end
try
    Fdem;
catch
    Fdem = spm_figure;
    set(Fdem,'name','Dynamic Expectation Maximisation')
end
 
% unpack structures and variables
%--------------------------------------------------------------------------
if nargin == 1
    DEM    = M;
    M      = DEM.M;
    Y      = DEM.Y;
    try  U = DEM.U; end
    try  X = DEM.X; end
end
if ~exist('U','var'), U = []; end
if ~exist('X','var'), X = []; end
 
% set model indices and missing fields
%--------------------------------------------------------------------------
M    = spm_M_set(M);

% tolerance for changes in norm
%--------------------------------------------------------------------------
TOL  = 1e-3;

% number of states and parameters
%--------------------------------------------------------------------------

nt   = size(Y,2);                       % number of samples
nl   = size(M,2);                       % number of levels
ne   = sum(cat(1,M.l));                 % number of e (errors)
nv   = sum(cat(1,M.m));                 % number of v (casual states)
nx   = sum(cat(1,M.n));                 % number of x (hidden states)
ny   = M(1).l;                          % number of y (inputs)
nc   = M(end).l;                        % number of c (prior causes)
 
% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d;                       % approximation order of q(x,v)
n    = M(1).E.n;                       % embedding order      (n >= d)
r    = M(1).E.r;                       % restriction order    (r <= d)
s    = M(1).E.s;                       % smoothness - s.d. of kernel
dt   = M(1).E.dt;                      % time bins {seconds}

% number of iterations
%--------------------------------------------------------------------------
try nD = M(1).E.nD; catch nD = 8; end
try nE = M(1).E.nE; catch nE = 8; end
try nM = M(1).E.nM; catch nM = 8; end
try nI = M(1).E.nI; catch nI = 16;end

% initialise regularisation parameters
%--------------------------------------------------------------------------
if nx
    kt = exp(8);                          % rate constant for D-Step
    td = dt/nD;                           % integration time for D-Step
    te = 64;
else
    kt = 1;
    td = {64};
    te = 64;
end

% check length of time-series
%--------------------------------------------------------------------------
if nt < n
    warndlg({'Please ensure time-series is longer than embedding order'})
    error(' ')
end
 
% check prior expectations are not all zero
%--------------------------------------------------------------------------
if ~any(spm_vec({M.P}')) && ~any(any(U))
    warndlg({'Please initialise'})
    error(' ')
end
 
% check prior expectation of causes (at level n) and confounds
%--------------------------------------------------------------------------
if ~nnz(U), U = sparse(nc,nt); end
if ~nnz(X), X = sparse(0 ,nt); end
 
% transpose causes and confounds, if specified in conventional fashion
%--------------------------------------------------------------------------
if size(U,2) < nt, U = U';    end
if size(X,2) < nt, X = X';    end
 
 
% precision components Q{} requiring [Re]ML estimators (M-Step)
%==========================================================================
Q     = {};
for i = 1:nl
    P{i,i} = sparse(M(i).l,M(i).l);
end
for i = 1:nl
    for j = 1:length(M(i).Q)
        q          = P;
        q{i,i}     = M(i).Q{j};
        Q{end + 1} = spm_cat(q);
    end
end
 
% and fixed components P
%--------------------------------------------------------------------------
P     = spm_cat(diag({M.V}));
nh    = length(Q);                      % number of hyperparameters
 
% hyperpriors
%--------------------------------------------------------------------------
try
    ph.h = spm_vec(M.hE);               % prior expectation of h
    ph.c = spm_cat(diag({M.hC}));       % prior covariances of h
catch
    ph.h = sparse(nh,1);                % prior expectation of h
    ph.c = speye(nh,nh)*16;             % prior covariances of h
end
ph.ic    = inv(ph.c);                   % prior precision
qh.h     = spm_vec(M.h);                % conditional expectation
qh.c     = ph.c;                        % conditional covariance
qh.e     = qh.h - ph.h;
 
 
% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
 
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC,0,1-8);              % basis for parameters
    qp.p{i}   = spm_vec(M(i).P) - spm_vec(M(i).pE);  % initial estimate
    qp.p{i}   = qp.u{i}'*qp.p{i};                    % projection
    pp.c{i,i} = qp.u{i}'*M(i).pC*qp.u{i};            % prior covariance
    
    M(i).p    = size(qp.u{i},2);
 
end
Up    = spm_cat(diag(qp.u));
 
% initialise and augment with confound parameters B; with flat priors
%--------------------------------------------------------------------------
np    = sum(cat(1,M.p));                    % number of model parameters
nb    = size(X,1);                          % number of confounds
nn    = nb*ny;                              % number of nuisance parameters
nf    = np + nn;                            % numer of free parameters
pp.c  = spm_cat(diag({pp.c,speye(nn,nn)*1e8}));
pp.ic = inv(pp.c);
 
% initialise conditional density q(p) (for D-Step)
%--------------------------------------------------------------------------
qp.e  = spm_vec({qp.p sparse(nn,1)});
qp.c  = sparse(nf,nf);

% initialise dedb
%--------------------------------------------------------------------------
for i = 1:nl
    dedbi{i} = sparse(M(i).l,nn);
end
 
 
% initialise cell arrays for D-Step; e{i + 1} = (d/dt)^i[e] = e[i]
%==========================================================================
qu.e      = cell(n    ,1);
qu.x      = cell(n + 1,1);
qu.v      = cell(n + 1,1);
qu.y      = cell(n + 1,1);
qu.u      = cell(n + 1,1);
[qu.e{:}] = deal(sparse(ne,1));
[qu.x{:}] = deal(sparse(nx,1));
[qu.v{:}] = deal(sparse(nv,1));
[qu.y{:}] = deal(sparse(ny,1));
[qu.u{:}] = deal(sparse(nc,1));
 
% initialise cell arrays for hierarchical structure of x[0] and v[0]
%--------------------------------------------------------------------------
v         = {M(2:end).v};
x         = {M.x};
e         = {M.v};
qu.x{1}   = spm_vec(x);
qu.v{1}   = spm_vec(v);
qU        = qu;
qU.c      = [];
qU.C      = [];

% fixed-form derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
dhdx      = cell(d,1);
dhdv      = cell(d,d);
dhdy      = cell(d,n);
dhdc      = cell(d,d);
dydy      = cell(n,n);
dcdc      = cell(d,d);
[dhdx{:}] = deal(sparse(nv,nx));
[dhdv{:}] = deal(sparse(nv,nv));
[dhdy{:}] = deal(sparse(nv,ny));
[dhdc{:}] = deal(sparse(nv,nc));
[dydy{:}] = deal(sparse(ny,ny));
[dcdc{:}] = deal(sparse(nc,nc));
 
% add constant terms
%--------------------------------------------------------------------------
for i = 2:d
    dhdv{i - 1,i}  =  speye(nv,nv);
    dcdc{i - 1,i}  =  speye(nc,nc);
end
for i = 2:n
    dydy{i - 1,i}  =  speye(ny,ny);
end
dhdy  = spm_cat(dhdy,2);
dhdc  = spm_cat(dhdc,2);
dydy  = spm_cat(dydy);
dcdc  = spm_cat(dcdc);
 
% Embedding operator (D) for updating response and related matrices
%--------------------------------------------------------------------------
pt    = spm_DEM_t(n,r,s,dt);
[R,V] = spm_DEM_P(n,r,s,dt,pt);
 
% projector for conditional covariance of states - cov{v(t + pt)}
%--------------------------------------------------------------------------
j     = 1:(d - 1);
Le    = [1 pt.^j./cumprod(j)];
Le    = kron(Le,speye(nv,nv));

% gradients and curvatures for conditional uncertainty
%--------------------------------------------------------------------------
dUdv       = sparse(nv,1);
dUdp       = sparse(nf,1);
dUdvx      = sparse(nv,nx);
dUdpp      = sparse(nf,nf);
dUdvv      = cell(1,d);
[dUdvv{:}] = deal(zeros(nv,nv));

% preclude unneceassry iterations
%--------------------------------------------------------------------------
if ~nh,        nM = 1; end
if ~nf,        nE = 1; end
if ~nf && ~nh, nI = 1; end

% Iterate DEM
%==========================================================================
Fm     = -Inf;
for iI = 1:nI
 
    % E-Step: (with embedded D-Step)
    %======================================================================
    mp     = zeros(nf,1);
    Fe     = -Inf;
    for iE = 1:nE
 
        % [re-]set states & their derivatives
        %------------------------------------------------------------------
        qu     = qU(1);
        
        % [re-]set accumulators for E-Step
        %------------------------------------------------------------------
        dFdp   = zeros(nf,1);
        dFdpp  = zeros(nf,nf);
        EE     = sparse(0);
        EJCJ   = sparse(0);
        qp.ic  = sparse(0);
        qu.c   = speye(1);
 
        
        % [re-]set hierarchical parameters
        %------------------------------------------------------------------
        qp.p   = spm_unvec(qp.e,qp.p);
        qp.b   = qp.e([1:nn] + np,1);
        
        
        % [re-]set precisions using ReML hyperparameter estimates
        %------------------------------------------------------------------
        iC    = P;
        for i = 1:nh
           iC = iC + Q{i}*exp(qh.h(i));
        end
        iS    = kron(R,iC);
 
        % [re-]set precision operator for D-Step
        %------------------------------------------------------------------
        for i = 1:d
            W{i} = kron(V{i},iC);
        end
 
        % [re-]adjust for confounds
        %------------------------------------------------------------------
        Y      = Y - reshape(qp.b,ny,nb)*X;
 
        % D-Step: (with embedded S-Step)
        %==================================================================
        for t  = 1:nt
            
            % [re-]set states & their derivatives for static systems
            %--------------------------------------------------------------
            if ~nx && iE > 1
                qu = qU(t);
            end

            % derivatives of responses and inputs,
            %--------------------------------------------------------------
            qu.y(1:n) = spm_DEM_embed(Y,t,n);
            qu.u(1:d) = spm_DEM_embed(U,t,d);
 
            % compute dEdb
            %--------------------------------------------------------------
            if nb
                for i = 1:n
                    dedbi{1,1} = -kron(D(i,:)*X(:,kn)',speye(ny,ny));
                    dedb{i,1}  = spm_cat(dedbi);
                end
                dEdb  = spm_cat(dedb);
            else
                dEdb  = sparse(ne*n,0);
            end
 
 
            % D-Step: until convergence for static systems
            %==============================================================
            Fd     = -Inf;
            for iD = 1:nD
 
                % evaluate functions:
                % e = v - g(x,v), dx/dt = f(x,v) and derivatives dE.dx, ...
                %==========================================================       
                [qu df dE] = spm_DEM_diff(M,qu,qp);
                
                % and conditional covariance [of states {v}]
                %----------------------------------------------------------
                qu.C       = inv(dE.dV'*iS*dE.dV);
                qU(t)      = qu;
                
                % vectorise error
                %----------------------------------------------------------
                dE.dP      = [dE.dp dEdb];
                E          = spm_vec(qu.e);
                JCJ        = dE.dV*qu.C*dE.dV' + dE.dP*qp.c*dE.dP';
          
                % and evaluate objective function F(t) (for static models)
                %----------------------------------------------------------
                if ~nx
                    
                    L = - trace(E'*iS*E)/(2*r) ...  % states (u)
                        - trace(iS*JCJ)/(2*r) ...      % expectation
                        + spm_logdet(qu.C)/(2*r);         % entropy q(u)

                    % if F is increasing, save expansion point and dervatives
                    %----------------------------------------------------------
                    if L > Fd
                        td  = {td{1}*2};
                        Fd  = L;
                        save tempD qu df dE
                    else

                        % otherwise, return to previous expansion point
                        %------------------------------------------------------
                        load tempD
                        td  = {min(td{1}/2,16)};
                    end
                end

                % conditional uncertainty about parameters
                %==========================================================
                if np
                    for i = 1:nv
                        
                        % 1st-order derivatives: dUdv, ... ; U = ln(|qp.c|)
                        %--------------------------------------------------
                        CJ                    = qp.c(1:np,1:np)*dE.dpv{i}';
                        dUdv(i,1)             = trace(CJ*W{1}*dE.dp);

                        % 2nd-order derivatives
                        %--------------------------------------------------
                        for j = 1:nx
                            dUdvx(i,j)        = trace(CJ*W{1}*dE.dpx{j});
                        end
                        for j = 1:nv
                            for k = 1:d
                                dUdvv{k}(i,j) = trace(CJ*W{k}*dE.dpv{j});
                            end
                        end
                    end
                end

 
                % D-step update: of causes v{i}, and other states u(i)
                %==========================================================

                % compute h = dv[d]/dt, dh/dx, ...
                %----------------------------------------------------------
                dE.dh    =  dE.dv'*W{1};
                h        = -kt*(dE.dh*E     + dUdv );
                dhdx{d}  = -kt*(dE.dh*dE.dx + dUdvx);
                dhdy{d}  = -kt*(dE.dh*dE.dy        );
                dhdc{d}  = -kt*(dE.dh*dE.dc        );
 
                % & d-th derivatives
                %----------------------------------------------------------
                for i = 1:d
                    dhdv{d,i} = -kt*(dE.dv'*W{i}*dE.dv + dUdvv{i});
                end
 
                % Curvature and gradients for states u = {v x y c}
                %----------------------------------------------------------
                dFduu = spm_cat({dhdv  dhdx  dhdy  dhdc  ;
                                 df.dv df.dx []    []    ;
                                 []    []    dydy  []    ;
                                 []    []    []    dcdc});                         
                
                dudt  = {qu.v{2:d} h qu.x{2} qu.y{2:n + 1} qu.u{2:d + 1}};
                dFdu  = spm_vec(dudt);
                
                           
                % update conditional expecatiions of states u = {x,v,y,u}
                %----------------------------------------------------------
                du    = spm_dx(dFduu,dFdu,td);
                dq    = spm_unvec(du,dudt);
                for i = 1:d
                    qu.v{i} = qu.v{i} + dq{i};
                    qu.u{i} = qu.u{i} + dq{d + 1 + n + i};
                end
                for i = 1:n
                    qu.y{i} = qu.y{i} + dq{d + 1 + i};
                end
                qu.x{1}     = qu.x{1} + dq{d + 1};
                
                % D-Step: break if convergence (for static models)
                %----------------------------------------------------------
                if ~nx && ((dFdu'*du < 1e-2) | (norm(du,1) < TOL))
                    break
                end

                % report (D-Steps)
                %----------------------------------------------------------
                str{1} = 'D-Step: ';
                str{2} = sprintf('%i',t);
                str{3} = sprintf('%i',iD);
                str{4} = sprintf('F:%.6e',full(Fd));
                str{5} = sprintf('u:%.2e',full(du'*du));
                if ~nx, fprintf('%-8s%-4s%-4s%-20s%-16s\n',str{1:5}), end

            end % D-Step
            
            % Gradients and curvatures for E-Step:
            %==============================================================
            qu.c  = qu.c*qu.C;
            EE    = E*E'+ EE;
            EJCJ  = EJCJ + JCJ;
            for i = 1:np
                
                % 1st-order derivatives: U = tr(C*J'*iS*J)
                %----------------------------------------------------------
                CJ        = qu.C*dE.dVp{i}'*iS;
                dUdp(i,1) = trace(CJ*dE.dV);
 
                % 2nd-order derivatives
                %----------------------------------------------------------
                for j = 1:np
                    dUdpp(i,j) = trace(CJ*dE.dVp{j});
                end
            end
 
            % Accumulate; dF/dP = <dL/dp>, dF/dpp = ...
            %--------------------------------------------------------------
            dFdp   = dFdp  - dUdp  - dE.dP'*iS*E;
            dFdpp  = dFdpp - dUdpp - dE.dP'*iS*dE.dP;
            qp.ic  = qp.ic         + dE.dP'*iS*dE.dP;
 
            
        end % sequence
 
        % augment with priors
        %------------------------------------------------------------------
        dFdp   = dFdp/r  - pp.ic*qp.e;
        dFdpp  = dFdpp/r - pp.ic;
        qp.ic  = qp.ic/r + pp.ic;
        qp.c   = inv(qp.ic);
            
        % evaluate objective function L(t)
        %==================================================================
        L = - trace(iS*EE)/(2*r)  ...            % states (u)
            - trace(qp.e'*pp.ic*qp.e)/2  ...     % parameters (p)
            - trace(qh.e'*ph.ic*qh.e)/2  ...     % hyperparameters (h)
            - trace(iS*EJCJ)/(2*r)   ...         % expectation under q
            + spm_logdet(qu.c)/(2*r) ...         % entropy q(u)
            + spm_logdet(qp.c)/2     ...         % entropy q(p)
            + spm_logdet(qh.c)/2     ...         % entropy q(h)
            - spm_logdet(pp.c)/2     ...         % entropy - prior p
            - spm_logdet(ph.c)/2     ...         % entropy - prior h     
            + spm_logdet(iS)*nt/(2*r) ...        % entropy - error
            - ne*nt*log(2*pi)/(2*r);
 
        % if F is increasing, save expansion point and dervatives
        %------------------------------------------------------------------
        if L > Fe

            Fe  = L;
            te  = te*2;
            save tempE dFdp dFdpp qp mp
            
        else
            
            % otherwise, return to previous expansion point
            %--------------------------------------------------------------
            load tempE
            te  = min(te/2,1/4);
        end
 
        % E-step: update expectation (p)
        %==================================================================

        % update conditional expectation
        %------------------------------------------------------------------
        dp   = spm_dx(dFdpp,dFdp,{te});
        qp.e = qp.e + dp;
        mp   = mp   + dp;
        
        % convergence (E-Step)
        %------------------------------------------------------------------
        if (dFdp'*dp < 1e-2) | (norm(dp,1) < TOL), break, end
        
        % report (E-Steps)
        %------------------------------------------------------------------
        str{1} = 'E-Step: ';
        str{2} = sprintf('%i',iE);
        str{3} = sprintf('%i',iD);
        str{4} = sprintf('F:%.6e',full(Fe));
        str{5} = sprintf('p:%.2e',full(dp'*dp));
        fprintf('%-8s%-4s%-4s%-20s%-16s\n',str{1:5})
        
    end % E-Step
    
  
    % M-step - hyperparameters (h = exp(l))
    %======================================================================
    mh     = zeros(nh,1);
    dFdh   = zeros(nh,1);
    dFdhh  = zeros(nh,nh);
    
    for iM = 1:nM
 
        % [re-]set precisions using ReML hyperparameter estimates
        %------------------------------------------------------------------
        iC    = P;
        for i = 1:nh
           iC = iC + Q{i}*exp(qh.h(i));
        end
        S     = kron(pinv(R),inv(iC));
 
        % 1st-order derivatives: dFdh = dF/dh
        %------------------------------------------------------------------
        for i = 1:nh
            dPdh{i}        =  kron(R,Q{i}*exp(qh.h(i)));
            dFdh(i,1)      = -trace(dPdh{i}*(EJCJ + EE - S*nt))/(2*r);
        end
 
        % 2nd-order derivatives: dFdhh
        %------------------------------------------------------------------
        for i = 1:nh
            for j = 1:nh
                dFdhh(i,j) =  -trace(dPdh{i}*S*dPdh{j}*S*nt)/(2*r);
            end
        end
        
        % hyperpriors
        %------------------------------------------------------------------
        qh.e  = qh.h  - ph.h;
        dFdh  = dFdh  - ph.ic*qh.e;
        dFdhh = dFdhh - ph.ic;
        
        % update ReML estimate of parameters
        %------------------------------------------------------------------
        dh    = spm_dx(dFdhh,dFdh);
        qh.h  = qh.h + dh;
        mh    = mh   + dh;
 
        % convergence (M-Step)
        %------------------------------------------------------------------
        if (dFdh'*dh < 1e-2) | (norm(dh,1) < TOL), break, end
        
    end % M-Step
 
    % evaluate objective function (F)
    %======================================================================
    qh.ic = -dFdhh;
    qh.c  = inv(qh.ic);
    iS    = kron(R,iC);

    L   = - trace(iS*EE)/(2*r)  ...            % states (u)
          - trace(qp.e'*pp.ic*qp.e)/2  ...     % parameters (p)
          - trace(qh.e'*ph.ic*qh.e)/2  ...     % hyperparameters (h)
          - trace(iS*EJCJ)/(2*r)   ...         % expectation under q
          + spm_logdet(qu.c)/(2*r) ...         % entropy q(u)
          + spm_logdet(qp.c)/2     ...         % entropy q(p)
          + spm_logdet(qh.c)/2     ...         % entropy q(h)
          - spm_logdet(pp.c)/2     ...         % entropy - prior p
          - spm_logdet(ph.c)/2     ...         % entropy - prior h
          + spm_logdet(iS)*nt/(2*r) ...        % entropy - error
          - ne*nt*log(2*pi)/(2*r);

    % if F is increasing, save expansion point and dervatives
    %----------------------------------------------------------------------
    if L > (Fm + 1e-2)

        Fm    = L;
        F(iI) = Fm;

        % save model-states (for each time point)
        %==================================================================
        for t = 1:length(qU)
            v     = spm_unvec(qU(t).v{1},v);
            x     = spm_unvec(qU(t).x{1},x);
            e     = spm_unvec(qU(t).e{1},e);
            for i = 1:(nl - 1)
                QU.v{i + 1}(:,t) = spm_vec(v{i});
                try
                    QU.x{i}(:,t) = spm_vec(x{i});
                end
                QU.e{i}(:,t)     = spm_vec(e{i});
            end
            QU.v{1}(:,t)         = spm_vec(qU(t).y{1} - e{1});
            QU.e{nl}(:,t)        = spm_vec(e{nl});

            % and conditional covariances
            %--------------------------------------------------------------
            QU.C{t} = Le*qU(t).C*Le';
        end

        save tempM

        % report and break if convergence
        %------------------------------------------------------------------
        figure(Fdem)
        spm_DEM_qU(QU)
        if np
            subplot(nl,4,4*nl)
            bar(full(Up*qp.e(1:np,1)))
            xlabel({'parameters';'{minus prior}'})
            axis square, grid on
        end
        if length(F) > 2
            subplot(nl,4,4*nl - 1)
            plot(F(2:end))
            xlabel('iteractions')
            title('Log-evidence')
            axis square, grid on
        end
        drawnow

        % report (EM-Steps)
        %------------------------------------------------------------------   
        str{1} = 'M-Step: ';
        str{2} = sprintf('%i',iI);
        str{3} = sprintf('%i',iE);
        str{4} = sprintf('F:%.6e',full(Fm));
        str{5} = sprintf('p:%.2e',full(mp'*mp));
        str{6} = sprintf('h:%.2e',full(mh'*mh));
        fprintf('%-8s%-4s%-4s%-20s%-16s%-16s\n',str{1:6})

    else

        % otherwise, return to previous expansion point and break
        %------------------------------------------------------------------
        load tempM
        break
    end
end
 
% Assemble output arguments
%==========================================================================

% conditional moments of model-parameters (rotated into original space)
%--------------------------------------------------------------------------
qP.P  = Up*qp.e(1:np) + spm_vec(M.pE);
qP.Pi = spm_unvec(qP.P,M.pE);
qP.C  = Up*qp.c(1:np,1:np)*Up';
 
% conditional moments of hyper-parameters (log-transformed)
%--------------------------------------------------------------------------
qH.h  = qh.h;
qH.hi = spm_unvec(qH.h,M.h);
qH.C  = qh.c;
 
qH.iC = iC;
qH.iV = R;
 
% assign output variables
%--------------------------------------------------------------------------
if nargout == 1
 
    DEM.U  = U;                   % causes
    DEM.X  = X;                   % confounds
 
    DEM.qU = QU;                  % conditional moments of states
    DEM.qP = qP;                  % conditional moments of model-parameters
    DEM.qH = qH;                  % conditional moments of hyper-parameters
 
    DEM.F  = F;                   % [-ve] Free energy
    QU     = DEM;
end
 
warning on
