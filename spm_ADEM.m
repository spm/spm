function [DEM] = spm_ADEM(DEM)
% Dynamic expectation maximisation:  Active inversion
% FORMAT DEM   = spm_ADEM(DEM)
%
% DEM.G  - generative process
% DEM.M  - recognition  model
% DEM.C  - causes
% DEM.U  - prior expectation of causes
%__________________________________________________________________________
%
% This implementation of DEM is the same as spm_DEM but integrates both the
% generative and inversion models in parallel. Its functionality is exactly
% the same apart from the fact that confounds are not accommodated
% explicitly.  The generative model is specified by DEM.G and the veridical
% causes by DEM.C; these may or may not be used as priors on the causes for
% the inversion model DEM.M (i.e., DEM.U = DEM.C).  Clearly, DEM.G does not
% require any priors or precision components; it will use the values of the
% parameters specified in the prior expectation fields.
%
% This routine is not used for model inversion per se but to simulate the
% dynamical inversion of models.  Critically, it includes action 
% variables a that couple the model back to the generative process 
% This enables active inverse (c.f., action-perception; see also spm_GDEM)
%
% hierarchical models G(i) and M(i)
%--------------------------------------------------------------------------
%   M(i).g  = y(t)  = g(x,v,[a],P)    {inline function, string or m-file}
%   M(i).f  = dx/dt = f(x,v,[a],P)    {inline function, string or m-file}
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h hyper-parameters (cause noise)
%   M(i).hC = prior covariances of h hyper-parameters (cause noise)
%   M(i).gE = prior expectation of g hyper-parameters (state noise)
%   M(i).gC = prior covariances of g hyper-parameters (state noise)
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%   M(i).xP = precision (states)
%   
%   M(1).Ra = indices of prediction errors driving action
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i)
%   M(i).l  = number of output v(i)
%  [M(i).k  = number of action a(i)]
%
% Returns the following fields of DEM
%--------------------------------------------------------------------------
%
% true model-states - u
%--------------------------------------------------------------------------
%   pU.x    = hidden states
%   pU.v    = causal states v{1} = response (Y)
%
% model-parameters - p
%--------------------------------------------------------------------------
%   pP.P    = parameters for each level
%
% hyper-parameters (log-transformed) - h,g
%--------------------------------------------------------------------------
%   pH.h    = cause noise
%   pH.g    = state noise
%
% conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.a    = Action
%   qU.x    = Conditional expectation of hidden states
%   qU.v    = Conditional expectation of causal states
%   qU.z    = Conditional prediction errors (v)
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
% represents an alternative to Kalman filtering or alternative Bayesian
% update procedures.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ADEM.m 4322 2011-05-04 15:28:08Z karl $
 
% check model, data, priors and unpack
%--------------------------------------------------------------------------
DEM   = spm_ADEM_set(DEM);
M     = DEM.M;
G     = DEM.G;
C     = DEM.C;
U     = DEM.U;
 
% ensure embedding dimensions are compatible
%--------------------------------------------------------------------------
G(1).E.n = M(1).E.n;
G(1).E.d = M(1).E.n;
 
% find or create a DEM figure
%--------------------------------------------------------------------------
Fdem = spm_figure('GetWin','DEM');
 
% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d + 1;                   % embedding order of q(v)
n    = M(1).E.n + 1;                   % embedding order of q(x) (n >= d)
s    = M(1).E.s;                       % smoothness - s.d. of kernel (bins)
 
% number of states and parameters - recognition model
%--------------------------------------------------------------------------
nY   = size(C,2);                      % number of samples
nl   = size(M,2);                      % number of levels
nv   = sum(spm_vec(M.m));              % number of v (causal states)
nx   = sum(spm_vec(M.n));              % number of x (hidden states)
ny   = M(1).l;                         % number of y (inputs)
nc   = M(end).l;                       % number of c (prior causes)
nu   = nv*d + nx*n;                    % number of generalised states q(u)
 
% number of states and parameters - generative model
%--------------------------------------------------------------------------
gr   = sum(spm_vec(G.l));              % number of v (outputs)
gv   = sum(spm_vec(G.m));              % number of v (causal states)
ga   = sum(spm_vec(G.k));              % number of a (active states)
gx   = sum(spm_vec(G.n));              % number of x (hidden states)
gy   = G(1).l;                         % number of y (inputs)
na   = ga;
 
% number of iterations
%--------------------------------------------------------------------------
try, nE = M(1).E.nE; catch, nE = 16; end
try, nM = M(1).E.nM; catch, nM = 8;  end
 
 
% initialise regularisation parameters
%--------------------------------------------------------------------------
td = 1;                               % integration time for D-Step
te = 2;                               % integration time for E-Step (log)
 
 
% precision (R) and covariance of generalised errors
%--------------------------------------------------------------------------
iV    = spm_DEM_R(n,s);
 
% precision components Q{} requiring [Re]ML estimators (M-Step)
%==========================================================================
Q     = {};
for i = 1:nl
    q0{i,i} = sparse(M(i).l,M(i).l);
    r0{i,i} = sparse(M(i).n,M(i).n);
end
Q0    = kron(iV,spm_cat(q0));
R0    = kron(iV,spm_cat(r0));
for i = 1:nl
    for j = 1:length(M(i).Q)
        q          = q0;
        q{i,i}     = M(i).Q{j};
        Q{end + 1} = blkdiag(kron(iV,spm_cat(q)),R0);
    end
    for j = 1:length(M(i).R)
        q          = r0;
        q{i,i}     = M(i).R{j};
        Q{end + 1} = blkdiag(Q0,kron(iV,spm_cat(q)));
    end
end
 
 
% and fixed components P
%--------------------------------------------------------------------------
Q0    = kron(iV,spm_cat(spm_diag({M.V})));
R0    = kron(iV,spm_cat(spm_diag({M.W})));
Qp    = blkdiag(Q0,R0);
nh    = length(Q);                         % number of hyperparameters
 
% fixed priors on states (u)
%--------------------------------------------------------------------------
xP    = spm_cat(spm_diag({M.xP}));
Px    = kron(iV(1:n,1:n),speye(nx,nx)*exp(-8) + xP);
Pv    = kron(iV(1:d,1:d),speye(nv,nv)*exp(-8));
Pa    = kron(iV(1:1,1:1),speye(na,na)*exp(-8));
Pu    = spm_cat(spm_diag({Px Pv}));
 
% hyperpriors
%--------------------------------------------------------------------------
ph.h  = spm_vec({M.hE M.gE});              % prior expectation of h
ph.c  = spm_cat(spm_diag({M.hC M.gC}));    % prior covariances of h
qh.h  = ph.h;                              % conditional expectation
qh.c  = ph.c;                              % conditional covariance
ph.ic = spm_inv(ph.c);                     % prior precision
 
% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
 
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC);                    % basis for parameters
    M(i).p    = size(qp.u{i},2);                     % number of qp.p
    qp.p{i}   = sparse(M(i).p,1);                    % initial qp.p
    pp.c{i,i} = qp.u{i}'*M(i).pC*qp.u{i};            % prior covariance
    
    try
        qp.e{i} = qp.p{i} + qp.u{i}'*(spm_vec(M(i).P) - spm_vec(M(i).pE));
    catch
        qp.e{i} = qp.p{i};                           % initial qp.e
    end
 
end
Up    = spm_cat(spm_diag(qp.u));
 
% initialise and augment with confound parameters B; with flat priors
%--------------------------------------------------------------------------
np    = sum(spm_vec(M.p));                  % number of model parameters
pp.c  = spm_cat(pp.c);
pp.ic = spm_inv(pp.c);
 
% initialise conditional density q(p) (for D-Step)
%--------------------------------------------------------------------------
qp.e  = spm_vec(qp.e);
qp.c  = sparse(np,np);
 
% initialise cell arrays for D-Step; e{i + 1} = (d/dt)^i[e] = e[i]
%==========================================================================
qu.x      = cell(n,1);
qu.v      = cell(n,1);
qu.a      = cell(n,1);
qu.y      = cell(n,1);
qu.u      = cell(n,1);
pu.v      = cell(n,1);
pu.x      = cell(n,1);
pu.z      = cell(n,1);
pu.w      = cell(n,1);
 
[qu.x{:}] = deal(sparse(nx,1));
[qu.v{:}] = deal(sparse(nv,1));
[qu.a{:}] = deal(sparse(na,1));
[qu.y{:}] = deal(sparse(ny,1));
[qu.u{:}] = deal(sparse(nc,1));
[pu.v{:}] = deal(sparse(gr,1));
[pu.x{:}] = deal(sparse(gx,1));
[pu.z{:}] = deal(sparse(gr,1));
[pu.w{:}] = deal(sparse(gx,1));
 
% initialise cell arrays for hierarchical structure of x[0] and v[0]
%--------------------------------------------------------------------------
qu.x{1}   = spm_vec({M(1:end - 1).x});
qu.v{1}   = spm_vec({M(1 + 1:end).v});
qu.a{1}   = spm_vec({G.a});
pu.x{1}   = spm_vec({G.x});
pu.v{1}   = spm_vec({G.v});
 
 
% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dx    = kron(spm_speye(n,n,1),spm_speye(nx,nx,0));
Dv    = kron(spm_speye(d,d,1),spm_speye(nv,nv,0));
Dc    = kron(spm_speye(d,d,1),spm_speye(nc,nc,0));
Da    = kron(spm_speye(1,1,1),spm_speye(na,na,0));
Du    = spm_cat(spm_diag({Dx,Dv}));
Dq    = spm_cat(spm_diag({Dx,Dv,Dc,Da}));
 
Dx    = kron(spm_speye(n,n,1),spm_speye(gx,gx,0));
Dv    = kron(spm_speye(n,n,1),spm_speye(gr,gr,0));
Dp    = spm_cat(spm_diag({Dv,Dx,Dv,Dx}));
dfdw  = kron(speye(n,n),speye(gx,gx));
dydv  = kron(speye(n,n),speye(gy,gr));

% reistiction matrix, mapping prediction errors to action
%--------------------------------------------------------------------------
try
    Ra = M(1).Ra;
catch
    Ra = 1:ny;
end
Ra    = kron(speye(n,n),sparse(Ra,Ra,1,ny,ny));
 
% and null blocks
%--------------------------------------------------------------------------
dVdc  = sparse(d*nc,1);
 
% gradients and curvatures for conditional uncertainty
%--------------------------------------------------------------------------
dWdu  = sparse(nu,1);
dWduu = sparse(nu,nu);

% preclude unnecessary iterations
%--------------------------------------------------------------------------
if ~np && ~nh, nE = 1; end
 
 
% create innovations (and add causes)
%--------------------------------------------------------------------------
[z w]  = spm_DEM_z(G,nY);
z{end} = C + z{end};
Z      = spm_cat(z(:));
W      = spm_cat(w(:));
 
% Iterate DEM
%==========================================================================
F      = -Inf;
for iE = 1:nE
    
    % get time and celar persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic; clear spm_DEM_eval
 
    % E-Step: (with embedded D-Step)
    %======================================================================
 
    % [re-]set accumulators for E-Step
    %----------------------------------------------------------------------
    dFdp  = zeros(np,1);
    dFdpp = zeros(np,np);
    EE    = sparse(0);
    ECE   = sparse(0);
    qp.ic = sparse(0);
    Hqu.c = sparse(0);

 
    % [re-]set precisions using ReML hyperparameter estimates
    %----------------------------------------------------------------------
    iS    = Qp;
    for i = 1:nh
       iS = iS + Q{i}*exp(qh.h(i));
    end
 
    % [re-]set states & their derivatives
    %----------------------------------------------------------------------
    try
        qu = qU(1);
        pu = pU(1);
    end
 
    % D-Step: (nD D-Steps for each sample)
    %======================================================================
    for iY = 1:nY
 
 
        % D-Step: until convergence for static systems
        %==================================================================
 
        % derivatives of responses and inputs
        %------------------------------------------------------------------
        pu.z = spm_DEM_embed(Z,n,iY);
        pu.w = spm_DEM_embed(W,n,iY);
        qu.u = spm_DEM_embed(U,n,iY);
        
        % pass action to pu.a
        %------------------------------------------------------------------
        pu.a = qu.a;
 
        % evaluate generative model
        %------------------------------------------------------------------
        [pu dg df] = spm_ADEM_diff(G,pu);
 
        % tensor products for Jacobian
        %------------------------------------------------------------------
        Dgda = kron(spm_speye(n,1,1),dg.da);
        Dgdv = kron(spm_speye(n,n,1),dg.dv);
        Dgdx = kron(spm_speye(n,n,1),dg.dx);
        dfda = kron(spm_speye(n,1,0),df.da);
        dfdv = kron(spm_speye(n,n,0),df.dv);
        dfdx = kron(spm_speye(n,n,0),df.dx);
        
        dgda = kron(spm_speye(n,1,0),dg.da);
        dgdx = kron(spm_speye(n,n,0),dg.dx);
        
        % and pass response to qu.y
        %------------------------------------------------------------------
        for i = 1:n
            y       = spm_unvec(pu.v{i},{G.v});
            qu.y{i} = y{1};
        end
 
        % evaluate recognition model
        %------------------------------------------------------------------       
        [E dE] = spm_DEM_eval(M,qu,qp);
 
        
        % conditional covariance [of states {u}]
        %------------------------------------------------------------------
        qu.c   = spm_inv(dE.du'*iS*dE.du + Pu);
        Hqu.c  = Hqu.c + spm_logdet(qu.c);
        
        % save at qu(t)
        %------------------------------------------------------------------
        qE{iY} = E;
        qC{iY} = qu.c;
        qU(iY) = qu;
        pU(iY) = pu;
 
        % and conditional covariance [of parameters {P}]
        %------------------------------------------------------------------
        ECEu   = dE.du*qu.c*dE.du';
        ECEp   = dE.dp*qp.c*dE.dp';
 
        
        % uncertainty about parameters dWdv, ... ; W = ln(|qp.c|)
        %==================================================================
        if np
            for i = 1:nu
                CJp(:,i)   = spm_vec(qp.c*dE.dpu{i}'*iS);
                dEdpu(:,i) = spm_vec(dE.dpu{i}');
            end
            dWdu  = CJp'*spm_vec(dE.dp');
            dWduu = CJp'*dEdpu;
        end
        
        % change in error w.r.t. action
        %------------------------------------------------------------------
        dyda  = dydv*dgda;
        dydx  = dydv*dgdx;
        for i = 1:n
            Dfdx = kron(spm_speye(n,n,-i),df.dx^(i - 1));
            dyda = dyda + dydx*Dfdx*dfda;
        end
        
        % dE/da with restriction
        %------------------------------------------------------------------
        dE.da = dE.dy*Ra*dyda;
        dE.dv = dE.dy*dydv;
        
        % first-order derivatives
        %------------------------------------------------------------------
        dVdu  = -dE.du'*iS*E - Pu*spm_vec({qu.x{1:n} qu.v{1:d}}) - dWdu/2;
        dVda  = -dE.da'*iS*E - Pa*spm_vec( qu.a{1:1});
        
        % and second-order derivatives
        %------------------------------------------------------------------
        dVduu = -dE.du'*iS*dE.du - Pu - dWduu/2 ;
        dVdaa = -dE.da'*iS*dE.da - Pa;
        dVduv = -dE.du'*iS*dE.dv;
        dVduc = -dE.du'*iS*dE.dc;
        dVdua = -dE.du'*iS*dE.da;
        dVdav = -dE.da'*iS*dE.dv;
        dVdau = -dE.da'*iS*dE.du;
        dVdac = -dE.da'*iS*dE.dc;

        
 
         
        % D-step update: of causes v{i}, and hidden states x(i)
        %==================================================================
 
        % states and conditional modes
        %------------------------------------------------------------------
        p     = {pu.v{1:n} pu.x{1:n} pu.z{1:n} pu.w{1:n}};
        q     = {qu.x{1:n} qu.v{1:d} qu.u{1:d} qu.a{1:1}};
        u     = [p q];  
        
        % gradient
        %------------------------------------------------------------------
        dFdu  = [                              Dp*spm_vec(p); 
                 spm_vec({dVdu; dVdc; dVda}) + Dq*spm_vec(q)];
 
 
        % Jacobian (variational flow)
        %------------------------------------------------------------------
        dFduu = spm_cat(...
                {Dgdv  Dgdx Dv   []   []       []    Dgda;
                 dfdv  dfdx []   dfdw []       []    dfda;
                 []    []   Dv   []   []       []    [];
                 []    []   []   Dx   []       []    [];
                 dVduv []   []   []   Du+dVduu dVduc dVdua;
                 []    []   []   []   []       Dc    []
                 dVdav []   []   []   dVdau    dVdac dVdaa});
 
 
        % update states q = {x,v,z,w} and conditional modes
        %==================================================================
        du    = spm_dx(dFduu,dFdu,td);
        u     = spm_unvec(spm_vec(u) + du,u);
 
        % and save them
        %------------------------------------------------------------------
        pu.v(1:n) = u([1:n]);
        pu.x(1:n) = u([1:n] + n);
        qu.x(1:n) = u([1:n] + n + n + n + n);
        qu.v(1:d) = u([1:d] + n + n + n + n + n);
        qu.a(1:1) = u([1:1] + n + n + n + n + n + d + d);
 
        % Gradients and curvatures for E-Step: W = tr(C*J'*iS*J)
        %==================================================================
        if np
            for i = 1:np
                CJu(:,i)   = spm_vec(qu.c*dE.dup{i}'*iS);
                dEdup(:,i) = spm_vec(dE.dup{i}');
            end
            dWdp  = CJu'*spm_vec(dE.du');
            dWdpp = CJu'*dEdup;

            % Accumulate; dF/dP = <dL/dp>, dF/dpp = ...
            %--------------------------------------------------------------
            dFdp  = dFdp  - dWdp/2  - dE.dp'*iS*E;
            dFdpp = dFdpp - dWdpp/2 - dE.dp'*iS*dE.dp;
            qp.ic = qp.ic           + dE.dp'*iS*dE.dp;
            
        end
 
        % and quantities for M-Step
        %------------------------------------------------------------------
        EE  = E*E'+ EE;
        ECE = ECE + ECEu + ECEp;
 
    end % sequence (nY)
 
    % augment with priors
    %----------------------------------------------------------------------
    dFdp   = dFdp  - pp.ic*qp.e;
    dFdpp  = dFdpp - pp.ic;
    qp.ic  = qp.ic + pp.ic;
    qp.c   = spm_inv(qp.ic);
 
 
    % E-step: update expectation (p)
    %======================================================================
 
    % update conditional expectation
    %----------------------------------------------------------------------
    dp     = spm_dx(dFdpp,dFdp,{te});
    qp.e   = qp.e + dp;
    qp.p   = spm_unvec(qp.e,qp.p);
 
 
    % M-step - hyperparameters (h = exp(l))
    %======================================================================
    mh     = zeros(nh,1);
    dFdh   = zeros(nh,1);
    dFdhh  = zeros(nh,nh);
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
            dPdh{i}        =  Q{i}*exp(qh.h(i));
            dFdh(i,1)      = -trace(dPdh{i}*dS)/2;
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
        dh    = spm_dx(dFdhh,dFdh);
        qh.h  = qh.h + dh;
        mh    = mh   + dh;
 
        % conditional covariance of hyperparameters
        %------------------------------------------------------------------
        qh.c  = -spm_inv(dFdhh);
 
        % convergence (M-Step)
        %------------------------------------------------------------------
        if (dFdh'*dh < 1e-2) || (norm(dh,1) < exp(-8)), break, end
 
    end % M-Step
 
    % evaluate objective function (F)
    %======================================================================
    L   = - trace(iS*EE)/2  ...              % states (u)
        - trace(qp.e'*pp.ic*qp.e)/2  ...     % parameters (p)
        - trace(qh.e'*ph.ic*qh.e)/2  ...     % hyperparameters (h)
        + Hqu.c/2             ...            % entropy q(u)
        + spm_logdet(qp.c)/2  ...            % entropy q(p)
        + spm_logdet(qh.c)/2  ...            % entropy q(h)
        - spm_logdet(pp.c)/2  ...            % entropy - prior p
        - spm_logdet(ph.c)/2  ...            % entropy - prior h
        + spm_logdet(iS)*nY/2 ...            % entropy - error
        - n*ny*nY*log(2*pi)/2;
    
    
    % if F is increasing, save expansion point and derivatives
    %----------------------------------------------------------------------
    if L > F(end) || iE < 3
   
        % save model-states (for each time point)
        %==================================================================
        for t = 1:length(qU)
 
            % states
            %--------------------------------------------------------------
            a     = spm_unvec(qU(t).a{1},{G.a});
            v     = spm_unvec(pU(t).v{1},{G.v});
            x     = spm_unvec(pU(t).x{1},{G.x});
            z     = spm_unvec(pU(t).z{1},{G.v});
            w     = spm_unvec(pU(t).w{1},{G.x});
            for i = 1:nl
                try
                    PU.v{i}(:,t) = spm_vec(v{i});
                    PU.z{i}(:,t) = spm_vec(z{i});
                end
                try
                    PU.x{i}(:,t) = spm_vec(x{i});
                    PU.w{i}(:,t) = spm_vec(w{i});
                end
                try
                    QU.a{i}(:,t) = spm_vec(a{i});
                end
            end
 
            % conditional modes
            %--------------------------------------------------------------
            v     = spm_unvec(qU(t).v{1},{M(1 + 1:end).v});
            x     = spm_unvec(qU(t).x{1},{M(1:end - 1).x});
            z     = spm_unvec(qE{t}(1:(ny + nv)),{M.v});
            w     = spm_unvec(qE{t}((1:nx) + (ny + nv)*n),{M.x});
            for i = 1:(nl - 1)
                if M(i).m, QU.v{i + 1}(:,t) = spm_vec(v{i}); end
                if M(i).l, QU.z{i}(:,t)     = spm_vec(z{i}); end
                if M(i).n, QU.x{i}(:,t)     = spm_vec(x{i}); end
                if M(i).n, QU.w{i}(:,t)     = spm_vec(w{i}); end
            end
            QU.v{1}(:,t)  = spm_vec(qU(t).y{1}) - spm_vec(z{1});
            QU.z{nl}(:,t) = spm_vec(z{nl});
 
            % and conditional covariances
            %--------------------------------------------------------------
            i       = (1:nx);
            QU.S{t} = qC{t}(i,i);
            i       = (1:nv) + nx*n;
            QU.C{t} = qC{t}(i,i);
        end
 
        % save conditional densities
        %------------------------------------------------------------------
        B.QU  = QU;
        B.PU  = PU;
        B.qp  = qp;
        B.qh  = qh;
 
        % decrease regularisation
        %------------------------------------------------------------------
        F(iE) = L;
        te    = min(te + 1,8);
 
    else
 
        % otherwise, return to previous expansion point and break
        %------------------------------------------------------------------
        QU    = B.QU;
        PU    = B.PU;
        qp    = B.qp;
        qh    = B.qh;
 
        % increase regularisation
        %------------------------------------------------------------------
        F(iE) = F(end);
        te    = min(te - 1,0);
        
    end
 
    % report and break if convergence
    %======================================================================
    figure(Fdem)
    spm_DEM_qU(QU)
    if np
        subplot(nl,4,4*nl)
        bar(full(Up*qp.e))
        xlabel({'parameters';'{minus prior}'})
        axis square, grid on
    end
    if length(F) > 2
        subplot(nl,4,4*nl - 1)
        plot(F - F(1))
        xlabel('updates')
        title('log-evidence')
        axis square, grid on
    end
    drawnow
 
    % report (EM-Steps)
    %----------------------------------------------------------------------
    str{1} = sprintf('ADEM: %i (%i)',iE,iM);
    str{2} = sprintf('F:%.4e',full(L - F(1)));
    str{3} = sprintf('p:%.2e',full(dp'*dp));
    str{4} = sprintf('h:%.2e',full(mh'*mh));
    str{5} = sprintf('(%.2e sec)',full(toc));
    fprintf('%-16s%-16s%-14s%-14s%-16s\n',str{:})
    
    if (norm(dp,1) < exp(-8)) && (norm(mh,1) < exp(-8)), break, end
 
end
 
% assemble output arguments
%==========================================================================
 
% Fill in DEM with response and its causes
%--------------------------------------------------------------------------
DEM.Y    = PU.v{1};
DEM.pU   = PU;
DEM.pP.P = {G.pE};
DEM.pH.h = {G.hE};
DEM.pH.g = {G.gE};
 
% conditional moments of model-parameters (rotated into original space)
%--------------------------------------------------------------------------
qP.P   = spm_unvec(Up*qp.e + spm_vec(M.pE),M.pE);
qP.C   = Up*qp.c*Up';
qP.V   = spm_unvec(diag(qP.C),M.pE);
 
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
 
DEM.qU = QU;                  % conditional moments of model-states
DEM.qP = qP;                  % conditional moments of model-parameters
DEM.qH = qH;                  % conditional moments of hyper-parameters
 
DEM.F  = F;                   % [-ve] Free energy

