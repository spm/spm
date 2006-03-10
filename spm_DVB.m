function [DEM] = spm_DVB(DEM)
% Dynamic varaitional Bayes (c.f., DEM but with no system noise)
% FORMAT DEM   = spm_DVB(DEM)
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
%   M(i).V  = fixed precision
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

% Check model, data, priros and confounds and unpack
%--------------------------------------------------------------------------
[M Y U X] = spm_DEM_set(DEM);

% find or create a DEM figure
%--------------------------------------------------------------------------
warning off
Fdem     = spm_figure('GetWin','DEM');
if isempty(Fdem)
    name = 'Dynamic Expectation Maximisation';
    Fdem = spm_figure('CreateWin','DEM',name,'on');
end

% tolerance for changes in norm
%--------------------------------------------------------------------------
TOL  = 1e-3;

% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d;                       % truncation order of q(x,v)
n    = M(1).E.n;                       % truncation order of e (n >= d)
s    = M(1).E.s;                       % smoothness - s.d. of kernel (secs)
dt   = M(1).E.dt;                      % time bins {secs}


% number of states and parameters
%--------------------------------------------------------------------------
ns   = size(Y,2);                       % number of samples
nl   = size(M,2);                       % number of levels
ne   = sum(cat(1,M.l));                 % number of e (errors)
nv   = sum(cat(1,M.m));                 % number of v (casual states)
nx   = sum(cat(1,M.n));                 % number of x (hidden states)
ny   = M(1).l;                          % number of y (inputs)
nc   = M(end).l;                        % number of c (prior causes)
nu   = nv*d + nx;                       % number of generalised states
kt   = 1;                               % rate constant for D-Step

% number of iterations
%--------------------------------------------------------------------------
try nD = M(1).E.nD; catch nD = 8; end
try nE = M(1).E.nE; catch nE = 4; end
try nM = M(1).E.nM; catch nM = 8; end
try nI = M(1).E.nI; catch nI = 16;end

% initialise regularisation parameters
%--------------------------------------------------------------------------
if nx
    td = dt/nD;                           % integration time for D-Step
    te = 64;                              % integration time for E-Step
else
    td = {64};
    te = 64;
end
 
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
nh    = length(Q);                         % number of hyperparameters

% hyperpriors
%--------------------------------------------------------------------------
try
    ph.h = spm_vec(M.hE);                  % prior expectation of h
    ph.c = spm_cat(diag({M.hC}));          % prior covariances of h
catch
    ph.h = sparse(nh,1);                   % prior expectation of h
    ph.c = speye(nh,nh)*16;                % prior covariances of h
end
ph.ic    = inv(ph.c);                      % prior precision
qh.h     = spm_vec(M.h);                   % conditional expectation
qh.c     = ph.c;                           % conditional covariance
qh.e     = qh.h - ph.h;
 

% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
 
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC);                    % basis for parameters
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
ip    = [1:np];
ib    = [1:nn] + np;
pp.c  = spm_cat(pp.c);
pp.ic = inv(pp.c);
 
% initialise conditional density q(p) (for D-Step)
%--------------------------------------------------------------------------
qp.e  = spm_vec(qp.p);
qp.c  = sparse(nf,nf);
qb    = sparse(nn,1 );

% initialise dedb
%--------------------------------------------------------------------------
for i = 1:nl
    dedbi{i,1} = sparse(M(i).l,nn);
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

dq        = {qu.v{1:d} qu.x{1} qu.y{1:n} qu.u{1:d}};
 
% initialise cell arrays for hierarchical structure of x[0] and v[0]
%--------------------------------------------------------------------------
v         = {M(2:end).v};
x         = {M.x};
e         = {M.v};
qu.x{1}   = spm_vec(x);
qu.v{1}   = spm_vec(v);
qU        = qu;
qU.c      = [];

% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dv        = cell(d,d);
Dx        = cell(1,1);
Dy        = cell(n,n);
Dc        = cell(d,d);
[Dv{:}]   = deal(sparse(nv,nv));
[Dx{:}]   = deal(sparse(nv,nx));
[Dy{:}]   = deal(sparse(ny,ny));
[Dc{:}]   = deal(sparse(nc,nc));
 
% add constant terms
%--------------------------------------------------------------------------
for i = 2:d
    Dv{i - 1,i} = speye(nv,nv);
    Dc{i - 1,i} = speye(nc,nc);
end
for i = 2:n
    Dy{i - 1,i} = speye(ny,ny);
end
Du        = spm_cat({Dv Dx});
Dy        = spm_cat( Dy);
Dc        = spm_cat( Dc);


% gradients and curvatures for conditional uncertainty
%--------------------------------------------------------------------------
dUdu      = sparse(nu,1);
dUdp      = sparse(nf,1);
dUduu     = sparse(nu,nu);
dUdpp     = sparse(nf,nf);

% preclude unneceassry iterations
%--------------------------------------------------------------------------
if ~nh,        nM = 1; end
if ~nf,        nE = 1; end
if ~nf && ~nh, nI = 1; end

%  Precision (R) and covariance of generalised errors
%--------------------------------------------------------------------------
[R V]  = spm_DEM_R(n,s);
if s < 1/2
    Rm = R;
else
    Rm = sparse(1,1,R(1),n,n);
end


% Iterate DEM
%==========================================================================
Fm     = -exp(64);
for iI = 1:nI
 
    % E-Step: (with embedded D-Step)
    %======================================================================
    mp     = zeros(nf,1);
    Fe     = -exp(64);
    for iE = 1:nE
 

        % [re-]set accumulators for E-Step
        %------------------------------------------------------------------
        dFdp  = zeros(nf,1);
        dFdpp = zeros(nf,nf);
        EE    = sparse(0);
        JCJ   = sparse(0);
        qp.ic = sparse(0);
        qu_c  = speye(1);
 
        
        % [re-]set hierarchical parameters
        %------------------------------------------------------------------
        qp.p  = spm_unvec(qp.e,qp.p);
        
        
        % [re-]set precisions using ReML hyperparameter estimates
        %------------------------------------------------------------------
        iC    = P;
        for i = 1:nh
           iC = iC + Q{i}*exp(qh.h(i));
        end
        iS    = kron(R,iC);
 
        % [re-]adjust for confounds
        %------------------------------------------------------------------
        Y     = Y - reshape(qb,ny,nb)*X;
        
        % [re-]set states & their derivatives
        %------------------------------------------------------------------
        qu    = qU(1);
        
        % D-Step: (nD D-Steps for each sample)
        %==================================================================
        for t = 1:ns

            % [re-]set states for static systems
            %--------------------------------------------------------------
            if ~nx
                try, qu = qU(t); end
            end

            % D-Step: until convergence for static systems
            %==============================================================
            Fd     = -exp(64);
            for iD = 1:nD
                
                % sampling time
                %----------------------------------------------------------
                ts        = (t + (iD - 1)/nD)*dt;
                
                % derivatives of responses and inputs
                %----------------------------------------------------------
                qu.y(1:n) = spm_DEM_embed(Y,n,ts,dt);
                qu.u(1:d) = spm_DEM_embed(U,d,ts,dt);

                % compute dEdb (derivatives of confounds)
                %----------------------------------------------------------
                b     = spm_DEM_embed(X,n,ts,dt);
                for i = 1:n
                    dedbi{1}  = -kron(b{i}',speye(ny,ny));
                    dEdb{i,1} =  spm_cat(dedbi);
                end

                % evaluate functions:
                % e = v - g(x,v), dx/dt = f(x,v) and derivatives dE.dx, ...
                %==========================================================       
                [qu df dE] = spm_DEM_diff(M,qu,qp);
                
                % and conditional covariance [of states {v}]
                %----------------------------------------------------------
                qu.c       = inv(dE.du'*iS*dE.du);
                
                % save states at qu(t)
                %----------------------------------------------------------
                if iD == 1, qU(t)  = qu; end
                
                % vectorise error
                %----------------------------------------------------------
                dE.dP  = [dE.dp spm_cat(dEdb)];
                E      = spm_vec(qu.e);
                JCJu   = dE.du*qu.c*dE.du';
                JCJp   = dE.dP*qp.c*dE.dP';
          
                % and evaluate objective function L(t) (for static models)
                %----------------------------------------------------------
                if ~nx
                    
                    L = - trace(E'*iS*E)/2 ...           % states (u)
                        - trace(iS*JCJp)/2;              % expectation q(p)

                    % if F is increasing, save expansion point
                    %------------------------------------------------------
                    if L > Fd
                        td  = {min(td{1}*2,256)};
                        Fd  = L;
                        save tempD qu df dE E JCJu JCJp
                    else
                        % otherwise, return to previous expansion point
                        %--------------------------------------------------
                        load tempD
                        td  = {min(td{1}/2,16)};
                    end
                end

                % conditional uncertainty about parameters
                %==========================================================
                if np
                    for i = 1:nu
                        
                        % 1st-order derivatives: dUdv, ... ; U = ln(|qp.c|)
                        %--------------------------------------------------
                        CJ             = qp.c(ip,ip)*dE.dpu{i}'*iS;
                        dUdu(i,1)      = trace(CJ*dE.dp);

                        % 2nd-order derivatives
                        %--------------------------------------------------
                        for j = 1:nu
                            dUduu(i,j) = trace(CJ*dE.dpu{j});
                        end
                    end
                end


                % D-step update: of causes v{i}, and other states u(i)
                %==========================================================

                % compute dqdt: q = {u y c}; and dudt: u = {v{1:d} x}
                %----------------------------------------------------------
                dIdu  = -dE.du'*iS*E - dUdu/2; 
                
                % and second-order derivatives
                %----------------------------------------------------------
                dIduu = -dE.du'*iS*dE.du - dUduu/2;
                dIduy = -dE.du'*iS*dE.dy;
                dIduc = -dE.du'*iS*dE.dc;
                
                % mean field effects
                %----------------------------------------------------------
                D     = spm_cat({Du; df.du});
                if nx
                    kt = exp(32)*normest(D)/normest(dIduu);
                end
                
                % gradient
                %----------------------------------------------------------
                dFdu  = spm_vec({qu.v(2:d + 1) qu.x(2)}) + kt*dIdu;
                dFdu  = spm_vec({dFdu           ;
                                 qu.y(2:n + 1)  ;
                                 qu.u(2:d + 1)});
                
                
                % Jacobian
                %----------------------------------------------------------
                dFduu = spm_cat({(D + kt*dIduu) kt*dIduy    kt*dIduc;
                                  []            Dy          []      ;
                                  []            []          Dc   }) ;                         
                
      
                % update conditional modes of states
                %----------------------------------------------------------
                du    = spm_dx(dFduu,dFdu,td);
                dq    = spm_unvec(du,dq);
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
                if ~nx, qU(t)  = qu; end
                if ~nx && ((dFdu'*du < 1e-2) | (norm(du,1) < TOL))
                    break
                end
                if ~nx && ns < 8
                    % report (D-Steps)
                    %------------------------------------------------------
                    str{1} = 'D-Step: ';
                    str{2} = sprintf('%i',t);
                    str{3} = sprintf('%i',iD);
                    str{4} = sprintf('I:%.6e',full(Fd));
                    str{5} = sprintf('u:%.2e',full(du'*du));
                    fprintf('%-8s%-4s%-4s%-20s%-16s\n',str{1:5})
                end

            end % D-Step
            
            % Gradients and curvatures for E-Step:
            %==============================================================
            for i = ip
                
                % 1st-order derivatives: U = tr(C*J'*iS*J)
                %----------------------------------------------------------
                CJ             = qu.c*dE.dup{i}'*iS;
                dUdp(i,1)      = trace(CJ*dE.du);
 
                % 2nd-order derivatives
                %----------------------------------------------------------
                for j = ip
                    dUdpp(i,j) = trace(CJ*dE.dup{j});
                end
            end
 
            % Accumulate; dF/dP = <dL/dp>, dF/dpp = ...
            %--------------------------------------------------------------
            dFdp  = dFdp  - dUdp/2  - dE.dP'*iS*E;
            dFdpp = dFdpp - dUdpp/2 - dE.dP'*iS*dE.dP;
            qp.ic = qp.ic           + dE.dP'*iS*dE.dP;
 
            % and quantities for M-Step
            %--------------------------------------------------------------
            qu_c  = qu_c*qu.c;
            EE    = E*E'+ EE;
            JCJ   = JCJ + JCJu + JCJp;
            
        end % sequence (t)
 
        % augment with priors
        %------------------------------------------------------------------
        dFdp(ip)     = dFdp(ip)     - pp.ic*qp.e;
        dFdpp(ip,ip) = dFdpp(ip,ip) - pp.ic;
        qp.ic(ip,ip) = qp.ic(ip,ip) + pp.ic;
        qp.c         = inv(qp.ic);
             
        % evaluate objective function <L(t)>
        %==================================================================
        L = - trace(iS*EE)/2  ...                    % states (u)
            - trace(qp.e'*pp.ic*qp.e)/2;             % parameters (p)

 
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
        qp.e = qp.e + dp(ip);
        db   = dp(ib);
        mp   = mp + dp;
        
        % convergence (E-Step)
        %------------------------------------------------------------------
        if (dFdp'*dp < 1e-2) | (norm(dp,1) < TOL), break, end
        
        % report (E-Steps)
        %------------------------------------------------------------------
        str{1} = 'E-Step: ';
        str{2} = sprintf('%i',iE);
        str{3} = sprintf('%i',iD);
        str{4} = sprintf('I:%.6e',full(Fe));
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
        S     = kron(V,inv(iC));
        dS    = JCJ + EE - S*ns;
         
        % 1st-order derivatives: dFdh = dF/dh
        %------------------------------------------------------------------
        for i = 1:nh
            dPdh{i}        =  kron(Rm,Q{i}*exp(qh.h(i)));
            dFdh(i,1)      = -trace(dPdh{i}*dS)/2;
        end
 
        % 2nd-order derivatives: dFdhh
        %------------------------------------------------------------------
        for i = 1:nh
            for j = 1:nh
                dFdhh(i,j) = -trace(dPdh{i}*S*dPdh{j}*S*ns)/2;
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
        dh    = min(dh, 8);
        dh    = max(dh,-8);
        qh.h  = qh.h + dh;
        mh    = mh   + dh;
        
        % conditional covariance of hyperparameters
        %------------------------------------------------------------------
        if iM == 1
            qh.c = -inv(dFdhh);
        end
        
        % convergence (M-Step)
        %------------------------------------------------------------------
        if (dFdh'*dh < 1e-2) | (norm(dh,1) < TOL), break, end
        
    end % M-Step
 
    % evaluate objective function (F)
    %======================================================================
    L   = - trace(iS*EE)/2  ...                % states (u)
          - trace(qp.e'*pp.ic*qp.e)/2  ...     % parameters (p)
          - trace(qh.e'*ph.ic*qh.e)/2  ...     % hyperparameters (h)
          + spm_logdet(qu_c)/2  ...            % entropy q(u)
          + spm_logdet(qp.c)/2  ...            % entropy q(p)
          + spm_logdet(qh.c)/2  ...            % entropy q(h)
          - spm_logdet(pp.c)/2  ...            % entropy - prior p
          - spm_logdet(ph.c)/2  ...            % entropy - prior h
          + spm_logdet(iS)*ns/2 ...            % entropy - error
          - ny*ns*log(2*pi)/2;

    
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
            i       = [1:nv];
            QU.C{t} = qU(t).c(i,i);
            i       = [1:nx] + nv*d;
            QU.S{t} = qU(t).c(i,i);
        end

        save tempM

        % report and break if convergence
        %------------------------------------------------------------------
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
        str{3} = sprintf('%i',iM);
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
qP.P  = Up*qp.e + spm_vec(M.pE);
qP.Pi = spm_unvec(qP.P,M.pE);
qP.C  = Up*qp.c(ip,ip)*Up';
 
% conditional moments of hyper-parameters (log-transformed)
%--------------------------------------------------------------------------
qH.h  = qh.h;
qH.hi = spm_unvec(qH.h,M.h);
qH.C  = qh.c;
 
qH.iC = iC;
qH.iV = R;
 
% assign output variables
%--------------------------------------------------------------------------
DEM.U  = U;                   % causes
DEM.X  = X;                   % confounds
 
DEM.qU = QU;                  % conditional moments of states
DEM.qP = qP;                  % conditional moments of model-parameters
DEM.qH = qH;                  % conditional moments of hyper-parameters
 
DEM.F  = F;                   % [-ve] Free energy

warning on
