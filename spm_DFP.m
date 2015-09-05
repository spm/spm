function [DEM] = spm_DFP(DEM)
% Dynamic free-energy Fokker-Planck free-form scheme
% FORMAT [DEM] = spm_DFP(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - output or data
% DEM.U  - inputs or prior expectation of causes
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
%   M(i).hE = prior expectation of h hyper-parameters (input noise)
%   M(i).hC = prior covariances of h hyper-parameters (input noise)
%   M(i).gE = prior expectation of g hyper-parameters (state noise)
%   M(i).gC = prior covariances of g hyper-parameters (state noise)
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%
% conditional moments of model-states - q(u)
%--------------------------------------------------------------------------
%   qU.x    = Conditional expectation of hidden states
%   qU.v    = Conditional expectation of causal states
%   qU.e    = Conditional residuals
%   qU.C    = Conditional covariance: cov(v)
%   qU.S    = Conditional covariance: cov(x)
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
% spm_DFP implements a variational Bayes (VB) scheme under the Laplace
% approximation to the conditional densities of the model's, parameters (p)
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
% The D-step is implemented with variational filtering, which does not
% assume a fixed form for the conditional density; it uses the sample
% density of an ensemble of particles that drift up free-energy gradients
% and 'explore' the local curvature though (Wiener) perturbations.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DFP.m 6540 2015-09-05 10:06:42Z karl $


% Check model, data, priros and confounds and unpack
%--------------------------------------------------------------------------
[M,Y,U,X] = spm_DEM_set(DEM);
MOVIE     = 0;

% find or create a DEM figure
%--------------------------------------------------------------------------
Fdem = spm_figure('GetWin','DEM');
Fdfp = spm_figure('GetWin','DFP');

% tolerance for changes in norm
%--------------------------------------------------------------------------
TOL  = 1e-2;

% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d    = M(1).E.d + 1;                   % embedding order of q(v)
n    = M(1).E.n + 1;                   % embedding order of q(x)
s    = M(1).E.s;                       % smoothness - s.d. of kernel (bins)
try
    N = M(1).E.N;                      % number of particles
catch
    N = 16;                            % number of particles
end

% number of states and parameters
%--------------------------------------------------------------------------
nY   = size(Y,2);                      % number of samples
nl   = size(M,2);                      % number of levels
ne   = sum(cat(1,M.l));                % number of e (errors)
nv   = sum(cat(1,M.m));                % number of v (casual states)
nx   = sum(cat(1,M.n));                % number of x (hidden states)
ny   = M(1).l;                         % number of y (inputs)
nc   = M(end).l;                       % number of c (prior causes)
nu   = nv*d + nx*n;                    % number of generalised states
kt   = 1;                              % rate constant for D-Step

% number of iterations
%--------------------------------------------------------------------------
if nx, nD = 1;       else   nD = 8; end
try, nE = M(1).E.nE; catch, nE = 1; end
try, nM = M(1).E.nM; catch, nM = 8; end
try, nN = M(1).E.nN; catch, nN = 8; end


% initialise regularisation parameters
%--------------------------------------------------------------------------
td    = 1/nD;                          % integration time for D-Step
te    = 2;                             % integration time for E-Step

%  Precision (R) and covariance of generalised errors
%--------------------------------------------------------------------------
[iV,V] = spm_DEM_R(n,s);

% precision components Q{} requiring [Re]ML estimators (M-Step)
%==========================================================================
Q     = {};
for i = 1:nl
    q0{i,i} = sparse(M(i).l,M(i).l);
end
for i = 1:nl - 1
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
end
for i = 1:nl - 1
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


% hyperpriors
%--------------------------------------------------------------------------
ph.h  = spm_vec({M.hE; M.gE});             % prior expectation of h
ph.c  = spm_cat(spm_diag({M.hC M.gC}));        % prior covariances of h
ph.ic = spm_pinv(ph.c);                    % prior precision
qh.h  = ph.h;                              % conditional expectation
qh.c  = ph.c;                              % conditional covariance
 

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

end
Up    = spm_cat(spm_diag(qp.u));
 
% initialise and augment with confound parameters B; with flat priors
%--------------------------------------------------------------------------
np    = sum(cat(1,M.p));                    % number of model parameters
nb    = size(X,1);                          % number of confounds
nn    = nb*ny;                              % number of nuisance parameters
nf    = np + nn;                            % numer of free parameters
ip    = [1:np];
ib    = [1:nn] + np;
pp.c  = spm_cat(pp.c);
pp.ic = spm_pinv(pp.c);
 
% initialise conditional density q(p) (for D-Step)
%--------------------------------------------------------------------------
qp.e  = spm_vec(qp.p);
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
qu.x      = cell(n + 1,1);
qu.v      = cell(n + 1,1);
qy        = cell(n + 1,1);
qc        = cell(n + 1,1);
[qu.x{:}] = deal(sparse(nx,1));
[qu.v{:}] = deal(sparse(nv,1));
[qy{:}  ] = deal(sparse(ny,1));
[qc{:}  ] = deal(sparse(nc,1));

% initialise cell arrays for hierarchical structure of x[0] and v[0]
%--------------------------------------------------------------------------
x         = {M(1:end - 1).x};
v         = {M(1 + 1:end).v};
qu.x{1}   = spm_vec(x);
qu.v{1}   = spm_vec(v);
qu(1:N)   = deal(qu);

try xp = DEM.M(1).E.xp; catch, xp = 1; end
try vp = DEM.M(1).E.vp; catch, vp = 1; end
for i  = 1:N
    qu(i).x{1} = qu(i).x{1} + randn(nx,1)/xp;
    qu(i).v{1} = qu(i).v{1} + randn(nv,1)/vp;
end

dq        = {qu(1).x{1:n} qu(1).v{1:d} qy{1:n} qc{1:d}};

% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dx        = cell(n,n);
Dv        = cell(d,d);
Dy        = cell(n,n);
Dc        = cell(d,d);
[Dx{:}]   = deal(sparse(nx,nx));
[Dv{:}]   = deal(sparse(nv,nv));
[Dy{:}]   = deal(sparse(ny,ny));
[Dc{:}]   = deal(sparse(nc,nc));

% Wiener process
%--------------------------------------------------------------------------
Ix        = Dx;
Iv        = Dv;
for i = 1:d
    Iv{i,i} = speye(nv,nv);
end
for i = 1:n
    Ix{i,i} = speye(nx,nx);
end
dfdw        = spm_cat(spm_diag({Ix,Iv,Dy,Dc}));

% add constant terms
%--------------------------------------------------------------------------
for i = 2:d
    Dv{i - 1,i} = speye(nv,nv);
    Dc{i - 1,i} = speye(nc,nc);
end
for i = 2:n
    Dx{i - 1,i} = speye(nx,nx);
    Dy{i - 1,i} = speye(ny,ny);
end
Du        = spm_cat(spm_diag({Dx,Dv}));
Dc        = spm_cat(Dc);
Dy        = spm_cat(Dy);

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
if ~nf && ~nh, nN = 1; end


% Iterate DEM
%==========================================================================
for iN = 1:nN
    
    % get time and celar persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic; clear spm_DEM_eval
 
    % E-Step: (with embedded D-Step)
    %======================================================================
    mp     = zeros(nf,1);
    for iE = 1:nE
 

        % [re-]set accumulators for E-Step
        %------------------------------------------------------------------
        dFdp  = zeros(nf,1);
        dFdpp = zeros(nf,nf);
        EE    = sparse(0);
        ECE   = sparse(0);
        qp.ic = sparse(0);
        qu_c  = speye(1);
 

        % [re-]set precisions using ReML hyperparameter estimates
        %------------------------------------------------------------------
        iS    = Qp;
        for i = 1:nh
           iS = iS + Q{i}*exp(qh.h(i));
        end
        
        % [re-]adjust for confounds
        %------------------------------------------------------------------
        Y     = Y - qp.b*X;
        
        % [re-]set states & their derivatives
        %------------------------------------------------------------------
        try
            qu = QU{1};
        end
        
        % D-Step: (nD D-Steps for each sample)
        %==================================================================
        for iY = 1:nY

            % [re-]set states for static systems
            %--------------------------------------------------------------
            if ~nx
                try, qu = QU{iY}; end
            end

            % D-Step: until convergence for static systems
            %==============================================================
            for iD = 1:nD
                
                % sampling time
                %----------------------------------------------------------
                ts      = iY + (iD - 1)/nD;
                
                % derivatives of responses and inputs
                %----------------------------------------------------------
                qy(1:n) = spm_DEM_embed(Y,n,ts);
                qc(1:d) = spm_DEM_embed(U,d,ts);

                % compute dEdb (derivatives of confounds)
                %----------------------------------------------------------
                b     = spm_DEM_embed(X,n,ts);
                for i = 1:n
                    dedbi{1}  = -kron(b{i}',speye(ny,ny));
                    dEdb{i,1} =  spm_cat(dedbi);
                end

                % moments of ensemble density
                %==========================================================
                q     = [qu.x];
                for i = 1:n + 1
                       qx{i} = mean([q{i,:}],2);
                end
                q     = [qu.v];
                for i = 1:n + 1
                       qv{i} = mean([q{i,:}],2);
                end
                
                % mean field effects
                %----------------------------------------------------------
                dudt   = spm_vec({qx(2:n + 1)
                                  qv(2:d + 1)
                                  qy(2:n + 1)
                                  qc(2:d + 1)});
                
                if iD == nD

                    % ensemble covariance
                    %------------------------------------------------------
                    ux     = [qu.x];
                    uv     = [qu.v];
                    ux     = ux(1:n,:);
                    uv     = uv(1:d,:);
                    c      = cov(spm_cat([ux; uv])');
                    
                    % quantities for M-Step
                    %------------------------------------------------------
                    quy.x  = qx;
                    quy.v  = qv;
                    quy.y  = qy;
                    quy.u  = qc;
                    [E,dE] = spm_DEM_eval(M,quy,qp);
                    dE.dP  = [dE.dp spm_cat(dEdb)];

                    qu_c   = qu_c*c;
                    EE     = E*E'+ EE;
                    ECE    = ECE + dE.du*c*dE.du'+ dE.dP*qp.c*dE.dP';

                    % save states for qu(iY)
                    %------------------------------------------------------
                    qU(iY).x = qx;
                    qU(iY).v = qv;
                    qU(iY).y = qy;
                    qU(iY).u = qc;
                    qU(iY).e = E;
                    qU(iY).c = c;
                    
                end

                % evaluate functions:
                % e = v - g(x,v), dx/dt = f(x,v) and derivatives dE.dx, ...
                %==========================================================
                for iP = 1:N
                    
                    quy.x  = qu(iP).x;
                    quy.v  = qu(iP).v;
                    quy.y  = qy;
                    quy.u  = qc;
                    [e,de] = spm_DEM_eval(M,quy,qp);

                    % conditional uncertainty about parameters
                    %======================================================
                    if np
                        for i = 1:nu

                            % 1st-order derivatives: dUdv, ... ;
                            %----------------------------------------------
                            CJ             = qp.c(ip,ip)*de.dpu{i}'*iS;
                            dUdu(i,1)      = trace(CJ*de.dp);

                            % 2nd-order derivatives
                            %----------------------------------------------
                            for j = 1:nu
                                dUduu(i,j) = trace(CJ*de.dpu{j});
                            end
                        end
                    end


                    % D-step update: of causes v{i}, and other states u(i)
                    %======================================================

                    % compute dqdt: q = {u y c}; and dudt: u = {v{1:d} x}
                    %------------------------------------------------------
                    dIdu  = -de.du'*iS*e     - dUdu/2;

                    % and second-order derivatives
                    %------------------------------------------------------
                    dIduu = -de.du'*iS*de.du - dUduu/2;
                    dIduy = -de.du'*iS*de.dy;
                    dIduc = -de.du'*iS*de.dc;

                    % gradient
                    %------------------------------------------------------
                    dFdu       = dudt;
                    dFdu(1:nu) = dIdu + dFdu(1:nu);

                    % Jacobian
                    %------------------------------------------------------
                    dFduu = spm_cat({dIduu dIduy dIduc;
                                     []    Dy    []   ;
                                     []    []    Dc}) ;


                    % update conditional modes of states
                    %------------------------------------------------------
                    du    = spm_sde_dx(dFduu,dfdw,dFdu,td);
                    dq    = spm_unvec(du,dq);
                    for i = 1:n
                        qu(iP).x{i} = qu(iP).x{i} + dq{i};
                    end
                    for i = 1:d
                        qu(iP).v{i} = qu(iP).v{i} + dq{i + n};
                    end

                end

            end % D-Step
            
            % D-Step: save ensemble density and plot (over samples)
            %--------------------------------------------------------------
            QU{iY} = qu;
            figure(Fdfp)
            spm_DFP_plot(QU,nY)
            if MOVIE
                subplot(2,1,1)
                set(gca,'YLim',[-0.4 1.2])
                drawnow
                MOV(iY) = getframe(gca); 
            end

            % Gradients and curvatures for E-Step:
            %==============================================================
            for i = ip
                
                % 1st-order derivatives: U = tr(C*J'*iS*J)
                %----------------------------------------------------------
                CJ             = c*dE.dup{i}'*iS;
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
            
        end % sequence (iY)
 
        % augment with priors
        %------------------------------------------------------------------
        dFdp(ip)     = dFdp(ip)     - pp.ic*qp.e;
        dFdpp(ip,ip) = dFdpp(ip,ip) - pp.ic;
        qp.ic(ip,ip) = qp.ic(ip,ip) + pp.ic;
        qp.c         = spm_pinv(qp.ic);
        
        
        % E-step: update expectation (p)
        %==================================================================

        % update conditional expectation
        %------------------------------------------------------------------
        dp   = spm_dx(dFdpp,dFdp,{te});
        qp.e = qp.e + dp(ip);
        qp.p = spm_unvec(qp.e,qp.p);
        qp.b = spm_unvec(dp(ib),qp.b);
        mp   = mp + dp;

        % convergence (E-Step)
        %------------------------------------------------------------------
        if (dFdp'*dp < 1e-2) || (norm(dp,1) < TOL), break, end
        
    end % E-Step
    
    
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
        S     = inv(iS);
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
        qh.c = -spm_pinv(dFdhh);
        
        % convergence (M-Step)
        %------------------------------------------------------------------
        if (dFdh'*dh < 1e-2) || (norm(dh,1) < TOL), break, end
        
    end % M-Step

    % evaluate objective function (F)
    %======================================================================
    F(iN) = - trace(iS*EE)/2  ...                % states (u)
            - trace(qp.e'*pp.ic*qp.e)/2  ...     % parameters (p)
            - trace(qh.e'*ph.ic*qh.e)/2  ...     % hyperparameters (h)
            + spm_logdet(qu_c)/2  ...            % entropy q(u)
            + spm_logdet(qp.c)/2  ...            % entropy q(p)
            + spm_logdet(qh.c)/2  ...            % entropy q(h)
            - spm_logdet(pp.c)/2  ...            % entropy - prior p
            - spm_logdet(ph.c)/2  ...            % entropy - prior h
            + spm_logdet(iS)*nY/2 ...            % entropy - error
            - n*ny*nY*log(2*pi)/2;

    % save model-states (for each time point)
    %======================================================================
    for t = 1:length(qU)
        v     = spm_unvec(qU(t).v{1},v);
        x     = spm_unvec(qU(t).x{1},x);
        e     = spm_unvec(qU(t).e,{M.v});
        for i = 1:(nl - 1)
            Qu.v{i + 1}(:,t) = spm_vec(v{i});
            try
                Qu.x{i}(:,t) = spm_vec(x{i});
            end
            Qu.z{i}(:,t)     = spm_vec(e{i});
        end
        Qu.v{1}(:,t)         = spm_vec(qU(t).y{1} - e{1});
        Qu.z{nl}(:,t)        = spm_vec(e{nl});

        % and conditional covariances
        %--------------------------------------------------------------
        i       = [1:nx];
        Qu.S{t} = qU(t).c(i,i);
        i       = [1:nv] + nx*n;
        Qu.C{t} = qU(t).c(i,i);
    end

    % report and break if convergence
    %------------------------------------------------------------------
    figure(Fdem)
    spm_DEM_qU(Qu)
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
    str{1} = sprintf('DEM: %i (%i:%i:%i)',iN,iD,iE,iM);
    str{2} = sprintf('F:%.6e',full(F(iN)));
    str{3} = sprintf('p:%.2e',full(mp'*mp));
    str{4} = sprintf('h:%.2e',full(mh'*mh));
    str{5} = sprintf('(%.2e sec)',full(toc));
    fprintf('%-16s%-16s%-14s%-14s%-16s\n',str{:})
    
    if norm(mp) < TOL && norm(mh) < TOL, break, end


end

% Assemble output arguments
%==========================================================================

% conditional moments of model-parameters (rotated into original space)
%--------------------------------------------------------------------------
qP.P   = spm_unvec(Up*qp.e + spm_vec(M.pE),M.pE);
qP.C   = Up*qp.c(ip,ip)*Up';
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
DEM.X  = X;                   % confounds

DEM.QU = QU;                  % sample density of model-states
DEM.qU = Qu;                  % conditional moments of model-states
DEM.qP = qP;                  % conditional moments of model-parameters
DEM.qH = qH;                  % conditional moments of hyper-parameters
 
DEM.F  = F;                   % [-ve] Free energy

% set ButtonDownFcn
%--------------------------------------------------------------------------
if MOVIE
    figure(Fdfp), subplot(2,1,1)
    set(gca,'Userdata',{MOV,16})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
end
