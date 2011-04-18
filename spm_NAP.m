function [DEM] = spm_NAP(DEM)
% Laplacian model inversion
% FORMAT DEM   = spm_NAP(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - response variable, output or data
% DEM.U  - explanatory variables, inputs or prior expectation of causes
%__________________________________________________________________________
%
% generative model
%--------------------------------------------------------------------------
%   M(i).g  = v     =  g(x,v,P)   {inline function, string or m-file}
%   M(i).f  = dx/dt =  f(x,v,P)   {inline function, string or m-file}
%
%   M(i).ph = pi(v) = ph(x,v,h,M) {inline function, string or m-file}
%   M(i).pg = pi(x) = pg(x,v,g,M) {inline function, string or m-file}
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h log-precision (cause noise)
%   M(i).hC = prior covariances of h log-precision (cause noise)
%   M(i).gE = prior expectation of g log-precision (state noise)
%   M(i).gC = prior covariances of g log-precision (state noise)
%   M(i).xP = precision (states)
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
% F         = log-evidence = log-marginal likelihood = negative free-energy
%__________________________________________________________________________
%
% spm_NAP implements a variational scheme under the Laplace
% approximation to the conditional joint density q on states (u), parameters 
% (p) and hyperparameters (h,g) of any analytic nonlinear hierarchical dynamic
% model, with additive Gaussian innovations.
%
%            q(u,p,h,g) = max <L(t)>q
%
% L is the ln p(y,u,p,h,g|M) under the model M. The conditional covariances
% obtain analytically from the curvature of L with respect to the unknowns.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_NAP.m 4310 2011-04-18 16:07:35Z guillaume $


% find or create a DEM figure
%--------------------------------------------------------------------------
Fdem = spm_figure('GetWin','DEM');

% check model, data and priors
%==========================================================================
[M Y U] = spm_DEM_set(DEM);


% number of iterations
%--------------------------------------------------------------------------
try nD = M(1).E.nD; catch nD = 1;   end
try nN = M(1).E.nN; catch nN = 16;  end


% ensure integration scheme evaluates gradients at each time-step
%--------------------------------------------------------------------------
M(1).E.linear = 4;

% assume predicions are a function of[f] hyperparameters
%--------------------------------------------------------------------------
try
    method = M(1).E.precision;
catch
    method = 1;
    M(1).E.precision = method;
end

% additional checks for Laplace models (precision functions; ph and pg)
%--------------------------------------------------------------------------
for i  = 1:length(M)
    try
        feval(M(i).ph,M(i).x,M(i).v,M(i),hE,M(i));
    catch
        M(i).ph = inline('spm_LAP_ph(x,v,h,M)','x','v','h','M');
    end
    try
        feval(M(i).pg,M(i).x,M(i).v,M(i),gE,M(i));
    catch
        M(i).pg = inline('spm_LAP_pg(x,v,h,M)','x','v','h','M');
    end
end

 
% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d   = M(1).E.d + 1;                          % embedding order of q(v)
n   = M(1).E.n + 1;                          % embedding order of q(x)
 
% number of states and parameters
%--------------------------------------------------------------------------
ns  = size(Y,2);                             % number of samples
nl  = size(M,2);                             % number of levels
nv  = sum(spm_vec(M.m));                     % number of v (casual states)
nx  = sum(spm_vec(M.n));                     % number of x (hidden states)
ny  = M(1).l;                                % number of y (inputs)
nc  = M(end).l;                              % number of c (prior causes)
nu  = nv*d + nx*n;                           % number of generalised states
ne  = nv*n + nx*n + ny*n;                    % number of generalised errors
 

% precision (R) of generalised errors and null matrices for concatenation
%==========================================================================
W   = sparse(nx*n,nx*n);
V   = sparse((ny + nv)*n,(ny + nv)*n);
 
% fixed priors on states (u)
%--------------------------------------------------------------------------
Px    = kron(sparse(1,1,1,n,n),spm_cat(spm_diag({M.xP})));
Pv    = kron(sparse(1,1,1,d,d),sparse(nv,nv));
pu.ic = spm_cat(spm_diag({Px Pv}));
 
% hyperpriors
%--------------------------------------------------------------------------
s     = M(1).E.s;
sh    = log(s);
sg    = log(s);
sC    = speye(2,2)/128;
ph.h  = spm_vec({M.hE M.gE sh sg});          % prior expectation of h,g
ph.c  = spm_cat(spm_diag({M.hC M.gC sC}));   % prior covariances of h,g
ph.ic = spm_pinv(ph.c);                      % prior precision of h,g
 
qh.h  = {M.hE};                              % conditional expectation h
qh.g  = {M.gE};                              % conditional expectation g
nh    = length(spm_vec(qh.h));               % number of hyperparameters h
ng    = length(spm_vec(qh.g));               % number of hyperparameters g
qh.sh = sh;                                  % conditional expectation sh
qh.sg = sg;                                  % conditional expectation sg
npp   = nh + ng + 2;                         % number of hyerparameters


% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
 
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC);                    % basis for parameters
    M(i).p    = size(qp.u{i},2);                     % number of qp.p
    qp.p{i}   = sparse(M(i).p,1);                    % initial deviates
    pp.c{i,i} = qp.u{i}'*M(i).pC*qp.u{i};            % prior covariance
 
end
Up    = spm_cat(spm_diag(qp.u));
 
% priors on parameters
%--------------------------------------------------------------------------
pp.p  = spm_vec(M.pE);
pp.c  = spm_cat(pp.c);
pp.ic = spm_inv(pp.c);
 
% initialise conditional density q(p)
%--------------------------------------------------------------------------
for i = 1:(nl - 1)
    try
        qp.p{i} = qp.p{i} + qp.u{i}'*(spm_vec(M(i).P) - spm_vec(M(i).pE));
    end
end
np    = size(Up,2);


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
Dx     = kron(spm_speye(n,n,1),spm_speye(nx,nx));
Dv     = kron(spm_speye(d,d,1),spm_speye(nv,nv));
Dy     = kron(spm_speye(n,n,1),spm_speye(ny,ny));
Dc     = kron(spm_speye(d,d,1),spm_speye(nc,nc));
Du     = spm_cat(spm_diag({Dx,Dv}));
Ip     = spm_speye(np, np );
Ih     = spm_speye(npp,npp);
qp.dp  = sparse(np,1);                   % conditional expectation of dp/dt
qh.dp  = sparse(npp,1);                  % conditional expectation of dh/dt


% gradients of generalised weighted errors
%--------------------------------------------------------------------------
dedh   = sparse(nh,ne);
dedg   = sparse(ng,ne);
dedv   = sparse(nv,ne);
dedx   = sparse(nx,ne);
            
% curvatures of Gibb's energy w.r.t. hyperparameters
%--------------------------------------------------------------------------
dedhh  = sparse(nh,nh);
dedgg  = sparse(ng,ng);
dedhsh = sparse(nh,1);
dedshh = sparse(1,nh);
dedgsg = sparse(ng,1);
dedsgg = sparse(1,ng);
dEdup  = sparse(nu,np);

% preclude unnecessary iterations
%--------------------------------------------------------------------------
if ~np && ~npp, nN = 1; end
if ~npp, method = 0;    end
 
% precision on parameter fluctuations
%--------------------------------------------------------------------------
kh  = ns*16;
kp  = ns*16;
 
% Iterate Lapalace scheme
%==========================================================================
Fs     = -Inf;
for iN = 1:nN
 
    % get time and clear persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic; clear spm_DEM_eval
 
    % [re-]set states & their derivatives
    %----------------------------------------------------------------------
    try, qu = Q(1).u; end
    
    % increase precision on parameter fluctuations
    %----------------------------------------------------------------------
    kp = kp + 64;
    kh = kh + 64;

 
    % D-Step: (nD D-Steps for each sample)
    %======================================================================
    for is = 1:ns
 
        % D-Step: until convergence for static systems
        %==================================================================
        for iD = 1:nD
 
            % sampling time
            %--------------------------------------------------------------
            ts = is + (iD - 1)/nD;
 
            % derivatives of responses and inputs
            %--------------------------------------------------------------
            try
                qu.y(1:n) = spm_DEM_embed(Y,n,ts,1,M(1).delays);
                qu.u(1:d) = spm_DEM_embed(U,d,ts);
            catch
                qu.y(1:n) = spm_DEM_embed(Y,n,ts);
                qu.u(1:d) = spm_DEM_embed(U,d,ts);
            end
            
            
            % evaluate functions and derivatives
            %==============================================================
            
            % prediction errors (E) and precision (p)
            %--------------------------------------------------------------
            [E dE]  = spm_DEM_eval(M,qu,qp);
            [p dp]  = spm_LAP_eval(M,qu,qh);
            
 
            % gradients of log(det(iS)) dDd...
            %==============================================================
            
            % get precision matrices
            %--------------------------------------------------------------
            [Rh Vh] = spm_DEM_R(n,exp(qh.sh));
            [Rg Vg] = spm_DEM_R(n,exp(qh.sg));
            iSh     = diag(exp(p.h));
            iSg     = diag(exp(p.g));
            iS      = blkdiag(kron(Rh,iSh),kron(Rg,iSg));
            
            
            % gradients of trace(diag(p)) = sum(p); p = precision vector
            %--------------------------------------------------------------
            dpdx    = n*sum(spm_cat({dp.h.dx; dp.g.dx}));
            dpdv    = n*sum(spm_cat({dp.h.dv; dp.g.dv}));
            dpdh    = n*sum(dp.h.dh);
            dpdg    = n*sum(dp.g.dg);
            dpdx    = kron(sparse(1,1,1,1,n),dpdx);
            dpdv    = kron(sparse(1,1,1,1,d),dpdv);
            dDdu    = [dpdx dpdv]';
            

            % gradients w.r.t. hyperparameters
            %--------------------------------------------------------------
            [dRdhh dRdh] = spm_diff('spm_DEM_R',n,exp(qh.sh),[2 2]);
            [dRdgg dRdg] = spm_diff('spm_DEM_R',n,exp(qh.sg),[2 2]);
            dRdh    = dRdh{1}*exp(qh.sh);
            dRdg    = dRdg{1}*exp(qh.sg);
            dRdhh   = dRdhh{1}{1}*exp(qh.sh) + dRdh;
            dRdgg   = dRdgg{1}{1}*exp(qh.sg) + dRdg;
            dDdh    = length(iSh)*trace(dRdh*Vh);
            dDdg    = length(iSg)*trace(dRdg*Vg);
            dDdh    = [dpdh dpdg dDdh dDdg]';
 
            
            % gradients precision-weighted generalised error dSd..
            %==============================================================
            diS     = blkdiag(kron(dRdh, iSh),W);
            dedsh   = E'*diS;
            diS     = blkdiag(kron(dRdhh,iSh),W);
            dedshsh = E'*diS*E;
            diS     = blkdiag(V,kron(dRdg, iSg));
            dedsg   = E'*diS;
            diS     = blkdiag(V,kron(dRdgg,iSg));
            dedsgsg = E'*diS*E;
            
            % gradients w.r.t. hyperparameters
            %--------------------------------------------------------------
            if method > 0
                for i = 1:nh
                    diS       = diag(dp.h.dh(:,i).*exp(p.h));
                    diS       = blkdiag(kron(Rh,diS),W);
                    dedh(i,:) = E'*diS;
                end
                for i = 1:ng
                    diS       = diag(dp.g.dg(:,i).*exp(p.g));
                    diS       = blkdiag(V,kron(Rg,diS));
                    dedg(i,:) = E'*diS;
                end
            end
            
            % gradients w.r.t. causal states
            %--------------------------------------------------------------
            if method > 1
                for i = 1:nv
                    diV       = diag(dp.h.dv(:,i).*exp(p.h));
                    diW       = diag(dp.g.dv(:,i).*exp(p.g));
                    diS       = blkdiag(kron(R,diV),kron(R,diW));
                    dedv(i,:) = E'*diS;
                end
            end
            
            % gradients w.r.t. hidden states
            %--------------------------------------------------------------
            if method > 2
                for i = 1:nx
                    diV       = diag(dp.h.dx(:,i).*exp(p.h));
                    diW       = diag(dp.g.dx(:,i).*exp(p.g));
                    diS       = blkdiag(kron(R,diV),kron(R,diW));
                    dedx(i,:) = E'*diS;
                end
            end

            dSdx  = kron(sparse(1,1,1,n,1),dedx);
            dSdv  = kron(sparse(1,1,1,d,1),dedv);
            dSdu  = [dSdx; dSdv];
            dEdh  = [dedh; dedg; dedsh; dedsg];
            dEdp  = dE.dp'*iS;
            dEdu  = dE.du'*iS;

            % curvatures w.r.t. hyperparameters
            %--------------------------------------------------------------
            for i = 1:nh
                for j = i:nh
                    diS        = diag(dp.h.dh(:,i).*dp.h.dh(:,j).*exp(p.h));
                    diS        = blkdiag(kron(Rh,diS),W);
                    dedhh(i,j) = E'*diS*E;
                    dedhh(j,i) = dedhh(i,j);
                end
                diS         = diag(dp.h.dh(:,i).*exp(p.h));
                diS         = blkdiag(kron(dRdh,diS),W);
                dedhsh(i,1) = E'*diS*E;
                dedshh(1,i) = dedhsh(i,1);
            end            
            for i = 1:ng
                for j = i:ng
                    diS        = diag(dp.g.dg(:,i).*dp.g.dg(:,j).*exp(p.g));
                    diS        = blkdiag(V,kron(Rg,diS));
                    dedgg(i,j) = E'*diS*E;
                    dedgg(j,i) = dedgg(i,j);
                end
                diS         = diag(dp.g.dg(:,i).*exp(p.g));
                diS         = blkdiag(V,kron(dRdg,diS));
                dedgsg(i,1) = E'*diS*E;
                dedsgg(1,i) = dedgsg(i,1);
            end
            
            % combined curvature
            %--------------------------------------------------------------
            dSdhh = spm_cat({dedhh  []     dedhsh  [];
                             []     dedgg  []      dedgsg;
                             dedshh []     dedshsh [];
                             []     dedsgg []      dedsgsg});
                 
            
            % curvatures w.r.t. parameters and states
            %--------------------------------------------------------------
            % for i = 1:np, dEdup(:,i) = dE.dup{i}'*iS*E; end
            
            % errors (from prior expectations) (NB pp.p = 0)
            %--------------------------------------------------------------
            Eu    = spm_vec(qu.x(1:n),qu.v(1:d));
            Ep    = spm_vec(qp.p);
            Eh    = spm_vec(qh.h,qh.g,qh.sh,qh.sg) - ph.h;
            
 
            % first-order derivatives of Gibb's Energy
            %==============================================================
            dLdu  = dEdu*E + dSdu*E/2 - dDdu/2 + pu.ic*Eu;
            dLdh  = dEdh*E/2          - dDdh/2 + ph.ic*Eh;
            dLdp  = dEdp*E                     + pp.ic*Ep;
            
 
            % and second-order derivatives of Gibb's Energy
            %--------------------------------------------------------------
            dLduu = dEdu*dE.du + dSdu*dE.du*2 + pu.ic;
            dLdup = dEdu*dE.dp + dSdu*dE.dp + dEdup;
            dLdpp = dEdp*dE.dp + pp.ic;
            dLdhh = dSdhh/2    + ph.ic;            
            dLdhu = dEdh*dE.du;
            dLduy = dEdu*dE.dy;
            dLduc = dEdu*dE.dc;
            dLdpy = dEdp*dE.dy;
            dLdpc = dEdp*dE.dc;
            dLdhy = dEdh*dE.dy;
            dLdhc = dEdh*dE.dc;
            dLdhp = dEdh*dE.dp;
            dLdpu = dLdup';
            dLduh = dLdhu';
            dLdph = dLdhp';
            
 
            % save conditional moments (and prediction error) at Q{t}
            %==============================================================
            if iD == 1
                
                % means
                %----------------------------------------------------------
                Q(is).e = E;
                Q(is).u = qu;
                Q(is).p = qp;
                Q(is).h = qh;
                
                % precision and covariances
                %----------------------------------------------------------                        
                iC = spm_cat({dLduu dLdup  [];
                              dLdpu dLdpp  [];
                              []    [] dLdhh});
                
                C  = spm_inv(iC);
                                
                % save conditional covariances
                %----------------------------------------------------------
                Q(is).u.s = C([1:nx],[1:nx]);
                Q(is).u.c = C([1:nv]  + nx*n,   [1:nv]  + nx*n);
                Q(is).p.c = C([1:np]  + nu,     [1:np]  + nu);
                Q(is).h.c = C([1:npp] + nu + np,[1:npp] + nu + np);

                % Free-energy (states)
                %----------------------------------------------------------                
                L(is) = ... 
                - E'*iS*E/2      + spm_logdet(iS)/2    - n*ny*log(2*pi)/2 ...          
                - Eu'*pu.ic*Eu/2 + spm_logdet(pu.ic)/2 - spm_logdet(dLduu)/2;
                    
                % Free-energy (parameters)
                %----------------------------------------------------------
                A(is) = - E'*iS*E/2        + spm_logdet(iS)/2    ...
                        - Eu'*pu.ic*Eu/2   + spm_logdet(pu.ic)/2 ...
                        - Ep'*pp.ic*Ep/2   + spm_logdet(pp.ic)/2 ...
                        - Eh'*ph.ic*Eh/2   + spm_logdet(ph.ic)/2 ...
                        - n*ny*log(2*pi)/2 + spm_logdet(C)/2;

 
            end
 
            % update conditional moments
            %==============================================================
            
            % precision of fluctuations
            %--------------------------------------------------------------
            Kp    = kp*Ip;
            Kh    = kh*Ih;
            
            % assemble conditional means
            %--------------------------------------------------------------
            q{1}  = qu.y(1:n);
            q{2}  = qu.x(1:n);
            q{3}  = qu.v(1:d);
            q{4}  = qu.u(1:d);
            q{5}  = qp.p;
            q{6}  = qh.h;
            q{7}  = qh.g;
            q{8}  = qh.sh;
            q{9}  = qh.sg;
            q{10} = qp.dp;
            q{11} = qh.dp;

            % flow
            %--------------------------------------------------------------
            f{1}  =  Dy*spm_vec(q{1});
            f{2}  =  Du*spm_vec(q{2:3}) - dLdu;
            f{3}  =  Dc*spm_vec(q{4});
            f{4}  =     spm_vec(q{10});
            f{5}  =     spm_vec(q{11});
            f{6}  = -kp*spm_vec(q{10})  - dLdp;
            f{7}  = -kh*spm_vec(q{11})  - dLdh;
 
            % and Jacobian
            %--------------------------------------------------------------
            dfdq  = spm_cat({Dy      []       []     []     []     []   [];
                            -dLduy  Du-dLduu -dLduc  []     []     []   [];
                             []      []       Dc     []     []     []   [];
                             []      []       []     []     []     Ip   [];
                             []      []       []     []     []     []   Ih;
                            -dLdpy  -dLdpu   -dLdpc -dLdpp -dLdph -Kp   [];
                            -dLdhy  -dLdhu   -dLdhc -dLdhp -dLdhh  []  -Kh});
 
 
            % update conditional modes of states
            %==============================================================
            dq    = spm_dx(dfdq, spm_vec(f), 1/nD);
            q     = spm_unvec(spm_vec(q) + dq,q);
            
            % unpack conditional means
            %--------------------------------------------------------------
            qu.x(1:n) = q{2};
            qu.v(1:d) = q{3};
            qp.p      = q{5};
            qh.h      = q{6};
            qh.g      = q{7};
            qp.sh     = q{8};
            qh.sg     = q{9};
            qp.dp     = q{10};
            qh.dp     = q{11};

 
        end % D-Step
 
    end % sequence (ns)
 
    
    % Bayesian parameter averaging
    %======================================================================

    % Conditional moments of time-averaged parameters
    %----------------------------------------------------------------------
    Pp  = 0;
    Ep  = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).p.c);
        Ep  = Ep + P*spm_vec(Q(i).p.p);
        Pp  = Pp + P;       
    end
    Cp  = spm_inv(Pp);
    Ep  = Cp*Ep;

    % conditional moments of hyper-parameters
    %----------------------------------------------------------------------
    Ph  = 0;
    Eh  = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).h.c);
        Ph  = Ph + P;
        Eh  = Eh + P*spm_vec({Q(i).h.h Q(i).h.g Q(i).h.sh Q(i).h.sg});
    end
    Ch  = spm_inv(Ph);
    Eh  = Ch*Eh - ph.h;

    % Free-action of states plus free-energy of parameters
    %======================================================================
    Fs  = sum(A);
    Fi  = sum(L) ...
          - Ep'*pp.ic*Ep/2 + spm_logdet(pp.ic)/2 - spm_logdet(Pp)/2 ...
          - Eh'*ph.ic*Eh/2 + spm_logdet(ph.ic)/2 - spm_logdet(Ph)/2;


    % if F is increasing terminate
    %----------------------------------------------------------------------
    if Fi < Fa && iN > 16
        break
    else
        Fa    = Fi;
        F(iN) = Fi;
        S(iN) = Fs;
    end
 
    % otherwise save conditional moments (for each time point)
    %======================================================================
    for t = 1:length(Q)
 
 
        % states and predictions
        %------------------------------------------------------------------
        v     = spm_unvec(Q(t).u.v{1},v);
        x     = spm_unvec(Q(t).u.x{1},x);
        z     = spm_unvec(Q(t).e(1:(ny + nv)),{M.v});
        w     = spm_unvec(Q(t).e([1:nx] + (ny + nv)*n),{M.x});
        for i = 1:(nl - 1)
            if M(i).m, qU.v{i + 1}(:,t) = spm_vec(v{i});  end
            if M(i).n, qU.x{i}(:,t)     = spm_vec(x{i});  end
            if M(i).n, qU.w{i}(:,t)     = spm_vec(w{i});  end
            if M(i).l, qU.z{i}(:,t)     = spm_vec(z{i});  end
        end
        if    M(nl).l, qU.z{nl}(:,t)    = spm_vec(z{nl}); end
        qU.v{1}(:,t)  = spm_vec(Q(t).u.y{1}) - spm_vec(z{1});
 
        % and conditional covariances
        %------------------------------------------------------------------
        qU.S{t} = Q(t).u.s;
        qU.C{t} = Q(t).u.c;
 
        % parameters
        %------------------------------------------------------------------
        qP.p{t} = spm_vec(Q(t).p.p);
        qP.c{t} = Q(t).p.c;
 
        % hyperparameters
        %------------------------------------------------------------------
        qH.p{t} = spm_vec({Q(t).h.h Q(t).h.g Q(t).h.sh Q(t).h.sg});
        qH.c{t} = Q(t).h.c;
 
    end
 
    % graphics (states)
    %----------------------------------------------------------------------
    figure(Fdem)
    spm_DEM_qU(qU)
    
    % determine graphs to plot
    %----------------------------------------------------------------------
    if np && nh
        ap = subplot(2*nl,2,4*nl - 2);
        ah = subplot(2*nl,2,4*nl);
    elseif np
        ap = subplot(nl,2,2*nl);
    elseif nh
        ah = subplot(nl,2,2*nl);
    end
    
    % graphics (parameters and log-precisions)
    %----------------------------------------------------------------------
    if np
        axes(ap)
        plot([1:ns],spm_cat(qP.p))
        set(gca,'XLim',[1 ns])
        title('parameters','FontSize',16)
    end
    if nh
        axes(ah)
        plot([1:ns],spm_cat(qH.p))
        set(gca,'XLim',[1 ns])
        title('log-precision','FontSize',16)
    end
    drawnow
 
    % report (EM-Steps)
    %----------------------------------------------------------------------
    try
        dF = F(end) - F(end - 1);
    catch
        dF = 0;
    end
    str{1} = sprintf('LAP: %i (%i)', iN,iD);
    str{2} = sprintf('F:%.4e',       full(Fs - F(1)));
    str{3} = sprintf('dF:%.2e',      full(dF));
    str{4} = sprintf('(%.2e sec)',   full(toc));
    fprintf('%-16s%-16s%-14s%-16s\n',str{:})
 
end
 
 
% Place Bayesian parameter averages in output arguments
%==========================================================================
 
% Conditional moments of time-averaged parameters
%--------------------------------------------------------------------------
Pp = 0;
Ep = 0;
for i = 1:ns
    
    % weight in proportion to precisions
    %----------------------------------------------------------------------
    P  = spm_inv(qP.c{i});
    Ep = Ep + P*qP.p{i};
    Pp = Pp + P;
 
end
Cp     = spm_inv(Pp);
Ep     = Cp*Ep;
P      = {M.pE};
qP.P   = spm_unvec(Up*Ep + pp.p,P);
qP.C   = Up*Cp*Up';
qP.V   = spm_unvec(diag(qP.C),P);
qP.U   = Up;
 
% conditional moments of hyper-parameters
%--------------------------------------------------------------------------
Ph = 0;
Eh = 0;
for i = 1:ns
    
    % weight in proportion to precisions
    %----------------------------------------------------------------------
    P  = spm_inv(qH.c{i});
    Ph = Ph + P;
    Eh = Eh + P*qH.p{i};
 
end
Ch     = spm_inv(Ph);
Eh     = Ch*Eh;
P      = {qh.h qh.g qh.sh qh.sg};
P      = spm_unvec(Eh,P);
qH.h   = P{1};
qH.g   = P{2};
qH.sh  = P{3};
qH.sg  = P{4};
qH.C   = Ch;
P      = spm_unvec(diag(qH.C),P);
qH.V   = P{1};
qH.W   = P{2};

 
 
% assign output variables
%--------------------------------------------------------------------------
DEM.M  = M;                   % model
DEM.U  = U;                   % causes
 
DEM.qU = qU;                  % conditional moments of model-states
DEM.qP = qP;                  % conditional moments of model-parameters
DEM.qH = qH;                  % conditional moments of hyper-parameters
 
DEM.F  = F;                   % [-ve] Free energy
DEM.S  = S;                   % [-ve] Free action

return



% Notes
%==========================================================================

                % numerical approximations
                %----------------------------------------------------------
                qq.x  = qu.x(1:n);
                qq.v  = qu.v(1:d);
                qq.p  = qp.p;
                qq.h  = qh.h;
                qq.g  = qh.g;
                qq.sh = qh.sh;
                qq.sg = qh.sg;

                dLdqq = spm_diff('spm_LAP_F',qq,qu,qp,qh,pu,pp,ph,M,[1 1]);
                dLdqq = spm_cat(dLdqq');
                
                subplot(2,2,1);imagesc(dLdqq);     axis square
                subplot(2,2,2);imagesc(iC);        axis square
                subplot(2,2,3);imagesc(dLdqq - iC);axis square
                subplot(2,2,4);plot(iC,':k');hold on;
                plot(dLdqq - iC,'r');hold off; axis square
                drawnow
