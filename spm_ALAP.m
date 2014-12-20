function [DEM] = spm_ALAP(DEM)
% Laplacian model inversion (see also spm_LAP) with action
% FORMAT DEM   = spm_ALAP(DEM)
%
% DEM.G  - generative process
% DEM.M  - recognition  model
% DEM.C  - causes (n x t)
% DEM.U  - prior expectation of causes
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
%   pi(v,x) = vectors of log-precisions; (h,g) = precision parameters
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h log-precision (cause noise)
%   M(i).hC = prior covariances of h log-precision (cause noise)
%   M(i).gE = prior expectation of g log-precision (state noise)
%   M(i).gC = prior covariances of g log-precision (state noise)
%
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%   M(i).xP = precision (states)
%
%   M(i).m  = number of hidden inputs v(i + 1);
%   M(i).n  = number of hidden states x(i);
%   M(i).l  = number of outputs v(i);
%
% or (inital values)
%
%   M(i).x  = hidden states
%   M(i).v  = hidden causes
%
% hierarchical process G(i)
%--------------------------------------------------------------------------
%   G(i).g  = y(t)  = g(x,v,[a],P)    {inline function, string or m-file}
%   G(i).f  = dx/dt = f(x,v,[a],P)    {inline function, string or m-file}
%
%   G(i).pE = model-parameters
%   G(i).U  = precision (on sensory prediction errors - for action)
%   G(i).V  = precision (input noise)
%   G(i).W  = precision (state noise)
%
%   G(i).m  = number of inputs v(i + 1);
%   G(i).n  = number of states x(i)
%   G(i).l  = number of output v(i)
%   G(i).k  = number of action a(i)
%
% or (inital values)
%
%   G(i).x  = states
%   G(i).v  = causes
%   G(i).a  = action
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
%   F       = log-evidence = log-marginal likelihood = negative free-energy
%
%__________________________________________________________________________
% Accelerated methods: To accelerate computations one can specify the 
% nature of the model equations using:
%
% M(1).E.linear = 0: full        - evaluates 1st and 2nd derivatives
% M(1).E.linear = 1: linear      - equations are linear in x and v
% M(1).E.linear = 2: bilinear    - equations are linear in x, v & x*v
% M(1).E.linear = 3: nonlinear   - equations are linear in x, v, x*v, & x*x
% M(1).E.linear = 4: full linear - evaluates 1st derivatives (for GF)
%
% similarly, for evaluating precisions:
%
% M(1).E.method.h = 0,1  switch for precision parameters (hidden causes)
% M(1).E.method.g = 0,1  switch for precision parameters (hidden states)
% M(1).E.method.x = 0,1  switch for precision (hidden causes)
% M(1).E.method.v = 0,1  switch for precision (hidden states)
%__________________________________________________________________________
%
%__________________________________________________________________________
%
% spm_ALAP implements a variational scheme under the Laplace
% approximation to the conditional joint density q on states u, parameters
% p and hyperparameters (h,g) of an analytic nonlinear hierarchical dynamic
% model, with additive Gaussian innovations.
%
%            q(u,p,h,g) = max E[L(t)] - H(q(u,p,h,g))
%
% L is the ln p(y,u,p,h,g|M) under the model M. The conditional covariances
% obtain analytically from the curvature of L with respect to the unknowns.
%
% This implementation is the same as spm_LAP but integrates both the
% generative process and model inversion in parallel. Its functionality is
% exactly the same apart from the fact that confounds are not accommodated
% explicitly.  The generative model is specified by DEM.G and the veridical
% causes by DEM.C; these may or may not be used as priors on the causes for
% the inversion model DEM.M (i.e., DEM.U = DEM.C).  Clearly, DEM.G does not
% require any priors or precision components; it will use the values of the
% parameters specified in its prior expectation fields.
%
% This routine is not used for model inversion per se but to simulate the
% dynamical inversion of models.  Critically, it includes action
% variables a - that couple the model back to the generative process
% This enables active inference (c.f., action-perception) or embodied
% inference.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ALAP.m 6290 2014-12-20 22:11:50Z karl $


% check model, data and priors
%==========================================================================
DEM   = spm_ADEM_set(DEM);
M     = DEM.M;
G     = DEM.G;
C     = DEM.C;
U     = DEM.U;

% ensure embedding dimensions are compatible
%--------------------------------------------------------------------------
G(1).E.n = M(1).E.n;
G(1).E.d = M(1).E.n;


% set regularisation
%--------------------------------------------------------------------------
try
    dt = DEM.M(1).E.v;
catch
    dt = 0;
    DEM.M(1).E.v = dt;
end


% number of iterations of active inference
%--------------------------------------------------------------------------
try, nN = M(1).E.nN; catch, nN = 16;  end

% ensure integration scheme evaluates gradients at each time-step
%--------------------------------------------------------------------------
M(1).E.linear = 4;

% assume precisions have a Gaussian autocorrelation function
%--------------------------------------------------------------------------
try
    form = M(1).E.form;
catch
    form = 'Gaussian';
end

% checks for state-dependent precision (precision functions; ph and pg)
%--------------------------------------------------------------------------
for i  = 1:length(M)
    try
        feval(M(i).ph,M(i).x,M(i + 1).v,M(i).hE,M(i));
    catch
        M(i).ph = inline('spm_LAP_ph(x,v,h,M)','x','v','h','M');
    end
    try
        feval(M(i).pg,M(i).x,M(i + 1).v,M(i).gE,M(i));
    catch
        M(i).pg = inline('spm_LAP_pg(x,v,h,M)','x','v','h','M');
    end
end


% order parameters (d = n = 1 for static models) and checks
%==========================================================================
d   = M(1).E.d + 1;                          % embedding order of q(v)
n   = M(1).E.n + 1;                          % embedding order of q(x)
s   = M(1).E.s;                              % smoothness - s.d. (bins)

% number of states and parameters - generative model
%--------------------------------------------------------------------------
ns  = size(C,2);                             % number of samples
nl  = size(M,2);                             % number of levels
nv  = sum(spm_vec(M.m));                     % number of v (casual states)
nx  = sum(spm_vec(M.n));                     % number of x (hidden states)
ny  = M(1).l;                                % number of y (inputs)
nc  = M(end).l;                              % number of c (prior causes)
nu  = nv*d + nx*n;                           % number of generalised states
ne  = nv*n + nx*n + ny*n;                    % number of generalised errors


% number of states and parameters - generative process
%--------------------------------------------------------------------------
gv  = sum(spm_vec(G.l));                     % number of v (outputs)
ga  = sum(spm_vec(G.k));                     % number of a (active states)
gx  = sum(spm_vec(G.n));                     % number of x (hidden states)
gy  = ny;                                    % number of y (inputs)
na  = ga;                                    % number of a (action)

% precision (R) of generalised errors and null matrices for concatenation
%==========================================================================
Rh  = spm_DEM_R(n,s,form);
Rg  = spm_DEM_R(n,s,form);

QW  = sparse(nx*n,nx*n);
QV  = sparse((ny + nv)*n,(ny + nv)*n);


% restriction matrix, mapping prediction errors to action
%--------------------------------------------------------------------------
for i = 1:nl
    Qh{i,i} = sparse(M(i).l,M(i).l);
    Qg{i,i} = sparse(M(i).n,M(i).n);
end
Qh{1} = G(1).U;
iG    = blkdiag(kron(Rh,spm_cat(Qh)),kron(Rg,spm_cat(Qg)));


% fixed priors on states (u)
%--------------------------------------------------------------------------
Px    = kron(spm_DEM_R(n,2),spm_cat(spm_diag({M(1:end).xP})));
Pv    = kron(spm_DEM_R(d,2),spm_cat(spm_diag({M(2:end).vP})));
Pu    = spm_cat(spm_diag({Px Pv}));
Pa    = spm_speye(na,na)*exp(-2);


% hyperpriors
%--------------------------------------------------------------------------
ph.h  = spm_vec({M.hE M.gE});                % prior expectation of h,g
ph.c  = spm_cat(spm_diag({M.hC M.gC}));      % prior covariances of h,g
Ph    = spm_inv(ph.c);                       % prior precision of h,g

qh.h  = {M.hE};                              % conditional expectation h
qh.g  = {M.gE};                              % conditional expectation g
nh    = length(spm_vec(qh.h));               % number of hyperparameters h
ng    = length(spm_vec(qh.g));               % number of hyperparameters g
nb    = nh + ng;                             % number of hyperparameters


% priors on parameters (in reduced parameter space)
%==========================================================================
pp.c  = cell(nl,nl);
qp.p  = cell(nl,1);
for i = 1:(nl - 1)
    
    % eigenvector reduction: p <- pE + qp.u*qp.p
    %----------------------------------------------------------------------
    qp.u{i}   = spm_svd(M(i).pC,0);          % basis for parameters
    M(i).p    = size(qp.u{i},2);             % number of qp.p
    qp.p{i}   = sparse(M(i).p,1);            % initial deviates
    pp.c{i,i} = qp.u{i}'*M(i).pC*qp.u{i};    % prior covariance
    
end
Up    = spm_cat(spm_diag(qp.u));

% priors on parameters
%--------------------------------------------------------------------------
pp.p  = spm_vec(M.pE);
pp.c  = spm_cat(pp.c);
Pp    = spm_inv(pp.c);


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
qu.a      = cell(1,1);
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
[pu.v{:}] = deal(sparse(gv,1));
[pu.x{:}] = deal(sparse(gx,1));
[pu.z{:}] = deal(sparse(gv,1));
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
Dw    = kron(spm_speye(n,n,1),spm_speye(gx,gx,0));
Dz    = kron(spm_speye(n,n,1),spm_speye(gv,gv,0));
Iw    = kron(spm_speye(n,n,0),spm_speye(gx,gx,0));
Du    = spm_cat(spm_diag({Dx,Dv}));
Ib    = spm_speye(np + nb,np + nb);

dbdt  = sparse(np + nb,1);
dydv  = kron(speye(n,n),speye(gy,gv));


% gradients of generalised weighted errors
%--------------------------------------------------------------------------
dedh  = sparse(nh,ne);
dedg  = sparse(ng,ne);
dedv  = sparse(nv,ne);
dedx  = sparse(nx,ne);
dedhh = sparse(nh,nh);
dedgg = sparse(ng,ng);

% curvatures of Gibb's energy w.r.t. hyperparameters
%--------------------------------------------------------------------------
dHdh  = sparse(nh,1);
dHdg  = sparse(ng,1);
dHdp  = sparse(np,1);
dHdu  = sparse(nu,1);


% test for dependency of precisions on hyperparameters and states
%--------------------------------------------------------------------------
[p,dp]        = spm_LAP_eval(M,qu,qh);
try method.h  = M(1).E.method.h; catch, method.h = any(dp.h.dh(:)); end
try method.g  = M(1).E.method.g; catch, method.g = any(dp.g.dg(:)); end
try method.x  = M(1).E.method.x; catch, method.x = any([dp.g.dx(:);dp.h.dx(:)]); end
try method.v  = M(1).E.method.v; catch, method.v = any([dp.g.dv(:);dp.h.dv(:)]); end
M(1).E.method = method;


% preclude unnecessary iterations and set switches
%--------------------------------------------------------------------------
mnh   = nh*method.h;
mng   = ng*method.g;
mnx   = nx*method.x;
mnv   = nv*method.v;
if ~np && ~logical(mnh) && ~logical(mng), nN = 1; end


% preclude very precise states from entering free-energy/action
%--------------------------------------------------------------------------
ih    = p.h < 8;
ig    = p.g < 8;
ie    = kron(ones(n,1),ih);
ix    = kron(ones(n,1),ig);
iv    = kron(ones(d,1),ih((1:nv) + ny));
je    = find([ie; ix]); ix(1:nx) = 1;
ju    = find([ix; iv]);

% and other useful indices
%--------------------------------------------------------------------------
ix    = (1:nx);
ih    = (1:nb);
iv    = (1:nv) + nx*n;
ip    = (1:np) + nu;

% create innovations (and add causes)
%--------------------------------------------------------------------------
[z w]  = spm_DEM_z(G,ns);
z{end} = C + z{end};
a      = {G.a};
Z      = spm_cat(z(:));
W      = spm_cat(w(:));
A      = spm_cat(a(:));


% number of iterations for convergence
%--------------------------------------------------------------------------
convergence = -4;

% Iterate Laplace scheme
%==========================================================================
F      = -Inf;
for iN = 1:nN
    
    % get time and clear persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic; clear spm_DEM_eval qa
    
    % [re-]set states & their derivatives
    %----------------------------------------------------------------------
    try
        pu = Q(1).r;
        qu = Q(1).u; 
    end
    
    
    % D-Step: (over time)
    %======================================================================
    for is = 1:ns
        
        
        % pass action to pu.a (external states)
        %==================================================================
        if exist('qa','var'), A = spm_cat({qa,qu.a}); end
        
        % derivatives of responses and random fluctuations
        %------------------------------------------------------------------
        pu.z = spm_DEM_embed(Z,n,is);
        pu.w = spm_DEM_embed(W,n,is);
        pu.a = spm_DEM_embed(A,n,is);
        qu.u = spm_DEM_embed(U,n,is);
                
        % evaluate generative process
        %------------------------------------------------------------------
        [pu dg df] = spm_ADEM_diff(G,pu);
        
        
        % and pass response to qu.y
        %==================================================================
        for i = 1:n
            y       = spm_unvec(pu.v{i},{G.v});
            qu.y{i} = y{1};
        end
        
        % evaluate recognition model functions and derivatives
        %==================================================================
        
        % prediction errors (E) and precision vectors (p)
        %------------------------------------------------------------------
        [E dE] = spm_DEM_eval(M,qu,qp);
        [p dp] = spm_LAP_eval(M,qu,qh);
        
        
        % gradients of log(det(iS)) dDd...
        %==================================================================
        
        % get precision matrices
        %------------------------------------------------------------------
        iSh   = diag(exp(p.h));
        iSg   = diag(exp(p.g));
        iS    = blkdiag(kron(Rh,iSh),kron(Rg,iSg));
        
        
        % gradients of trace(diag(p)) = sum(p); p = precision vector
        %------------------------------------------------------------------
        dpdx  = n*sum(spm_cat({dp.h.dx; dp.g.dx}));
        dpdv  = n*sum(spm_cat({dp.h.dv; dp.g.dv}));
        dpdh  = n*sum(dp.h.dh);
        dpdg  = n*sum(dp.g.dg);
        dpdx  = kron(sparse(1,1,1,1,n),dpdx);
        dpdv  = kron(sparse(1,1,1,1,d),dpdv);
        dDdu  = [dpdx dpdv]';
        dDdh  = [dpdh dpdg]';
        
        
        % gradients precision-weighted generalised error dSd..
        %==================================================================
        
        % gradients w.r.t. hyperparameters
        %------------------------------------------------------------------
        for i = 1:nh
            diS       = diag(dp.h.dh(:,i).*exp(p.h));
            diSdh{i}  = blkdiag(kron(Rh,diS),QW);
            dedh(i,:) = E'*diSdh{i};
        end
        for i = 1:ng
            diS       = diag(dp.g.dg(:,i).*exp(p.g));
            diSdg{i}  = blkdiag(QV,kron(Rg,diS));
            dedg(i,:) = E'*diSdg{i};
        end
        
        % gradients w.r.t. hidden states
        %------------------------------------------------------------------
        for i = 1:mnx
            diV       = diag(dp.h.dx(:,i).*exp(p.h));
            diW       = diag(dp.g.dx(:,i).*exp(p.g));
            diSdx{i}  = blkdiag(kron(Rh,diV),kron(Rg,diW));
            dedx(i,:) = E'*diSdx{i};
        end
        
        % gradients w.r.t. causal states
        %------------------------------------------------------------------
        for i = 1:mnv
            diV       = diag(dp.h.dv(:,i).*exp(p.h));
            diW       = diag(dp.g.dv(:,i).*exp(p.g));
            diSdv{i}  = blkdiag(kron(Rh,diV),kron(Rg,diW));
            dedv(i,:) = E'*diSdv{i};
        end
        
        dSdx  = kron(sparse(1,1,1,n,1),dedx);
        dSdv  = kron(sparse(1,1,1,d,1),dedv);
        dSdu  = [dSdx; dSdv];
        dEdh  = [dedh; dedg];
        dEdp  = dE.dp'*iS;
        dEdu  = dE.du'*iS;
        
        % curvatures w.r.t. hyperparameters
        %------------------------------------------------------------------
        for i = 1:nh
            for j = i:nh
                diS        = diag(dp.h.dh(:,i).*dp.h.dh(:,j).*exp(p.h));
                diS        = blkdiag(kron(Rh,diS),QW);
                dedhh(i,j) = E'*diS*E;
                dedhh(j,i) = dedhh(i,j);
            end
        end
        for i = 1:ng
            for j = i:ng
                diS        = diag(dp.g.dg(:,i).*dp.g.dg(:,j).*exp(p.g));
                diS        = blkdiag(QV,kron(Rg,diS));
                dedgg(i,j) = E'*diS*E;
                dedgg(j,i) = dedgg(i,j);
            end
        end
        
        % combined curvature
        %------------------------------------------------------------------
        dSdhh = spm_cat({dedhh  [] ;
                         [] dedgg});
        
        
        % errors (from prior expectations) (NB pp.p = 0)
        %------------------------------------------------------------------
        Eu    = spm_vec(qu.x(1:n),qu.v(1:d));
        Ep    = spm_vec(qp.p);
        Eh    = spm_vec(qh.h,qh.g) - ph.h;
        
        
        % first-order derivatives of Gibb's Energy
        %==================================================================
        dLdu  = dEdu*E + dSdu*E/2 - dDdu/2 + Pu*Eu;
        dLdh  = dEdh*E/2          - dDdh/2 + Ph*Eh;
        dLdp  = dEdp*E                     + Pp*Ep;
        
        
        % and second-order derivatives of Gibb's Energy
        %------------------------------------------------------------------
        dLduu = dEdu*dE.du + Pu;
        dLdpp = dEdp*dE.dp + Pp;
        dLdhh = dSdhh/2    + Ph;
        dLdup = dEdu*dE.dp;
        dLdhp = dEdh*dE.dp;
        dLdpu = dLdup';
        dLdph = dLdhp';
 
 
        % precision and covariances for entropy
        %------------------------------------------------------------------                      
        dLdaa = spm_cat({dLduu dLdup  ;
                         dLdpu dLdpp});
        dLdbb = spm_cat({dLdpp  dLdph ;
                         dLdhp  dLdhh});
            
        Cup   = spm_inv(dLdaa);
        Chh   = spm_inv(dLdhh);
            
            
        % first-order derivatives of Entropy term
        %==================================================================
        
        % log-precision
        %------------------------------------------------------------------
        for i = 1:nh
            Luub    = dE.du'*diSdh{i}*dE.du;
            Lpub    = dE.dp'*diSdh{i}*dE.du;
            Lppb    = dE.dp'*diSdh{i}*dE.dp;
            diCdh   = spm_cat({Luub Lpub';
                               Lpub Lppb});
            dHdh(i) = spm_trace(diCdh,Cup)/2;
        end
        for i = 1:ng
            Luub    = dE.du'*diSdg{i}*dE.du;
            Lpub    = dE.dp'*diSdg{i}*dE.du;
            Lppb    = dE.dp'*diSdg{i}*dE.dp;
            diCdg   = spm_cat({Luub Lpub';
                               Lpub Lppb});
            dHdg(i) = spm_trace(diCdg,Cup)/2;
        end
        
        % parameters
        %------------------------------------------------------------------
        for i = 1:np
            Luup    = dE.dup{i}'*dEdu';
            Lpup    = dEdp*dE.dup{i};
            Luup    = Luup + Luup';
            diCdp   = spm_cat({Luup Lpup';
                               Lpup [] });
            dHdp(i) = spm_trace(diCdp,Cup)/2;
        end
        
        % and concatenate
        %------------------------------------------------------------------
        dHdb  = [dHdh; dHdg];
        dHdb  = [dHdp; dHdb];
        dLdb  = [dLdp; dLdh];
        
        
        % whiten generalised ascent on parameters and hyperparameters
        %==================================================================
        
        % accumulate curvatures of [hyper] parameters
        %------------------------------------------------------------------
        try
            dLdB  = dLdB*(1 - 1/ns)  + dLdb/ns;
            dLdBB = dLdBB*(1 - 1/ns) + dLdbb/ns;
        catch
            dLdB  = dLdb - dLdb;
            dLdBB = Ib*32;
        end
        
        % whiten gradient (and curvatures) with regularised precision
        %------------------------------------------------------------------
        Cb    = spm_inv(dLdBB + Ib*exp(dt));
        dLdb  = Cb*dLdB;
        dHdb  = Cb*dHdb;
        
        % prior precision of fluctuations on [hyper] parameters
        %------------------------------------------------------------------
        Kb    = ns*Ib;  
        
        % derivatives and curvature generative process (and action)
        %==================================================================
        
        % tensor products for Jacobian
        %------------------------------------------------------------------
        Dgda  = kron(spm_speye(n,1,1),dg.da);
        Dgdv  = kron(spm_speye(n,n,1),dg.dv);
        Dgdx  = kron(spm_speye(n,n,1),dg.dx);
        dfda  = kron(spm_speye(n,1,0),df.da);
        dfdv  = kron(spm_speye(n,n,0),df.dv);
        dfdx  = kron(spm_speye(n,n,0),df.dx);
        
        dgda  = kron(spm_speye(n,1,0),dg.da);
        dgdx  = kron(spm_speye(n,n,0),dg.dx);
        
        % change in error w.r.t. action
        %------------------------------------------------------------------
        Dfdx  = 0;
        for i = 1:n
            Dfdx = Dfdx + kron(spm_speye(n,n,-i),df.dx^(i - 1));
        end
        
        % dE/da with restriction
        %------------------------------------------------------------------
        dE.dv = dE.dy*dydv;
        dE.da = dE.dv*(dgda + dgdx*Dfdx*dfda);
        
        % first-order derivatives
        %------------------------------------------------------------------
        dVda  = -dE.da'*iG*E - Pa*spm_vec(qu.a{1:1});
        

        % and second-order derivatives
        %------------------------------------------------------------------
        dVdaa = -dE.da'*iG*dE.da - Pa;
        dVduv = -dE.du'*iS*dE.dv;
        dVduc = -dE.du'*iS*dE.dc;
        dVdua = -dE.du'*iS*dE.da;
        dVdav = -dE.da'*iG*dE.dv;
        dVdau = -dE.da'*iG*dE.du;
        dVdac = -dE.da'*iG*dE.dc;
        
  
        % save conditional moments (and prediction error) at Q{t}
        %==================================================================
        
        % save means
        %------------------------------------------------------------------
        Q(is).E   = diag(diag(iS))*E;
        Q(is).e   = E;
        Q(is).u   = qu;
        Q(is).p   = qp;
        Q(is).h   = qh;
        Q(is).r   = pu;
        
        % and action
        %------------------------------------------------------------------
        if na, qa(:,is) = qu.a{1}; end
        
        % and conditional covariances
        %------------------------------------------------------------------
        Q(is).u.s = Cup(ix,ix);
        Q(is).u.c = Cup(iv,iv);
        Q(is).p.c = Cup(ip,ip);
        Q(is).h.c = Chh(ih,ih);
        
        % Free-energy (components)
        %------------------------------------------------------------------
        Fc(is,1)  = - E(je)'*iS(je,je)*E(je)/2;
        Fc(is,2)  = - Eu(ju)'*Pu(ju,ju)*Eu(ju)/2;
        Fc(is,3)  = - n*ny*log(2*pi)/2;
        Fc(is,4)  = spm_logdet(iS(je,je))/2;
        Fc(is,5)  = spm_logdet(Pu(ju,ju)*Cup(ju,ju))/2;
        
        
        % update conditional moments
        %==================================================================
        
        % assemble true states and conditional means
        %------------------------------------------------------------------
        r.v  = pu.v(1:n);
        r.x  = pu.x(1:n);
        r.z  = pu.z(1:n);
        r.w  = pu.w(1:n);
        
        q.x  = qu.x(1:n);
        q.v  = qu.v(1:d);
        q.c  = qu.u(1:d);
        q.a  = qu.a(1:1);
        q.p  = qp.p;
        q.h  = qh.h;
        q.g  = qh.g;
        q.d  = dbdt;
        
        
        % flow
        %------------------------------------------------------------------
        g.v  =  Dz*spm_vec(r.v)                  ;
        g.x  =  Dw*spm_vec(r.x)                  ;
        g.z  =  Dz*spm_vec(r.z)                  ;
        g.w  =  Dw*spm_vec(r.w)                  ;
        
        f.u  =  Du*spm_vec(q.x,q.v) - dLdu - dHdu;
        f.c  =  Dc*spm_vec(q.c)                  ;
        f.a  =                        dVda       ;
        f.b  =     spm_vec(q.d)                  ;
        f.d  = -Kb*spm_vec(q.d)     - dLdb - dHdb;
           
        
        % Jacobian (variational flow)
        %------------------------------------------------------------------
        dfdq = {...
            Dgdv  Dgdx Dz   []   []       []    Dgda   []  [];
            dfdv  dfdx []   Iw   []       []    dfda   []  [];
            []    []   Dz   []   []       []    []     []  [];
            []    []   []   Dw   []       []    []     []  [];
            dVduv []   []   []   Du-dLduu dVduc dVdua  []  [];
            []    []   []   []   []       Dc    []     []  [];
            dVdav []   []   []   dVdau    dVdac dVdaa  []  [];
            []    []   []   []   []       []    []     []  Ib;
            []    []   []   []   []       []    []    -Ib -Kb};
        
        
        % update conditional modes of states
        %==================================================================
        dq        = spm_dx(spm_cat(dfdq),spm_vec(g,f),1);
        [r,q]     = spm_unvec(spm_vec(r,q) + dq,r,q);

        % unpack conditional means
        %------------------------------------------------------------------
        pu.v(1:n) = r.v;
        pu.x(1:n) = r.x;
        
        qu.x(1:n) = q.x;
        qu.v(1:d) = q.v;
        qu.a(1:1) = q.a;
        qp.p      = q.p;
        qh.h      = q.h;
        qh.g      = q.g;
        dbdt      = q.d;
        
        
    end % sequence (ns)
    
    
    % Bayesian parameter averaging
    %======================================================================
    
    % Conditional moments of time-averaged parameters
    %----------------------------------------------------------------------
    Ep    = 0;
    Qp    = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).p.c);
        Ep  = Ep  + P*spm_vec(Q(i).p.p);
        Qp  = Qp + P;
    end
    Ep    = spm_inv(Qp)*Ep;
    Cp    = spm_inv(Qp + (1 - ns)*Pp);
    
    % conditional moments of hyper-parameters
    %----------------------------------------------------------------------
    Eh    = 0;
    Qh    = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).h.c);
        Eh  = Eh  + P*spm_vec({Q(i).h.h Q(i).h.g});
        Qh  = Qh + P;
    end
    Eh    = spm_inv(Qh)*Eh - ph.h;
    Ch    = spm_inv(Qh + (1 - ns)*Ph);
    
    
    
    % Free-action of states plus free-energy of parameters
    %======================================================================
    FC(1) = sum(Fc(:,1));       % - E'*iS*E/2;
    FC(2) = sum(Fc(:,2));       % - Eu'*Pu*Eu/2;
    FC(3) = sum(Fc(:,3));       % - n*ny*log(2*pi)/2;
    FC(4) = sum(Fc(:,4));       %   spm_logdet(iS)/2;
    FC(5) = sum(Fc(:,5));       %   spm_logdet(Pu*Cu)/2;
    FC(6) = -Ep'*Pp*Ep/2;
    FC(7) = -Eh'*Ph*Eh/2;
    FC(8) = spm_logdet(Pp*Cp)/2;
    FC(9) = spm_logdet(Ph*Ch)/2;
    
    Fe    = sum(FC);
    
    % if F is decreasing, revert [hyper] parameters and slow down
    %----------------------------------------------------------------------
    if Fe < F(iN) && iN > 4
        
        % save free-energy
        %------------------------------------------------------------------
        F(iN + 1) = F(iN);
        
        % load current MAP estimates
        %------------------------------------------------------------------
        qp = PQ.qp;
        qh = PQ.qh;
        
        % decrease update time
        %------------------------------------------------------------------
        dt = max(dt + 2,2);
        
        % convergence
        %------------------------------------------------------------------
        if dt > 6; convergence = 1; end
        
    else
        
        % convergence
        %------------------------------------------------------------------
        if Fe - F(iN) < 1, convergence = convergence + 1; end
        
        % save free-energy
        %------------------------------------------------------------------
        F(iN)     = Fe;
        F(iN + 1) = Fe;
        
        % save current MAP estimates
        %------------------------------------------------------------------
        PQ.qp = qp;
        PQ.qh = qh;
        
        % increase update time
        %------------------------------------------------------------------
        dt = max(dt - 1,-8);
        
    end
    
    % Convergence
    %======================================================================
    if convergence > 0; break, end
    
    % otherwise save conditional moments (for each time point)
    %======================================================================
    for t = 1:length(Q)
        
        
        % states
        %------------------------------------------------------------------
        a     = spm_unvec(Q(t).u.a{1},{G.a});
        v     = spm_unvec(Q(t).r.v{1},{G.v});
        x     = spm_unvec(Q(t).r.x{1},{G.x});
        z     = spm_unvec(Q(t).r.z{1},{G.v});
        w     = spm_unvec(Q(t).r.w{1},{G.x});
        for i = 1:nl
            try
                pU.v{i}(:,t) = spm_vec(v{i});
                pU.z{i}(:,t) = spm_vec(z{i});
            end
            try
                pU.x{i}(:,t) = spm_vec(x{i});
                pU.w{i}(:,t) = spm_vec(w{i});
            end
            try
                qU.a{i}(:,t) = spm_vec(a{i});
            end
        end
 
        % states and predictions
        %------------------------------------------------------------------
        v     = spm_unvec(Q(t).u.v{1},{M(1 + 1:end).v});
        x     = spm_unvec(Q(t).u.x{1},{M(1:end - 1).x});
        z     = spm_unvec(Q(t).e(1:(ny + nv)),{M.v});
        e     = spm_unvec(Q(t).E(1:(ny + nv)),{M.v});
        w     = spm_unvec(Q(t).e((1:nx) + (ny + nv)*n),{M.x});
        u     = spm_unvec(Q(t).E((1:nx) + (ny + nv)*n),{M.x});
        for i = 1:(nl - 1)
            if M(i).m, qU.v{i + 1}(:,t) = spm_vec(v{i});  end
            if M(i).n, qU.x{i}(:,t)     = spm_vec(x{i});  end
            if M(i).n, qU.w{i}(:,t)     = spm_vec(w{i});  end
            if M(i).l, qU.z{i}(:,t)     = spm_vec(z{i});  end
            if M(i).n, qU.W{i}(:,t)     = spm_vec(u{i});  end
            if M(i).l, qU.Z{i}(:,t)     = spm_vec(e{i});  end
        end
        if    M(nl).l, qU.z{nl}(:,t)    = spm_vec(z{nl}); end
        if    M(nl).l, qU.Z{nl}(:,t)    = spm_vec(e{nl}); end
        
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
        qH.p{t} = spm_vec({Q(t).h.h Q(t).h.g});
        qH.c{t} = Q(t).h.c;
        
    end
    
    % graphics (states)
    %----------------------------------------------------------------------
    spm_figure('GetWin','GF');
    spm_DEM_qU(qU)
    
    % graphics (parameters and log-precisions)
    %----------------------------------------------------------------------
    if np && nb
        subplot(2*nl,2,4*nl - 2)
        plot(1:ns,spm_cat(qP.p))
        set(gca,'XLim',[1 ns])
        title('parameters (modes)','FontSize',16)
        
        subplot(2*nl,2,4*nl)
        plot(1:ns,spm_cat(qH.p))
        set(gca,'XLim',[1 ns])
        title('log-precision','FontSize',16)
        
    elseif nb
        subplot(nl,2,2*nl)
        plot(1:ns,spm_cat(qH.p))
        set(gca,'XLim',[1 ns])
        title('log-precision','FontSize',16)
        
    elseif np
        subplot(nl,2,2*nl)
        plot(1:ns,spm_cat(qP.p))
        set(gca,'XLim',[1 ns])
        title('parameters (modes)','FontSize',16)
        
    end
    drawnow

    
    % report (EM-Steps)
    %----------------------------------------------------------------------
    try
        dF = F(iN) - F(iN - 1);
    catch
        dF = 0;
    end
    str{1} = sprintf('LAP: %i', iN);
    if iN == 1
        str{2} = sprintf('  F0:%.4e', full(F(iN)));
    else
        str{2} = sprintf('F-F0:%.4e', full(F(iN) - F(1)));
    end
    str{3} = sprintf('dF:%.2e',      full(dF));
    str{4} = sprintf('(%.2e sec)',   full(toc));
    fprintf('%-16s%-20s%-14s%-16s\n',str{:})
    
end
fprintf('%-19sF:%.4e\n', 'LAP: Converged', F(end));

% Place Bayesian parameter averages in output arguments
%==========================================================================

% Conditional moments of time-averaged parameters
%--------------------------------------------------------------------------
Qp = 0;
Ep = 0;
for i = 1:ns
    
    % weight in proportion to precisions
    %----------------------------------------------------------------------
    P  = spm_inv(qP.c{i});
    Ep = Ep + P*qP.p{i};
    Qp = Qp + P;
    
end
Ep     = spm_inv(Qp)*Ep;
Cp     = spm_inv(Qp + (1 - ns)*Pp);
qP.P   = spm_unvec(Up*Ep + pp.p,{M.pE});
qP.C   = Up*Cp*Up';
qP.V   = spm_unvec(diag(qP.C),{M.pE});
qP.U   = Up;

% conditional moments of hyper-parameters
%--------------------------------------------------------------------------
Qh = 0;
Eh = 0;
for i = 1:ns
    
    % weight in proportion to precisions
    %----------------------------------------------------------------------
    P  = spm_inv(qH.c{i});
    Eh = Eh + P*qH.p{i};
    Qh = Qh + P;
    
end
Eh     = spm_inv(Qh)*Eh;
Ch     = spm_inv(Qh + (1 - ns)*Ph);
P      = spm_unvec(Eh,{qh.h qh.g});
qH.h   = P{1};
qH.g   = P{2};
qH.C   = Ch;
P      = spm_unvec(diag(qH.C),P);
qH.V   = P{1};
qH.W   = P{2};

% Fill in DEM with response and its causes
%--------------------------------------------------------------------------
DEM.M    = M;                   % model
DEM.U    = U;                   % causes
DEM.Y    = pU.v{1};             % simulated response variables

DEM.qU   = qU;                  % conditional moments of model-states
DEM.qP   = qP;                  % conditional moments of model-parameters
DEM.qH   = qH;                  % conditional moments of hyper-parameters
DEM.pU   = pU;                  % true states
DEM.pP.P = {G.pE};              % true parameters

DEM.F    = F(1:iN);             % [-ve] Free-energy

