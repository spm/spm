function [DEM] = spm_LAPF(DEM)
% Laplacian model inversion (see also spm_LAPS)
% FORMAT DEM   = spm_LAPF(DEM)
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
% spm_LAPF implements a variational scheme under the Laplace
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
% $Id: spm_LAPF.m 6018 2014-05-25 09:24:14Z karl $


% find or create a DEM figure
%--------------------------------------------------------------------------
try
    DEM.M(1).nograph;
catch
    DEM.M(1).nograph = 0;
end
if ~DEM.M(1).nograph
    Fdem = spm_figure('GetWin','DEM');
end


% check model, data and priors
%==========================================================================
[M,Y,U] = spm_DEM_set(DEM);


% number of iterations
%--------------------------------------------------------------------------
try, nD = M(1).E.nD; catch, nD = 1;   end
try, nN = M(1).E.nN; catch, nN = 16;  end


% ensure integration scheme evaluates gradients at each time-step
%--------------------------------------------------------------------------
M(1).E.linear = 4;

% assume precisions are a function of, and only of, hyperparameters
%--------------------------------------------------------------------------
try
    method = M(1).E.method;
catch
    method.h = 1;
    method.g = 1;
    method.x = 0;
    method.v = 0;
end
try method.h; catch, method.h = 0; end
try method.g; catch, method.g = 0; end
try method.x; catch, method.x = 0; end
try method.v; catch, method.v = 0; end


% assume precisions are a function of, and only of, hyperparameters
%--------------------------------------------------------------------------
try
    form = M(1).E.form;
catch
    form = 'Gaussian';
end

% checks for Laplace models (precision functions; ph and pg)
%--------------------------------------------------------------------------
for i  = 1:length(M)
    try
        feval(M(i).ph,M(i).x,M(i + 1).v,M(i).hE,M(i)); method.v = 1;
    catch
        M(i).ph = inline('spm_LAP_ph(x,v,h,M)','x','v','h','M');
    end
    try
        feval(M(i).pg,M(i).x,M(i + 1).v,M(i).gE,M(i)); method.x = 1;
    catch
        M(i).pg = inline('spm_LAP_pg(x,v,h,M)','x','v','h','M');
    end
end

M(1).E.method = method;

 
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
s     = M(1).E.s;
Rh    = spm_DEM_R(n,s,form);
Rg    = spm_DEM_R(n,s,form);

W     = sparse(nx*n,nx*n);
V     = sparse((ny + nv)*n,(ny + nv)*n);


% fixed priors on states (u)
%--------------------------------------------------------------------------
Px    = kron(sparse(1,1,1,n,n),spm_cat(spm_diag({M.xP})));
Pv    = kron(sparse(1,1,1,d,d),sparse(nv,nv));
pu.ic = spm_cat(spm_diag({Px Pv}));
 
% hyperpriors
%--------------------------------------------------------------------------
ph.h  = spm_vec({M.hE M.gE});                % prior expectation of h,g
ph.c  = spm_cat(spm_diag({M.hC M.gC}));      % prior covariances of h,g
ph.ic = spm_pinv(ph.c);                      % prior precision of h,g
 
qh.h  = {M.hE};                              % conditional expectation h
qh.g  = {M.gE};                              % conditional expectation g
nh    = length(spm_vec(qh.h));               % number of hyperparameters h
ng    = length(spm_vec(qh.g));               % number of hyperparameters g
nb    = nh + ng;                             % number of hyerparameters


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
Ip     = spm_speye(np,np);
Ih     = spm_speye(nb,nb);
qp.dp  = sparse(np,1);                   % conditional expectation of dp/dt
qh.dp  = sparse(nb,1);                   % conditional expectation of dh/dt

% precision of fluctuations on parameters of hyperparameters
%--------------------------------------------------------------------------  
Kp     = ns*Ip;
Kh     = ns*Ih;

% gradients of generalised weighted errors
%--------------------------------------------------------------------------
dedh   = sparse(nh,ne);
dedg   = sparse(ng,ne);
dedv   = sparse(nv,ne);
dedx   = sparse(nx,ne);
dedhh  = sparse(nh,nh);
dedgg  = sparse(ng,ng);
            
% curvatures of Gibb's energy w.r.t. hyperparameters
%--------------------------------------------------------------------------
dHdh   = sparse(nh,  1);
dHdg   = sparse(ng,  1);
dHdp   = sparse(np,  1);
dHdx   = sparse(nx*n,1);
dHdv   = sparse(nv*d,1);

% preclude unnecessary iterations and set switchs
%--------------------------------------------------------------------------
if ~np && ~nh && ~ng, nN = 1; end
mnx = nx*~~method.x;
mnv = nv*~~method.v;


% Iterate Lapalace scheme
%==========================================================================
Fa     = -Inf;
for iN = 1:nN
 
    % get time and clear persistent variables in evaluation routines
    %----------------------------------------------------------------------
    tic; clear spm_DEM_eval
 
    % [re-]set states & their derivatives
    %----------------------------------------------------------------------
    try, qu = Q(1).u; end
    
    
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
            
            % prediction errors (E) and precision vectors (p)
            %--------------------------------------------------------------
            [E,dE]  = spm_DEM_eval(M,qu,qp);
            [p,dp]  = spm_LAP_eval(M,qu,qh);
            
 
            % gradients of log(det(iS)) dDd...
            %==============================================================
            
            % get precision matrices
            %--------------------------------------------------------------
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
            dDdh    = [dpdh dpdg]';
 
            
            % gradients precision-weighted generalised error dSd..
            %==============================================================

            % gradients w.r.t. hyperparameters
            %--------------------------------------------------------------
            for i = 1:nh
                diS       = diag(dp.h.dh(:,i).*exp(p.h));
                diSdh{i}  = blkdiag(kron(Rh,diS),W);
                dedh(i,:) = E'*diSdh{i};
            end
            for i = 1:ng
                diS       = diag(dp.g.dg(:,i).*exp(p.g));
                diSdg{i}  = blkdiag(V,kron(Rg,diS));
                dedg(i,:) = E'*diSdg{i};
            end

            % gradients w.r.t. hidden states
            %--------------------------------------------------------------
            for i = 1:mnx
                diV       = diag(dp.h.dx(:,i).*exp(p.h));
                diW       = diag(dp.g.dx(:,i).*exp(p.g));
                diSdx{i}  = blkdiag(kron(Rh,diV),kron(Rg,diW));
                dedx(i,:) = E'*diSdx{i};
            end
            
            % gradients w.r.t. causal states
            %--------------------------------------------------------------
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
            %--------------------------------------------------------------
            for i = 1:nh
                for j = i:nh
                    diS        = diag(dp.h.dh(:,i).*dp.h.dh(:,j).*exp(p.h));
                    diS        = blkdiag(kron(Rh,diS),W);
                    dedhh(i,j) = E'*diS*E;
                    dedhh(j,i) = dedhh(i,j);
                end
            end            
            for i = 1:ng
                for j = i:ng
                    diS        = diag(dp.g.dg(:,i).*dp.g.dg(:,j).*exp(p.g));
                    diS        = blkdiag(V,kron(Rg,diS));
                    dedgg(i,j) = E'*diS*E;
                    dedgg(j,i) = dedgg(i,j);
                end
            end
            
            % combined curvature
            %--------------------------------------------------------------
            dSdhh = spm_cat({dedhh  []     ;
                             []     dedgg});
                 
            
            % errors (from prior expectations) (NB pp.p = 0)
            %--------------------------------------------------------------
            Eu    = spm_vec(qu.x(1:n),qu.v(1:d));
            Ep    = spm_vec(qp.p);
            Eh    = spm_vec(qh.h,qh.g) - ph.h;
            

            % first-order derivatives of Gibb's Energy
            %==============================================================
            dLdu  = dEdu*E + dSdu*E/2 - dDdu/2 + pu.ic*Eu;
            dLdh  = dEdh*E/2          - dDdh/2 + ph.ic*Eh;
            dLdp  = dEdp*E                     + pp.ic*Ep;
            
 
            % and second-order derivatives of Gibb's Energy
            %--------------------------------------------------------------
            % dLduu = dEdu*dE.du + dSdu*dE.du + dE.du'*dSdu' + pu.ic;
            % dLdup = dEdu*dE.dp + dSdu*dE.dp;
            dLduu = dEdu*dE.du + pu.ic;
            dLdpp = dEdp*dE.dp + pp.ic;
            dLdhh = dSdhh/2    + ph.ic;
            dLdup = dEdu*dE.dp;
            dLdhu = dEdh*dE.du;
            dLduy = dEdu*dE.dy;
            dLduc = dEdu*dE.dc;
            dLdpy = dEdp*dE.dy;
            dLdpc = dEdp*dE.dc;
            dLdhy = dEdh*dE.dy;
            dLdhc = dEdh*dE.dc;
            dLdhp = dEdh*dE.dp;
            dLdpu = dLdup';
            dLdph = dLdhp';
            
            % precision and covariances
            %--------------------------------------------------------------                        
            iC    = spm_cat({dLduu dLdup;
                             dLdpu dLdpp});
            
            C     = spm_inv(iC);
            
            % first-order derivatives of Entropy term
            %==============================================================
            
            % log-precision
            %--------------------------------------------------------------
            for i = 1:nh
                Luub    = dE.du'*diSdh{i}*dE.du;
                Lpub    = dE.dp'*diSdh{i}*dE.du;
                Lppb    = dE.dp'*diSdh{i}*dE.dp;
                diCdh   = spm_cat({Luub Lpub';
                                   Lpub Lppb});
                dHdh(i) = sum(sum(diCdh.*C))/2;
            end
            for i = 1:ng
                Luub    = dE.du'*diSdg{i}*dE.du;
                Lpub    = dE.dp'*diSdg{i}*dE.du;
                Lppb    = dE.dp'*diSdg{i}*dE.dp;
                diCdg   = spm_cat({Luub Lpub';
                                   Lpub Lppb});
                dHdg(i) = sum(sum(diCdg.*C))/2;
            end
            
            % parameters
            %--------------------------------------------------------------
            for i = 1:np
                Luup    = dE.dup{i}'*dEdu';
                Lpup    = dEdp*dE.dup{i};
                Luup    = Luup + Luup';
                diCdp   = spm_cat({Luup Lpup';
                                   Lpup [] });
                dHdp(i) = sum(sum(diCdp.*C))/2;
            end

%             % hidden and causal states
%             %--------------------------------------------------------------
%             for i = 1:mnx
%                 Luux    = dE.du'*diSdx{i}*dE.du;
%                 Lpux    = dE.dp'*diSdx{i}*dE.du;
%                 Lppx    = dE.dp'*diSdx{i}*dE.dp;
%                 diCdx   = spm_cat({Luux Lpux';
%                                    Lpux Lppx});
%                 dHdx(i) = sum(sum(diCdx.*C))/2;
%                                 
%             end
%             for i = 1:mnv
%                 Luuv    = dE.du'*diSdv{i}*dE.du;
%                 Lpuv    = dE.dp'*diSdv{i}*dE.du;
%                 Lppv    = dE.dp'*diSdv{i}*dE.dp;
%                 diCdv   = spm_cat({Luuv Lpuv';
%                                    Lpuv Lppv});
%                 dHdv(i) = sum(sum(diCdv.*C))/2;
%             end

            dHdb  = [dHdh; dHdg];
            dHdu  = [dHdx; dHdv];


            % save conditional moments (and prediction error) at Q{t}
            %==============================================================
            if iD == 1
                
                % save means
                %----------------------------------------------------------
                Q(is).e = E;
                Q(is).E = iS*E;
                Q(is).u = qu;
                Q(is).p = qp;
                Q(is).h = qh;
                                
                % and conditional covariances
                %----------------------------------------------------------
                Q(is).u.s = C((1:nx),(1:nx));
                Q(is).u.c = C((1:nv) + nx*n, (1:nv) + nx*n);
                Q(is).p.c = C((1:np) + nu,   (1:np) + nu);
                Q(is).h.c = spm_inv(dLdhh);
                Cu        = C(1:nu,1:nu);

                % Free-energy (states)
                %----------------------------------------------------------                
                L(is) = ... 
                - E'*iS*E/2      + spm_logdet(iS)/2    - n*ny*log(2*pi)/2 ...          
                - Eu'*pu.ic*Eu/2 + spm_logdet(pu.ic)/2 + spm_logdet(Cu)/2;
                    
                % Free-energy (states and parameters)
                %----------------------------------------------------------
                A(is) = - E'*iS*E/2        + spm_logdet(iS)/2    ...
                        - Eu'*pu.ic*Eu/2   + spm_logdet(pu.ic)/2 ...
                        - Ep'*pp.ic*Ep/2   + spm_logdet(pp.ic)/2 ...
                        - Eh'*ph.ic*Eh/2   + spm_logdet(ph.ic)/2 ...
                        - n*ny*log(2*pi)/2 - spm_logdet(iC)/2 - spm_logdet(dLdhh)/2;
            end
 
            % update conditional moments
            %==============================================================
            
            % uopdate curvatures of [hyper]paramters
            %--------------------------------------------------------------
            try
                dLdPP = dLdPP*(1 - 1/ns) + dLdpp/ns;
                dLdHH = dLdHH*(1 - 1/ns) + dLdhh/ns;
            catch
                dLdPP = dLdpp;
                dLdHH = dLdhh;
            end

            % rotate and scale gradient (and curvatures)
            %--------------------------------------------------------------
            [Vp,Sp] = spm_svd(dLdPP,0);
            [Vh,Sh] = spm_svd(dLdHH,0);
            Sp      = diag(1./(diag(sqrt(Sp))));
            Sh      = diag(1./(diag(sqrt(Sh))));
            
            dLdp  = Sp*Vp'*dLdp;
            dHdp  = Sp*Vp'*dHdp;
            dLdpy = Sp*Vp'*dLdpy;
            dLdpu = Sp*Vp'*dLdpu;
            dLdpc = Sp*Vp'*dLdpc;
            dLdph = Sp*Vp'*dLdph;
            dLdpp = Sp*Vp'*dLdpp*Vp;
            dLdhp =        dLdhp*Vp;
            
            dLdh  = Sh*Vh'*dLdh;
            dHdb  = Sh*Vh'*dHdb;
            dLdhy = Sh*Vh'*dLdhy;
            dLdhu = Sh*Vh'*dLdhu;
            dLdhc = Sh*Vh'*dLdhc;
            dLdhp = Sh*Vh'*dLdhp;
            dLdhh = Sh*Vh'*dLdhh*Vh;
            dLdph =        dLdph*Vh;
            
            % assemble conditional means
            %--------------------------------------------------------------
            q{1}  = qu.y(1:n);
            q{2}  = qu.x(1:n);
            q{3}  = qu.v(1:d);
            q{4}  = qu.u(1:d);
            q{5}  = spm_unvec(Vp'*spm_vec(qp.p),qp.p);
            qb    = spm_unvec(Vh'*spm_vec({qh.h qh.g}),{qh.h qh.g});
            q{6}  = qb{1};
            q{7}  = qb{2};
            q{8}  = Vp'*qp.dp;
            q{9}  = Vh'*qh.dp;
            


            % flow
            %--------------------------------------------------------------
            f{1}  =  Dy*spm_vec(q{1});
            f{2}  =  Du*spm_vec(q{2:3}) - dLdu - dHdu;
            f{3}  =  Dc*spm_vec(q{4});
            f{4}  =     spm_vec(q{8});
            f{5}  =     spm_vec(q{9});
            f{6}  = -Kp*spm_vec(q{8})   - dLdp - dHdp;
            f{7}  = -Kh*spm_vec(q{9})   - dLdh - dHdb;
            
 
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
            qp.p      = spm_unvec(Vp*spm_vec(q{5}),qp.p);
            qb        = spm_unvec(Vh*spm_vec(q{6:7}),{qh.h qh.g});
            qh.h      = qb{1};
            qh.g      = qb{2};
            qp.dp     = Vp*q{8};
            qh.dp     = Vh*q{9};

 
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
        Eh  = Eh + P*spm_vec({Q(i).h.h Q(i).h.g});
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
    if Fi < Fa && iN > 4
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
        Z     = spm_unvec(Q(t).E(1:(ny + nv)),{M.v});
        w     = spm_unvec(Q(t).e((1:nx) + (ny + nv)*n),{M.x});
        X     = spm_unvec(Q(t).E((1:nx) + (ny + nv)*n),{M.x});
        for i = 1:(nl - 1)
            if M(i).m, qU.v{i + 1}(:,t) = spm_vec(v{i});  end
            if M(i).n, qU.x{i}(:,t)     = spm_vec(x{i});  end
            if M(i).n, qU.w{i}(:,t)     = spm_vec(w{i});  end
            if M(i).l, qU.z{i}(:,t)     = spm_vec(z{i});  end
            if M(i).n, qU.W{i}(:,t)     = spm_vec(X{i});  end
            if M(i).l, qU.Z{i}(:,t)     = spm_vec(Z{i});  end
        end
        if    M(nl).l, qU.z{nl}(:,t)    = spm_vec(z{nl}); end
        if    M(nl).l, qU.Z{nl}(:,t)    = spm_vec(Z{nl}); end

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
    figure(Fdem)
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
        dF = F(end) - F(end - 1);
    catch
        dF = 0;
    end
    str{1} = sprintf('LAP: %i (%i)', iN,iD);
    str{2} = sprintf('F:%.4e',       full(F(iN) - F(1)));
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
P      = {qh.h qh.g};
P      = spm_unvec(Eh,P);
qH.h   = P{1};
qH.g   = P{2};
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



% Notes (check on curvature)
%==========================================================================


                % analytic form
                %----------------------------------------------------------
                iC = spm_cat({dLduu dLdup dLduh;
                              dLdpu dLdpp dLdph;
                              dLdhu dLdhp dLdhh});
                          

                % numerical approximations
                %----------------------------------------------------------
                qq.x  = qu.x(1:n);
                qq.v  = qu.v(1:d);
                qq.p  = qp.p;
                qq.h  = qh.h;
                qq.g  = qh.g;

                dLdqq = spm_diff('spm_LAP_F',qq,qu,qp,qh,pu,pp,ph,M,[1 1]);
                dLdqq = spm_cat(dLdqq');
                
                subplot(2,2,1);imagesc(dLdqq);     axis square
                subplot(2,2,2);imagesc(iC);        axis square
                subplot(2,2,3);imagesc(dLdqq - iC);axis square
                subplot(2,2,4);plot(iC,':k');hold on;
                plot(dLdqq - iC,'r');hold off; axis square
                drawnow
                
% Notes (descent on parameters
%==========================================================================
I     = eye(length(dLdpp));
k     = kp;
Luu   = dLdpp;
J     = spm_cat({[]    I;
               -Luu -k*I});
[u s] = eig(full(J));
max(diag(s))

[uj sj] = eig(full(dLdpp));
Luu     = min(diag(sj));
% Luu     = max(diag(sj));
k       = kp;

ss(1) = -(k + sqrt(k^2 - 4*Luu))/2;
ss(2) = -(k - sqrt(k^2 - 4*Luu))/2;
max(ss)


k    = (1:128);
s    = -(k - sqrt(k.^2 - 4*Luu))/2;

plot(k,-1./real(s))









