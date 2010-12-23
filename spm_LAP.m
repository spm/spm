function [DEM] = spm_LAP(DEM)
% Laplacian model inversion (see also spm_LAPS)
% FORMAT DEM   = spm_LAP(DEM)
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
% spm_LAP implements a variational scheme under the Laplace
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
% $Id: spm_LAP.m 4146 2010-12-23 21:01:39Z karl $
 
 
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
[M Y U] = spm_DEM_set(DEM);


% set regularisation
%--------------------------------------------------------------------------
try
    dt = DEM.M(1).E.v;
catch
    dt = 0;
    DEM.M(1).E.v = dt;
end

 
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
 
M(1).E.method = method;
 
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
        feval(M(i).ph,M(i).x,M(i).v,M(i).hE,M(i));
    catch
        M(i).ph = inline('spm_LAP_ph(x,v,h,M)','x','v','h','M');
    end
    try
        feval(M(i).pg,M(i).x,M(i).v,M(i).gE,M(i));
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
ph.ic = spm_inv(ph.c);                       % prior precision of h,g
 
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
    qp.u{i}   = spm_svd(M(i).pC,exp(-32));           % basis for parameters
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
Ib     = spm_speye(np + nb,np + nb);
dbdt   = sparse(np + nb,1);
 
 
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
dHdh   = sparse(nh,1);
dHdg   = sparse(ng,1);
dHdp   = sparse(np,1);
dHdu   = sparse(nu,1);

 
% preclude unnecessary iterations and set switches
%--------------------------------------------------------------------------
if ~np && ~nh && ~ng, nN = 1; end
mnx    = nx*~~method.x;
mnv    = nv*~~method.v;
 
% preclude very precise states from entering free-energy/action
%--------------------------------------------------------------------------
p      = spm_LAP_eval(M,qu,qh);
ih     = p.h < 16;
ig     = p.g < 16;
ie     = kron(ones(n,1),ih);
ix     = kron(ones(n,1),ig);
iv     = kron(ones(d,1),ih((1:nv) + ny));
je     = find([ie; ix]); ix(1:nx) = 1;
ju     = find([ix; iv]);
jub    = find([ix; iv; ones(np + nb,1)]);
 
% and other useful indices
%--------------------------------------------------------------------------
ix     = (1:nx);
ih     = (1:nb);
iv     = (1:nv) + nx*n;
ip     = (1:np) + nu;
iup    = (1:(nu + np));
ib     = (1:(np + nb)) + nu;
 
 
 
% Iterate Laplace scheme
%==========================================================================
F      = -Inf;
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
            [E dE] = spm_DEM_eval(M,qu,qp);
            [p dp] = spm_LAP_eval(M,qu,qh);
            
 
            % gradients of log(det(iS)) dDd...
            %==============================================================
            
            % get precision matrices
            %--------------------------------------------------------------
            iSh   = diag(exp(p.h));
            iSg   = diag(exp(p.g));
            iS    = blkdiag(kron(Rh,iSh),kron(Rg,iSg));
 
            
            % gradients of trace(diag(p)) = sum(p); p = precision vector
            %--------------------------------------------------------------
            dpdx  = n*sum(spm_cat({dp.h.dx; dp.g.dx}));
            dpdv  = n*sum(spm_cat({dp.h.dv; dp.g.dv}));
            dpdh  = n*sum(dp.h.dh);
            dpdg  = n*sum(dp.g.dg);
            dpdx  = kron(sparse(1,1,1,1,n),dpdx);
            dpdv  = kron(sparse(1,1,1,1,d),dpdv);
            dDdu  = [dpdx dpdv]';
            dDdh  = [dpdh dpdg]';
 
            
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
            dSdhh = spm_cat({dedhh  [] ;
                             [] dedgg});
                 
            
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
            dLduu = dEdu*dE.du + pu.ic;
            dLdpp = dEdp*dE.dp + pp.ic;
            dLdhh = dSdhh/2    + ph.ic;
            dLdup = dEdu*dE.dp;
            dLdhu = dEdh*dE.du;
            dLduy = dEdu*dE.dy;
            dLduc = dEdu*dE.dc;
            dLdhp = dEdh*dE.dp;
            dLdpu = dLdup';
            dLduh = dLdhu';
            dLdph = dLdhp';
 
 
            % precision and covariances
            %--------------------------------------------------------------                      
            iC    = spm_cat({dLduu dLdup dLduh  ;
                             dLdpu dLdpp dLdph  ;
                             dLdhu dLdhp dLdhh});
            
            dLdbb = spm_sqrtm(iC(ib,ib)*iC(ib,ib)');
            Cup   = spm_inv(iC(iup,iup));
            Chh   = spm_inv(dLdhh);
            
            
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
                dHdh(i) = sum(sum(diCdh.*Cup))/2;
            end
            for i = 1:ng
                Luub    = dE.du'*diSdg{i}*dE.du;
                Lpub    = dE.dp'*diSdg{i}*dE.du;
                Lppb    = dE.dp'*diSdg{i}*dE.dp;
                diCdg   = spm_cat({Luub Lpub';
                                   Lpub Lppb});
                dHdg(i) = sum(sum(diCdg.*Cup))/2;
            end
            
            % parameters
            %--------------------------------------------------------------
            for i = 1:np
                Luup    = dE.dup{i}'*dEdu';
                Lpup    = dEdp*dE.dup{i};
                Luup    = Luup + Luup';
                diCdp   = spm_cat({Luup Lpup';
                                   Lpup [] });
                dHdp(i) = sum(sum(diCdp.*Cup))/2;
            end
            
            % hidden states and causes (disabled for stability)
            %--------------------------------------------------------------
            for i = 1:(nu - nu)
                Lppu    = dE.dpu{i}'*dEdp';
                Lupu    = dEdu*dE.dpu{i};
                Lppu    = Lppu + Lppu';
                diCdu   = spm_cat({[]    Lupu;
                                   Lupu' Lppu});
                dHdu(i) = sum(sum(diCdu.*Cup))/2;
            end
 
            % and concatenate
            %--------------------------------------------------------------
            dHdb  = [dHdh; dHdg];
            dHdb  = [dHdp; dHdb];
            dLdb  = [dLdp; dLdh];

                        
            % save conditional moments (and prediction error) at Q{t}
            %==============================================================
            if iD == 1
                
                % save means
                %----------------------------------------------------------
                Q(is).e   = E;
                Q(is).E   = iS*E;
                Q(is).u   = qu;
                Q(is).p   = qp;
                Q(is).h   = qh;
                
                % and conditional covariances
                %----------------------------------------------------------
                Q(is).u.s = Cup(ix,ix);
                Q(is).u.c = Cup(iv,iv);
                Q(is).p.c = Cup(ip,ip);
                Q(is).h.c = Chh(ih,ih);              
                
                % Free-energy (components)
                %----------------------------------------------------------
                Fc(is,1)  = - E(je)'*iS(je,je)*E(je)/2;
                Fc(is,2)  = - Eu(ju)'*pu.ic(ju,ju)*Eu(ju)/2;
                Fc(is,3)  = - n*ny*log(2*pi)/2;
                Fc(is,4)  = spm_logdet(iS(je,je))/2;
                Fc(is,5)  = spm_logdet(pu.ic(ju,ju)*Cup(ju,ju))/2;

                                
                % Free-action (states and parameters)
                %----------------------------------------------------------
                AC(is)    = sum(Fc(is,:))       ...
                          - Ep'*pp.ic*Ep/2      ...
                          - Eh'*ph.ic*Eh/2      ...
                          + spm_logdet(pp.ic)/2 ...
                          + spm_logdet(ph.ic)/2 ...
                          - spm_logdet(dLdbb)/2;
  
            end
 
            % update conditional moments
            %==============================================================
            
            % prior precision of fluctuations on [hyper] parameters
            %--------------------------------------------------------------
            Kb    = ns*Ib;
            
            % accumulate curvatures of [hyper] parameters
            %--------------------------------------------------------------
            try
                dLdBB = dLdBB*(1 - 1/ns) + dLdbb/ns;

            catch
                dLdBB = dLdbb + Ib*32;
            end

            % whiten gradient (and curvatures) with regularised precision
            %--------------------------------------------------------------
            IB    = diag(diag(dLdBB));
            Cb    = spm_inv(dLdBB + IB*exp(dt));
            dLdb  = Cb*dLdb;
            dHdb  = Cb*dHdb;
            
            % assemble conditional means
            %--------------------------------------------------------------
            q.y  = qu.y(1:n);
            q.x  = qu.x(1:n);
            q.v  = qu.v(1:d);
            q.c  = qu.u(1:d);
            q.p  = qp.p;
            q.h  = qh.h;
            q.g  = qh.g;
            q.d  = dbdt;
                        
            % flow
            %--------------------------------------------------------------
            f.y  =  Dy*spm_vec(q.y)                  ;
            f.u  =  Du*spm_vec(q.x,q.v) - dLdu - dHdu;
            f.c  =  Dc*spm_vec(q.c)                  ;
            f.b  =     spm_vec(q.d)                  ;
            f.d  = -Kb*spm_vec(q.d)     - dLdb - dHdb;

            % and Jacobian
            %--------------------------------------------------------------
            dfdq = {Dy     []        []     []    [] ;
                   -dLduy  Du-dLduu -dLduc  []    [] ;
                    []     []        Dc     []    [] ;
                    []     []        []     []    Ib ;
                    []     []        []    -Ib   -Kb};
          
 
            % update conditional modes of states
            %==============================================================
            dq   = spm_dx(spm_cat(dfdq),spm_vec(f),1/nD);
            q    = spm_unvec(spm_vec(q) + dq,q);
            
            % unpack conditional means
            %--------------------------------------------------------------
            qu.x(1:n) = q.x;
            qu.v(1:d) = q.v;
            qp.p      = q.p;
            qh.h      = q.h;
            qh.g      = q.g;
            dbdt      = q.d;
 
        end % D-Step
 
    end % sequence (ns)
 
    
    % Bayesian parameter averaging
    %======================================================================
 
    % Conditional moments of time-averaged parameters
    %----------------------------------------------------------------------
    Ep    = 0;
    Pp    = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).p.c);
        Ep  = Ep  + P*spm_vec(Q(i).p.p);
        Pp  = Pp + P;       
    end
    Ep    = spm_inv(Pp)*Ep;
    Cp    = spm_inv(Pp + (1 - ns)*pp.ic);
 
    % conditional moments of hyper-parameters
    %----------------------------------------------------------------------
    Eh    = 0;
    Ph    = 0;
    for i = 1:ns
        P   = spm_inv(Q(i).h.c);
        Eh  = Eh  + P*spm_vec({Q(i).h.h Q(i).h.g});
        Ph  = Ph + P;
    end
    Eh    = spm_inv(Ph)*Eh - ph.h;
    Ch    = spm_inv(Ph + (1 - ns)*ph.ic);
 
    % Free-action of states plus free-energy of parameters
    %======================================================================
    FC(1) = sum(Fc(:,1));       % - E'*iS*E/2;
    FC(2) = sum(Fc(:,2));       % - Eu'*pu.ic*Eu/2;
    FC(3) = sum(Fc(:,3));       % - n*ny*log(2*pi)/2;
    FC(4) = sum(Fc(:,4));       %   spm_logdet(iS)/2;
    FC(5) = sum(Fc(:,5));       %   spm_logdet(pu.ic*Cu)/2;
    FC(6) = -Ep'*pp.ic*Ep/2;
    FC(7) = -Eh'*ph.ic*Eh/2;
    FC(8) = spm_logdet(pp.ic*Cp)/2;
    FC(9) = spm_logdet(ph.ic*Ch)/2;
    
    CC(iN,:) = FC;
    S(iN)    = sum(AC);
    Fe       = sum(FC);
 
    % if F is decreasing, revert [hyper] parameters and slow down
    %----------------------------------------------------------------------
    if Fe < F(iN)
        
        % start again if F never increased
        %------------------------------------------------------------------
        if iN == 2 && dt < 8
            DEM.M(1).E.v = DEM.M(1).E.v + 2;
            DEM = spm_LAP(DEM);
            return
        end
                
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
        if Fe - F(iN) < 1e-2 && iN > 4
            convergence = 1; 
        else
            convergence = 0; 
        end
            
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
        dt    = max(dt - 1,-8);
        
    end
 
    % Convergence
    %======================================================================
    if convergence; break, end
    
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
        dF = F(iN) - F(iN - 1);
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
Ep     = spm_inv(Pp)*Ep;
Cp     = spm_inv(Pp + (1 - ns)*pp.ic);
qP.P   = spm_unvec(Up*Ep + pp.p,{M.pE});
qP.C   = Up*Cp*Up';
qP.V   = spm_unvec(diag(qP.C),{M.pE});
qP.U   = Up;
 
% conditional moments of hyper-parameters
%--------------------------------------------------------------------------
Ph = 0;
Eh = 0;
for i = 1:ns
    
    % weight in proportion to precisions
    %----------------------------------------------------------------------
    P  = spm_inv(qH.c{i});
    Eh = Eh + P*qH.p{i};
    Ph = Ph + P;
 
end
Eh     = spm_inv(Ph)*Eh;
Ch     = spm_inv(Ph + (1 - ns)*ph.ic);
P      = spm_unvec(Eh,{qh.h qh.g});
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
 
DEM.F  = F(1:iN);             % [-ve] Free-energy
DEM.S  = S(1:iN);             % [-ve] Free-action
DEM.FC = FC;                  % Free-energy components
DEM.CC = CC;                  % over iterations

return

 
 
