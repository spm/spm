function [V,X,Z,W] = spm_DEM_int(M,z,w,c)
% Integrates/evaluates a hierarchical model given innovations z{i} and w{i}
% FORMAT [V,X,Z,W] = spm_DEM_int(M,z,w,c);
%
% M{i}    - model structure
% z{i}    - innovations (causes)
% w{i}    - innovations (states)
% c{i}    - exogenous causes
%
% V{i}    - causal states (V{1} = y = response)
% X{i}    - hidden states
% Z{i}    - fluctuations (causes)
% W{i}    - fluctuations (states)
%
% The system is evaluated at the prior expectation of the parameters
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_int.m 3715 2010-02-08 13:57:26Z karl $

% set model indices and missing fields
%--------------------------------------------------------------------------
M      = spm_DEM_M_set(M);

% innovations
%--------------------------------------------------------------------------
try, z = spm_cat(z(:)); end
try, w = spm_cat(w(:)); end
try, c = spm_cat(c(:)); end

% number of states and parameters
%--------------------------------------------------------------------------
nt   = size(z,2);                        % number of time steps
nl   = size(M,2);                        % number of levels
nv   = sum(spm_vec(M.l));                % number of v (casual states)
nx   = sum(spm_vec(M.n));                % number of x (hidden states)

% order parameters (n= 1 for static models)
%==========================================================================
dt   = M(1).E.dt;                        % time step
n    = M(1).E.n + 1;                     % order of embedding
nD   = M(1).E.nD;                        % number of iterations per sample
td   = dt/nD;                            % integration time for D-Step

% initialize cell arrays for derivatives z{i} = (d/dt)^i[z], ...
%--------------------------------------------------------------------------
u.v      = cell(n,1);
u.x      = cell(n,1);
u.z      = cell(n,1);
u.w      = cell(n,1);
[u.v{:}] = deal(sparse(nv,1));
[u.x{:}] = deal(sparse(nx,1));
[u.z{:}] = deal(sparse(nv,1));
[u.w{:}] = deal(sparse(nx,1));

% hyperparameters
%--------------------------------------------------------------------------
ph.h   = {M.hE};
ph.g   = {M.gE};

% initialize with starting conditions
%--------------------------------------------------------------------------
vi     = {M.v};
xi     = {M.x};
u.v{1} = spm_vec(vi);
u.x{1} = spm_vec(xi);


% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dx    = kron(spm_speye(n,n,1),spm_speye(nx,nx,0));
Dv    = kron(spm_speye(n,n,1),spm_speye(nv,nv,0));
D     = spm_cat(spm_diag({Dv,Dx,Dv,Dx}));
dfdw  = kron(eye(n,n),eye(nx,nx));

% initialize conditional estimators of states to be saved (V and X)
%--------------------------------------------------------------------------
for i = 1:nl
    V{i} = sparse(M(i).l,nt);
    X{i} = sparse(M(i).n,nt);
    Z{i} = sparse(M(i).l,nt);
    W{i} = sparse(M(i).n,nt);
end

% method for state-dependent precision
%------------------------------------------------------------------
if isempty(spm_vec(M.ph)) && isempty(spm_vec(M.pg))
    state_dependent = 0;
    Sz = 1;
    Sw = 1;
else
    state_dependent = 1;
end



% iterate over sequence (t) and within for static models
%==========================================================================
for t  = 1:nt
    for iD = 1:nD


        % Get generalised motion of random fluctuations
        %==================================================================

        % sampling time
        %------------------------------------------------------------------
        ts   = (t + (iD - 1)/nD)*dt;

        % evaluate state-dependent precision
        %------------------------------------------------------------------
        if state_dependent
            pu.x = {spm_vec(xi(1:end - 1))};
            pu.v = {spm_vec(vi(1 + 1:end))};
            p    = spm_LAP_eval(M,pu,ph);
            Sz   = sparse(diag(exp(-p.h/2)));
            Sw   = sparse(diag(exp(-p.g/2)));
        end

        % derivatives of innovations (and exogenous input)
        %------------------------------------------------------------------
        u.z  = spm_DEM_embed(Sz*z + c,n,ts,dt);
        u.w  = spm_DEM_embed(Sw*w,    n,ts,dt);


        % Evaluate and update states
        %==================================================================

        % evaluate functions
        %------------------------------------------------------------------
        [u dg df] = spm_DEM_diff(M,u);

        % tensor products for Jacobian
        %------------------------------------------------------------------
        dgdv = kron(spm_speye(n,n,1),dg.dv);
        dgdx = kron(spm_speye(n,n,1),dg.dx);
        dfdv = kron(spm_speye(n,n,0),df.dv);
        dfdx = kron(spm_speye(n,n,0),df.dx);

        % Save realization
        %==================================================================
        vi   = spm_unvec(u.v{1},{M.v});
        xi   = spm_unvec(u.x{1},{M.x});
        zi   = spm_unvec(u.z{1},{M.v});
        wi   = spm_unvec(u.w{1},{M.x});
        if iD == 1
            for i = 1:nl
                if M(i).l, V{i}(:,t) = spm_vec(vi{i}); end
                if M(i).n, X{i}(:,t) = spm_vec(xi{i}); end
                if M(i).l, Z{i}(:,t) = spm_vec(zi{i}); end
                if M(i).n, W{i}(:,t) = spm_vec(wi{i}); end
            end
        end


        % Jacobian for update
        %------------------------------------------------------------------
        J      = spm_cat({dgdv dgdx Dv  []  ;
                          dfdv dfdx []  dfdw;
                          []   []   Dv  []  ;
                          []   []   []  Dx});

        % update states u = {x,v,z,w}
        %------------------------------------------------------------------
        du     = spm_dx(J,D*spm_vec(u),td);

        % and unpack
        %------------------------------------------------------------------
        u      = spm_unvec(spm_vec(u) + du,u);

    end % iterations over iD

end % iterations over t
