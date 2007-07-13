function [V,X] = spm_DEM_int(M,Z,W)
% Integrates/evaluates a hierarchical model given innovations z{i} and w{i}
% FORMAT [V,X] = spm_DEM_int(M,Z);
%
% M{i}    - model structure
% Z{i}    - innovations (casues)
% W{i}    - innovations (states)
%
% V{i}    - causal states (V{1} = y = response)
% X{i}    - hidden states
%
% The system is evaluated at the prior expectation of the parameters
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% set model indices and missing fields
%--------------------------------------------------------------------------
M    = spm_DEM_M_set(M);

% innovations
%--------------------------------------------------------------------------
try, Z = spm_cat(Z(:)); end
try, W = spm_cat(W(:)); end

% number of states and parameters
%--------------------------------------------------------------------------
nt   = size(Z,2);                        % number of time steps
nl   = size(M,2);                        % number of levels
nv   = sum(spm_vec(M.l));                % number of v (casual states)
nx   = sum(spm_vec(M.n));                % number of x (hidden states)

% estimation parameters (d = m = 1 for static models)
%==========================================================================
dt   = M(1).E.dt;                        % time step
d    = M(1).E.d;                         % order of embedding (causes)
n    = M(1).E.n;                         % order of embedding (states)
nD   = M(1).E.nD;                        % number of iterations per sample
td   = dt/nD;                            % integration time for D-Step

% initialise cell arrays for derivatives z{i} = (d/dt)^i[z], ...
%--------------------------------------------------------------------------
v         = cell(n,1);
x         = cell(n,1);
z         = cell(n,1);
w         = cell(n,1);
[v{:}]    = deal(sparse(nv,1));
[x{:}]    = deal(sparse(nx,1));
[z{:}]    = deal(sparse(nv,1));
[w{:}]    = deal(sparse(nx,1));

% initialise arrays for hierarchical form
%--------------------------------------------------------------------------
vi    = {M.v};
xi    = {M.x};
gi    = {M.v};
fi    = {M.x};
v{1}  = spm_vec(vi);
x{1}  = spm_vec(xi);

dfdvi = cell(nl,nl);
dfdxi = cell(nl,nl);
dgdvi = cell(nl,nl);
dgdxi = cell(nl,nl);
for i = 1:nl
    dgdvi{i,i} = sparse(M(i).l,M(i).l);
    dgdxi{i,i} = sparse(M(i).l,M(i).n);
    dfdvi{i,i} = sparse(M(i).n,M(i).l);
    dfdxi{i,i} = sparse(M(i).n,M(i).n);
end

% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dv        = cell(d,d);
Dx        = cell(n,n);
[Dv{:}]   = deal(sparse(nv,nv));
[Dx{:}]   = deal(sparse(nx,nx));

 
% add constant terms
%--------------------------------------------------------------------------
for i = 2:d
    Dv{i - 1,i} = speye(nv,nv);
end
for i = 2:n
    Dx{i - 1,i} = speye(nx,nx);
end
Dv     = spm_cat(Dv);
Dx     = spm_cat(Dx);
D      = spm_cat(diag({Dv,Dx,Dv,Dx}));


% initialise conditional estimators of states (V and X)
%--------------------------------------------------------------------------
for i = 1:nl
    V{i} = sparse(M(i).l,nt);
    X{i} = sparse(M(i).n,nt);
end


% iterate over sequence (t) and within for static models
%==========================================================================
for t  = 1:nt
    for iD = 1:nD

        % sampling time
        %------------------------------------------------------------------
        ts     = (t + (iD - 1)/nD)*dt;

        % derivatives of innovations
        %------------------------------------------------------------------
        z(1:d) = spm_DEM_embed(Z,d,ts,dt);
        w(1:n) = spm_DEM_embed(W,n,ts,dt);

        % partition states {x,v} into distinct vector arrays v{i}, ...
        %------------------------------------------------------------------
        vi     = spm_unvec(v{1},vi);
        xi     = spm_unvec(x{1},xi);
        zi     = spm_unvec(z{1},vi);
        wi     = spm_unvec(w{1},xi);

        % Derivatives for Jacobian
        %==================================================================
        vi{nl} = zi{nl};
        for  i = (nl - 1):-1:1

            % g(x,v) & f(x,v)
            %--------------------------------------------------------------
            gi{i}          = feval(M(i).g,xi{i},vi{i + 1},M(i).pE);
            fi{i}          = feval(M(i).f,xi{i},vi{i + 1},M(i).pE);
            vi{i}          = gi{i} + zi{i};
            
            % and partial derivatives
            %--------------------------------------------------------------
            dgdxi{i,    i} = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).pE,1);
            dgdvi{i,i + 1} = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).pE,2);
            dfdxi{i,    i} = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).pE,1);
            dfdvi{i,i + 1} = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).pE,2);

        end
        
        % concatenate hierarchical arrays
        %------------------------------------------------------------------
        dgdv = spm_cat(dgdvi);
        dgdx = spm_cat(dgdxi);
        dfdv = spm_cat(dfdvi);
        dfdx = spm_cat(dfdxi);
        
        % update generalised coordinates
        %------------------------------------------------------------------
        v{1}  = spm_vec(vi);
        x{2}  = spm_vec(fi) + w{1};
        for i = 2:(n - 1)
            v{i}     = dgdv*v{i} + dgdx*x{i} + z{i};
            x{i + 1} = dfdv*v{i} + dfdx*x{i} + w{i};
        end
        
        % tensor products for Jabobian
        %------------------------------------------------------------------
        dgdv = kron(eye(d,d),dgdv);
        dgdx = kron(eye(d,n),dgdx);
        dfdv = kron(eye(n,d),dfdv);
        dfdx = kron(eye(n,n),dfdx);
        dfdw = kron(eye(n,n),eye(nx,nx));

        % Jacobian
        %------------------------------------------------------------------
        Jq     = spm_cat({Dv*dgdv Dv*dgdx Dv    []  ;
                          dfdv    dfdx    []    dfdw;
                          []      []      Dv    []  ;
                          []      []      []    Dx});

        q      = {v{1:d} x{1:n} z{1:d} w{1:n}};
        
        % update states q = {x,v,z}
        %------------------------------------------------------------------
        dq     = spm_dx(Jq,D*spm_vec(q),td);

        % and unpack
        %------------------------------------------------------------------
        q      = spm_unvec(spm_vec(q) + dq,q);
        v(1:d)   = q([1:d]);
        x(1:n)   = q([1:n] + d);

        % Save realizarion 
        %==================================================================
        if iD == 1
            for i = 1:nl
                V{i}(:,t) = spm_vec(vi{i});
                if M(i).n
                    X{i}(:,t) = spm_vec(xi{i});
                end
            end
        end
    end
    
end % iterations over t
