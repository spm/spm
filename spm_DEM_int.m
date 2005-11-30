function [V,X,E] = spm_DEM_int(M,Z)
% Integrates/evaluates a hierarchical model given some innovations z{i}
% FORMAT [V,X,E] = spm_DEM_int(M,Z);
%
% M{i}    - HDM
% Z{i}    - innovations
%
% V{i}    - causal states (V{1} = y = response)
% X{i}    - hidden states
% E{i}    - innovations (integrated)
%
% The system is evaluated at the prior expectation of the parameters
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$
 
% set model indices and missing fields
%--------------------------------------------------------------------------
M    = spm_M_set(M);
 
% innovations
%--------------------------------------------------------------------------
try
    Z = spm_cat(Z(:));
end
 
 
% number of states and parameters
%--------------------------------------------------------------------------
nt   = size(Z,2);                        % number of time steps
nl   = size(M,2);                        % number of levels
nv   = sum(cat(1,M.l));                  % number of y (casual states)
nx   = sum(cat(1,M.n));                  % number of x (hidden states)
nz   = nv;                               % number of z (innovations)
 
% estimation parameters (d = m = 1 for static models)
%==========================================================================
dt   = M(1).E.dt;                        % time step
n    = M(1).E.d;                         % order of embedding
nq   = 1 + nv + nx + n*nz;               % number of q (augmented states)
 
% initialise cell arrays for derivatives z{i} = (d/dt)^i[z], ...
%---------------------------------------------------------------------------
z         = cell(n + 1,1);
x         = cell(1 + 1,1);
v         = cell(n + 1,1);
[z{:}]    = deal(sparse(nz,1));
[x{:}]    = deal(sparse(nx,1));
[v{:}]    = deal(sparse(nv,1));
q         = {v{1} x{1} z{1:n}}; 
 
% derivatives; dfdv{i,j} = dDi(x)/dDj(v), ...  Di(x) = (d/dt)^i[x], ...
%---------------------------------------------------------------------------
dhdz      = cell(1,n);
dzdz      = cell(n,n);
[dhdz{:}] = deal(sparse(nv,nz));
[dzdz{:}] = deal(sparse(nz,nz));
 
% initialise arrays for hierarchical form
%--------------------------------------------------------------------------
xi     = {M.x};
vi     = {M.v};
fi     = {M.x};
gi     = {M.v};
 
dfdvi  = cell(nl,nl);
dfdxi  = cell(nl,nl);
dgdvi  = cell(nl,nl);
dgdxi  = cell(nl,nl);
 
for i = 1:nl
    dgdvi{i,i} = sparse(M(i).l,M(i).l);
    dgdxi{i,i} = sparse(M(i).l,M(i).n);
    dfdvi{i,i} = sparse(M(i).n,M(i).l);
    dfdxi{i,i} = sparse(M(i).n,M(i).n);
end
 
 
% & add constant terms
%--------------------------------------------------------------------------
for i = 2:n
    dhdz{1,2}     = speye(nv,nz);
    dzdz{i - 1,i} = speye(nz,nz);
end
dhdz  = spm_cat(dhdz);
dzdz  = spm_cat(dzdz);
 
% initialise conditional estimators of states (V and X)
%--------------------------------------------------------------------------
for i = 1:nl
    X{i} = sparse(M(i).n,nt);
    V{i} = sparse(M(i).l,nt);
    E{i} = sparse(M(i).l,nt);
end
 
% initialize hidden states
%--------------------------------------------------------------------------
x{1}   = spm_vec(xi);
 
 
% iterate over sequence (t) and within for static models
%==========================================================================
for t  = 1:nt
 
    % derivatives of innovations
    %----------------------------------------------------------------------
    z(1:n) = spm_DEM_embed(Z,t,n);
 
    % partition states {x,v} into distinct vector arrays v{i}, ...
    %----------------------------------------------------------------------
    vi     = spm_unvec(v{1},vi);
    xi     = spm_unvec(x{1},xi);
    zi     = spm_unvec(z{1},vi);
 
    % Derivatives for Jacobian
    %======================================================================
    vi{nl} = zi{nl};
    for  i = (nl - 1):-1:1
 
        % g(x,v) & f(x,v)
        %------------------------------------------------------------------
        gi{i}          = feval(M(i).g,xi{i},vi{i + 1},M(i).P);
        fi{i}          = feval(M(i).f,xi{i},vi{i + 1},M(i).P);
        vi{i}          = gi{i} + zi{i};
 
        % and partial derivatives
        %------------------------------------------------------------------
        dgdxi{i,    i} = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).P,1);
        dgdvi{i,i + 1} = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).P,2);
        dfdxi{i,    i} = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).P,1);
        dfdvi{i,i + 1} = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).P,2);
 
    end
 
    % concatenate hierarchical arrays
    %----------------------------------------------------------------------
    dgdv = spm_cat(dgdvi);
    dgdx = spm_cat(dgdxi);
    dfdv = spm_cat(dfdvi);
    dfdx = spm_cat(dfdxi);
    dhdv = dgdx*dfdv;
    dhdx = dgdx*dfdx;
 
    v{1} = spm_vec(vi);
    x{2} = spm_vec(fi);
    v{2} = inv(speye(nv) - dgdv)*(dgdx*x{2} + z{2});
 
    % Jacobian
    %----------------------------------------------------------------------
    Jq     = spm_cat({dhdv  dhdx  dhdz  ;
                      dfdv  dfdx  []    ;
                      []    []    dzdz});
 
    dqdt   = spm_vec({v{2} x{2} z{[1:n] + 1}});
 
    % update states q = {x,v,z}
    %----------------------------------------------------------------------
    dq     = spm_dx(Jq,dqdt,dt);
 
    % and unpack
    %----------------------------------------------------------------------
    q      = spm_unvec(spm_vec(q) + dq,q);
    v(1)   = q(1);
    x(1)   = q(2);
    z(1:n) = q([1:n] + 2);
 
    % Save conditional moments
    %======================================================================
 
    % expectations
    %----------------------------------------------------------------------
    for i = 1:nl
        V{i}(:,t) = spm_vec(vi{i});
        if M(i).n
            X{i}(:,t) = spm_vec(xi{i});
        end
        E{i}(:,t) = spm_vec(zi{i});
    end
 
end % iterations over t
