function [E dE] = spm_DEM_eval(M,qu,qp)
% evaluates state equations and derivatives for DEM schemes
% FORMAT [E dE] = spm_DEM_eval(M,qu,qp)
%
% M  - model structure
% qu - conditional mode of states
%  qu.v{i} - casual states
%  qu.x(i) - hidden states
%  qu.y(i) - response
%  qu.u(i) - input
% qp - conditional density of parameters
%  qp.p{i} - parameter deviates for i-th level
%  qp.u(i) - basis set
%  qp.x(i) - expansion point ( = prior expectation)
%
% E  - generalised errors  (i.e.., y - g(x,v,P); x[1] - f(x,v,P))
%
% dE:
%  dE.du   - de[1:n]/du
%  dE.dy   - de[1:n]/dy[1:n]
%  dE.dc   - de[1:n]/dc[1:d]
%  dE.dp   - de[1:n]/dp
%  dE.dup  - d/dp[de[1:n]/du
%  dE.dpu  - d/du[de[1:n]/dp
%
% where u = x{1:d]; v[1:d]
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

%==========================================================================
nl   = size(M,2);                       % number of levels
ne   = sum(cat(1,M.l));                 % number of e (errors)
nv   = sum(cat(1,M.m));                 % number of v (casual states)
nx   = sum(cat(1,M.n));                 % number of x (hidden states)
np   = sum(cat(1,M.p));                 % number of p (parameters)
ny   = M(1).l;                          % number of y (inputs)
nc   = M(end).l;                        % number of c (prior causes)

% order parameters (d = n = 1 for static models)
%==========================================================================
d    = M(1).E.d;                        % generalisation order of q(x,v)
n    = M(1).E.n;                        % embedding order       (n >= d)

% derivatives: dgdv{i,j} = dDi(e)/dDj(v), ...  Di(e) = (d/dt)^i[e], ...
%--------------------------------------------------------------------------
fe        = cell(n,1);
ge        = cell(n,1);

dfdv      = cell(n,d);
dfdx      = cell(n,n);
dfdp      = cell(n,1);
dgdv      = cell(n,d);
dgdx      = cell(n,n);
dgdp      = cell(n,1);

dedy      = cell(n,n);
dedc      = cell(n,d);
dfdy      = cell(n,n);
dfdc      = cell(n,d);

[fe{:}]   = deal(sparse(nx,1));
[ge{:}]   = deal(sparse(ne,1));

[dfdv{:}] = deal(sparse(nx,nv));
[dfdx{:}] = deal(sparse(nx,nx));
[dfdp{:}] = deal(sparse(nx,np));
[dgdv{:}] = deal(sparse(ne,nv));
[dgdx{:}] = deal(sparse(ne,nx));
[dgdp{:}] = deal(sparse(ne,np));

[dedy{:}] = deal(sparse(ne,ny));
[dedc{:}] = deal(sparse(ne,nc));
[dfdy{:}] = deal(sparse(nx,ny));
[dfdc{:}] = deal(sparse(nx,nc));

% derivatives w.r.t. parameters
%--------------------------------------------------------------------------
if np
    dgdvp = cell(np,1);
    dgdxp = cell(np,1);
    dfdvp = cell(np,1);
    dfdxp = cell(np,1);
    [dgdvp{:}] = deal(dgdv);
    [dgdxp{:}] = deal(dgdx);
    [dfdvp{:}] = deal(dfdv);
    [dfdxp{:}] = deal(dfdx);
end


% initialise cell arrays for hierarchical structure
%--------------------------------------------------------------------------
dfdvi  = cell(nl - 1,nl - 1);
dfdxi  = cell(nl - 1,nl - 1);
dfdpi  = cell(nl - 1,nl - 1);
dgdvi  = cell(nl    ,nl - 1);
dgdxi  = cell(nl    ,nl - 1);
dgdpi  = cell(nl    ,nl - 1);

% & fill in hierarchical forms
%--------------------------------------------------------------------------
for i = 1:(nl - 1)
    dgdvi{i + 1,i} = sparse(M(i).m,M(i).m);
    dgdxi{i + 1,i} = sparse(M(i).m,M(i).n);
    dgdpi{i + 1,i} = sparse(M(i).m,M(i).p);
    dgdvi{i    ,i} = sparse(M(i).l,M(i).m);
    dgdxi{i    ,i} = sparse(M(i).l,M(i).n);
    dgdpi{i    ,i} = sparse(M(i).l,M(i).p);
    dfdvi{i    ,i} = sparse(M(i).n,M(i).m);
    dfdxi{i    ,i} = sparse(M(i).n,M(i).n);
    dfdpi{i    ,i} = sparse(M(i).n,M(i).p);
end

if np
    dgdvpi      = cell(np,1);
    dgdxpi      = cell(np,1);
    dfdvpi      = cell(np,1);
    dfdxpi      = cell(np,1);
    [dgdvpi{:}] = deal(dgdvi);
    [dgdxpi{:}] = deal(dgdxi);
    [dfdvpi{:}] = deal(dfdvi);
    [dfdxpi{:}] = deal(dfdxi);
end

% & add constant terms
%--------------------------------------------------------------------------
for i = 1:n
    dedy{i,i}   =  speye(ne,ny);
end
for i = 1:d
    dedc{i,i}   = -flipdim(flipdim(speye(ne,nc),1),2);
end


% un-concatenate states {v,x} into hierarchical form
%--------------------------------------------------------------------------
x     = qu.x;
v     = qu.v;
y     = qu.y;
u     = qu.u;
vi    = spm_unvec(v{1},{M(1 + 1:end).v});
xi    = spm_unvec(x{1},{M(1:end - 1).x});

% inline function for evaluating projected parameters
%--------------------------------------------------------------------------
h     = 'feval(f,x,v,spm_unvec(spm_vec(p) + u*q,p))';
h     = inline(h,'f','x','v','q','u','p');


% Derivatives at each hierarchical level
%==========================================================================
ix    = 1;
iv    = 1;
ip    = 1;
for i = 1:(nl - 1)

    % states and parameters for level i
    %----------------------------------------------------------------------
    xvp        = {xi{i},vi{i},qp.p{i},qp.u{i},M(i).pE};

    % g(x,v) & f(x,v)
    %----------------------------------------------------------------------
    g{i,1}     = h(M(i).g,xvp{:});
    f{i,1}     = h(M(i).f,xvp{:});

    % 1st-order partial derivatives
    %----------------------------------------------------------------------
    dgdxi{i,i} =  spm_diff(h,M(i).g,xvp{:},2);
    dgdvi{i,i} =  spm_diff(h,M(i).g,xvp{:},3);
    dgdpi{i,i} =  spm_diff(h,M(i).g,xvp{:},4);
    dfdxi{i,i} =  spm_diff(h,M(i).f,xvp{:},2);
    dfdvi{i,i} =  spm_diff(h,M(i).f,xvp{:},3);
    dfdpi{i,i} =  spm_diff(h,M(i).f,xvp{:},4);

    % add constant terms
    %----------------------------------------------------------------------
    dgdvi{i + 1,i} = -speye(M(i).m,M(i).m);

    % d/dP[de/dx], d/dP[de/dv],, ...
    %----------------------------------------------------------------------
    dgdxpj     =  spm_diff(h,M(i).g,xvp{:},[2 4]);
    dgdvpj     =  spm_diff(h,M(i).g,xvp{:},[3 4]);
    dfdxpj     =  spm_diff(h,M(i).f,xvp{:},[2 4]);
    dfdvpj     =  spm_diff(h,M(i).f,xvp{:},[3 4]);

    for j = 1:M(i).p
        dgdxpi{ip}{i,i} = dgdxpj{j};
        dgdvpi{ip}{i,i} = dgdvpj{j};
        dfdxpi{ip}{i,i} = dfdxpj{j};
        dfdvpi{ip}{i,i} = dfdvpj{j};
        ip              = ip + 1;
    end

end

% concatenate hierarchical forms
%--------------------------------------------------------------------------
dgdv{1} = spm_cat(dgdvi);
dgdx{1} = spm_cat(dgdxi);
dgdp{1} = spm_cat(dgdpi);
dfdv{1} = spm_cat(dfdvi);
dfdx{1} = spm_cat(dfdxi);
dfdp{1} = spm_cat(dfdpi);

for j = 1:np
    dgdvp{j}{1} = spm_cat(dgdvpi{j});
    dgdxp{j}{1} = spm_cat(dgdxpi{j});
    dfdvp{j}{1} = spm_cat(dfdvpi{j});
    dfdxp{j}{1} = spm_cat(dfdxpi{j});
end


% prediction error (E) - response
%--------------------------------------------------------------------------
ge{1} = [y{1}; v{1}] - [spm_vec(g); u{1}];
for i = 2:n
    ge{i} = dedy{1}*y{i} + dedc{1}*u{i} ...  % generalised response
          - dgdx{1}*x{i} - dgdv{1}*v{i};     % and prediction
end

% prediction error (E) - states
%--------------------------------------------------------------------------
fe{1} = x{2} - spm_vec(f);
for i = 2:n - 1
    fe{i} = x{i + 1} ...                    % generalised motion
          - dfdx{1}*x{i} - dfdv{1}*v{i};    % and prediction
end

% generalised temporal derivatives (dE.du) (states)
%--------------------------------------------------------------------------
for i = 2:n
    dgdx{i,i} = dgdx{1};
    dfdx{i,i} = dfdx{1};
    for k = 1:np
        dgdxp{k}{i,i} = dgdxp{k}{1};
        dfdxp{k}{i,i} = dfdxp{k}{1};
    end
    dfdx{i - 1,i}     = -speye(nx,nx);
end
for i = 2:d
    dgdv{i,i} = dgdv{1};
    dfdv{i,i} = dfdv{1};
    for k = 1:np
        dgdvp{k}{i,i} = dgdvp{k}{1};
        dfdvp{k}{i,i} = dfdvp{k}{1};
    end
end

% dE.dp (parameters)
%--------------------------------------------------------------------------
for p = 1:np
    for i = 2:n
        dgdp{i}(:,p) = dgdxp{p}{1}*x{i} + dgdvp{p}{1}*v{i};
        dfdp{i}(:,p) = dfdxp{p}{1}*x{i} + dfdvp{p}{1}*v{i};
    end
end


% error
%--------------------------------------------------------------------------
E      =  spm_vec({ge,   fe});

% 1st error derivatives
%--------------------------------------------------------------------------
dE.du  = -spm_cat({dgdx, dgdv;
                   dfdx, dfdv});
dE.dy  =  spm_cat({dedy; dfdy});
dE.dc  =  spm_cat({dedc; dfdc});
dE.dp  = -spm_cat({dgdp; dfdp});

% 2nd error derivatives
%--------------------------------------------------------------------------
for i = 1:np
    dE.dup{i} = -spm_cat({dgdxp{i}, dgdvp{i};
                          dfdxp{i}, dfdvp{i}});
end
if np
    dE.dpu    =  spm_cell_swap(dE.dup);
end

