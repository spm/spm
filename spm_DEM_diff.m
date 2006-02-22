function [qu df dE] = spm_DEM_diff(M,qu,qp)
% evaluates state equations and derivatives for DEM schemes
% FORMAT [qu df dE] = spm_DEM_diff(M,qu,qp)
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
% qu - evaluates:
%  qu.e{i} - generalised errors  (i.e.., g(x,v,P))
%  qu.x(i) - hidden states       (i.e.., f(x,v,P))
%
% df:
%   df.dx  - df/dx
%   df.dv  - df/dv
%   df.du  - [df/dv df/dx]       (u = [v[1:d] x])
%
%
% dE:
%  dE.dv   - de[1:n]/dv[1]
%  dE.du   - de[1:n]/du
%  dE.dy   - de[1:n]/dy[1:n]
%  dE.dc   - de[1:n]/dc[1:d]
%  dE.dp   - de[1:n]/dp
%  dE.dup  - d/dp[de[1:n]/du
%  dE.dpu  - d/du[de[1:n]/dp
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
d    = M(1).E.d;                         % approximation order of q(x,v)
n    = M(1).E.n;                         % embedding order      (n >= d)
w    = length(qu.y) - 1;
 
% derivatives: dedv{i,j} = dDi(e)/dDj(v), ...  Di(e) = (d/dt)^i[e], ...
%--------------------------------------------------------------------------
dfdv      = cell(1,d);
dxdv      = cell(1,d);
dxdx      = cell(1,1);
dxdp      = cell(1,1);
dedx      = cell(n,1);
dedv      = cell(n,d);
dedy      = cell(n,w);
dedc      = cell(n,d);
dedp      = cell(n,1);
[dfdv{:}] = deal(sparse(nx,nv));
[dxdv{:}] = deal(sparse(nx,nv));
[dxdx{:}] = deal(sparse(nx,nx));
[dxdp{:}] = deal(sparse(nx,np));
[dedx{:}] = deal(sparse(ne,nx));
[dedv{:}] = deal(sparse(ne,nv));
[dedy{:}] = deal(sparse(ne,ny));
[dedc{:}] = deal(sparse(ne,nc));
[dedp{:}] = deal(sparse(ne,np));

% derivatives w.r.t. parameters
%--------------------------------------------------------------------------
if np
    dedvp  = cell(np,1);
    dedxp  = cell(np,1);
    dxdvp  = cell(np,1);
    dxdxp  = cell(np,1);
    [dedvp{:}] = deal(dedv);
    [dedxp{:}] = deal(dedx);
    [dxdvp{:}] = deal(dxdv);
    [dxdxp{:}] = deal(dxdx);
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
v     = qu.v;
x     = qu.x;
y     = qu.y;
u     = qu.u;
vi    = spm_unvec(v{1},{M(2:end).v});
xi    = spm_unvec(x{1},{M.x});

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
    dfdxi{i,i} =  spm_diff(h,M(i).f,xvp{:},2);
    dfdvi{i,i} =  spm_diff(h,M(i).f,xvp{:},3);
    dgdpi{i,i} =  spm_diff(h,M(i).g,xvp{:},4);
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
        ip              =  ip + 1;
    end
 
end
 
% concatenate hierarchical forms
%--------------------------------------------------------------------------
e{1}    = [y{1}; v{1}] - [spm_vec(g); u{1}];
x{2}    = spm_vec(f);

dxdx{1} =  speye(nx,nx);
dedv{1} = -spm_cat(dgdvi);
dedx{1} = -spm_cat(dgdxi);
dedp{1} = -spm_cat(dgdpi);
dxdv{2} =  spm_cat(dfdvi);
dxdx{2} =  spm_cat(dfdxi);
dxdp{2} =  spm_cat(dfdpi);

for j = 1:np
    gvp{j}      = -spm_cat(dgdvpi{j});
    gxp{j}      = -spm_cat(dgdxpi{j});
    fvp{j}      =  spm_cat(dfdvpi{j});
    fxp{j}      =  spm_cat(dfdxpi{j});
    dxdvp{j}{2} =  fvp{j};
    dxdxp{j}{2} =  fxp{j};
    dedvp{j}{1} =  gvp{j};
    dedxp{j}{1} =  gxp{j};
end

% 1st-order derivatives
%--------------------------------------------------------------------------
fx    = dxdx{2};
fv    = dxdv{2};
gx    = dedx{1};
gv    = dedv{1};
gy    = dedy{1};
gc    = dedc{1};

% and compute temporal derivatives
%--------------------------------------------------------------------------
for i = 2:M(1).E.n
 
    % 1st-order terms
    %----------------------------------------------------------------------
    x{i + 1}    = fx*x{i} + fv*v{i};
    e{i}        = gx*x{i} + gv*v{i} + gy*y{i} + gc*u{i};
    
    % error derivatives w.r.t. states
    %----------------------------------------------------------------------
    dxdx{i + 1} = fx*dxdx{i};
    dxdv{i + 1} = fx*dxdv{i};
    dedx{i}     = gx*dxdx{i};
    dedv{i}     = gx*dxdv{i};

    % parameters
    %----------------------------------------------------------------------
    for p = 1:np

        % 1st-order terms
        %------------------------------------------------------------------
        dxdp{i + 1}(:,p) = fx*dxdp{i}(:,p) + fxp{p}*x{i} + fvp{p}*v{i};
        dedp{i}(:,p)     = gx*dxdp{i}(:,p) + gxp{p}*x{i} + gvp{p}*v{i};
        
        dxdxp{p}{i + 1}  = fx*dxdxp{p}{i} + fxp{p}*dxdx{i};
        dxdvp{p}{i + 1}  = fx*dxdvp{p}{i} + fxp{p}*dxdv{i};
        dedxp{p}{i}      = gx*dxdxp{p}{i} + gxp{p}*dxdx{i};
        dedvp{p}{i}      = gx*dxdvp{p}{i} + gxp{p}*dxdv{i};
        
    end
    
    % & higher temporal derivatives
    %----------------------------------------------------------------------
    for j = 2:d
        dedv{i,j} = dedv{i - 1,j - 1};
        for k = 1:np
            dedvp{k}{i,j} = dedvp{k}{i - 1,j - 1};
        end
    end
 
end
 
% concatenate over derivatives
%--------------------------------------------------------------------------
qu.e    = e;
qu.x    = x;

% 1st state derivatives
%--------------------------------------------------------------------------
dfdv{1} = fv;
df.dx   = fx;
df.dv   = fv;
df.du   = spm_cat({dfdv fx});

% 1st error derivatives
%--------------------------------------------------------------------------
dE.dx   = dedx(:,1);
dE.dv   = spm_cat(dedv(:,1));
dE.du   = spm_cat({dedv(:,1:d) dedx});
dE.dy   = spm_cat(dedy);
dE.dc   = spm_cat(dedc);
dE.dp   = spm_cat(dedp);

% 2nd error derivatives
%--------------------------------------------------------------------------
for i = 1:np
    dE.dup{i} = spm_cat({dedvp{i}(:,1:d) dedxp{i}});
end
if np
    dE.dpu    = spm_cell_swap(dE.dup);
end
