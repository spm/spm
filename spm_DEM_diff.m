function [qu df dE] = spm_DEM_diff(M,qu,qp)
% evaluates state equations and derivatives for DEM schemes
% FORMAT [qu df dE] = spm_DEM_diff(M,qu,qp)
%
% M  - model structure
% qu - conditional mode of states
%  qu.e{i} - embedded errors e[i]: (i - 1)-th temporal derivative of e(t)
%  qu.v{i} - casual states
%  qu.x(i) - hidden states
%  qu.y(i) - response
%  qu.u(i) - input
% qp - conditional density of parameters
%  qp.p{i} - parameter deviates for i-th level
%  qp.u(i) - basis set
%  qp.x(i) - expansion point (prior expectation)
%
% qu - evaluates:
%  qu.e{i} - embedded errors  (i.e.., g(x,v,P))
%  qu.x(i) - hidden states    (i.e.., f(x,v,P))
% df:
%  df.dv   - df/dv
%  df.dx   - df/dx
% dE:
%  dE.dv   - de[1:n}/dv[1]
%  dE.dV   - de[1:n}/dv[1:d]
%  dE.dx   - de[1:n}/dx[1]
%  dE.dy   - de[1:n}/dy[1:n]
%  dE.dc   - de[1:n}/dc[1:d]
%  dE.dp   - de[1:n}/dp
%  dE.dvp  - d/dp[de[1:n}/dv[1]]
%  dE.dVp  - d/dp[de[1:n}/dv[1:d]]
%  dE.dxp  - d/dp[de[1:n}/dx[1]]
%  dE.dpv  - d/dv[de[1:n}/dp]
%  dE.dpx  - d/dx[de[1:n}/dp]
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
 
% derivatives: dedv{i,j} = dDi(e)/dDj(v), ...  Di(e) = (d/dt)^i[e], ...
%--------------------------------------------------------------------------
dfdv      = cell(1,d);
dfdx      = cell(1,1);
dfdp      = cell(1,1);
dedx      = cell(n,1);
dedv      = cell(n,d);
dedy      = cell(n,n);
dedc      = cell(n,d);
dedp      = cell(n,1);
[dfdv{:}] = deal(sparse(nx,nv));
[dfdx{:}] = deal(sparse(nx,nx));
[dfdp{:}] = deal(sparse(nx,np));
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
    dfdvp  = cell(np,1);
    dfdxp  = cell(np,1);
    [dedvp{:}] = deal(dedv);
    [dedxp{:}] = deal(dedx);
    [dfdvp{:}] = deal(dfdv);
    [dfdxp{:}] = deal(dfdx);
end
 
 
% initialise cell arrays for hierarchical structure
%--------------------------------------------------------------------------
dfdvi  = cell(nl - 1,nl - 1);
dfdxi  = cell(nl - 1,nl - 1);
dfdpi  = cell(nl - 1,nl - 1);
dedvi  = cell(nl    ,nl - 1);
dedxi  = cell(nl    ,nl - 1);
dedpi  = cell(nl    ,nl - 1);
 
% & fill in hierarchical forms
%--------------------------------------------------------------------------
for i = 1:(nl - 1)
    dedvi{i + 1,i} = sparse(M(i).m,M(i).m);
    dedxi{i + 1,i} = sparse(M(i).m,M(i).n);
    dedpi{i + 1,i} = sparse(M(i).m,M(i).p);
    dedvi{i    ,i} = sparse(M(i).l,M(i).m);
    dedxi{i    ,i} = sparse(M(i).l,M(i).n);
    dedpi{i    ,i} = sparse(M(i).l,M(i).p);
    dfdvi{i    ,i} = sparse(M(i).n,M(i).m);
    dfdxi{i    ,i} = sparse(M(i).n,M(i).n);
    dfdpi{i    ,i} = sparse(M(i).n,M(i).p);
end
 
if np
    dedvpi      = cell(np,1);
    dedxpi      = cell(np,1);
    dfdvpi      = cell(np,1);
    dfdxpi      = cell(np,1);
    [dedvpi{:}] = deal(dedvi);
    [dedxpi{:}] = deal(dedxi);
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
k     = 1;
for i = 1:(nl - 1)
 
    % states and parameters for level i
    %----------------------------------------------------------------------
    xvp        = {xi{i},vi{i},qp.p{i},qp.u{i},M(i).pE};
 
    % g(x,v) & f(x,v)
    %----------------------------------------------------------------------
    g{i,1}     = h(M(i).g,xvp{:});
    f{i,1}     = h(M(i).f,xvp{:});
 
    % and partial derivatives
    %----------------------------------------------------------------------
    dedxi{i,i} = -spm_diff(h,M(i).g,xvp{:},2);
    dedvi{i,i} = -spm_diff(h,M(i).g,xvp{:},3);
    dfdxi{i,i} =  spm_diff(h,M(i).f,xvp{:},2);
    dfdvi{i,i} =  spm_diff(h,M(i).f,xvp{:},3);
    dedpi{i,i} = -spm_diff(h,M(i).g,xvp{:},4);
    dfdpi{i,i} =  spm_diff(h,M(i).f,xvp{:},4);
    
    % add constant terms
    %----------------------------------------------------------------------
    dedvi{i + 1,i} =  speye(M(i).m,M(i).m);

    % d/dP[de/dx], d/dP[de/dv],, ...
    %----------------------------------------------------------------------
    dgdxpj     =  spm_diff(h,M(i).g,xvp{:},[2 4]);
    dgdvpj     =  spm_diff(h,M(i).g,xvp{:},[3 4]);
    dfdxpj     =  spm_diff(h,M(i).f,xvp{:},[2 4]);
    dfdvpj     =  spm_diff(h,M(i).f,xvp{:},[3 4]);

    for j = 1:M(i).p
        dedxpi{k}{i,i} = -dgdxpj{j};
        dedvpi{k}{i,i} = -dgdvpj{j};
        dfdxpi{k}{i,i} =  dfdxpj{j};
        dfdvpi{k}{i,i} =  dfdvpj{j};
        k              = k + 1;
    end
 
end
 
% concatenate hierarchical forms
%--------------------------------------------------------------------------
dedv{1} = spm_cat(dedvi);
dedx{1} = spm_cat(dedxi);
dedp{1} = spm_cat(dedpi);
dfdv{1} = spm_cat(dfdvi);
dfdx{1} = spm_cat(dfdxi);
dfdp{1} = spm_cat(dfdpi);
for   j = 1:np
    dedvp{j}{1} = spm_cat(dedvpi{j});
    dedxp{j}{1} = spm_cat(dedxpi{j});
    dfdvp{j}{1} = spm_cat(dfdvpi{j});
    dfdxp{j}{1} = spm_cat(dfdxpi{j});
end
 
% and compute temporal derivatives
%--------------------------------------------------------------------------
e{1}  = [y{1}; v{1}] - [spm_vec(g); u{1}];
x{2}  = spm_vec(f);
for i = 2:M(1).E.n
 
    x{i + 1} = dfdv{1}*v{i} + dfdx{1}*x{i};
    e{i}     = dedy{1}*y{i} + dedc{1}*u{i} + ...
               dedv{1}*v{i} + dedx{1}*x{i};
 
    % error derivatives w.r.t. states
    %----------------------------------------------------------------------
    dedx{i}  = dedx{i - 1}*dfdx{1};
    dedv{i}  = dedx{i - 1}*dfdv{1};
 
    % parameters
    %----------------------------------------------------------------------
    for j = 1:np
        dedp{i}(:,j) = dedxp{j}{1}*x{i} + ...
                       dedvp{j}{1}*v{i} + ...
                       dedx{i - 1}*dfdp{1}(:,j);
        dedxp{j}{i}  = dedxp{j}{i - 1}*dfdx{1} + ...
                       dedx{i - 1}*dfdxp{j}{1};
        dedvp{j}{i}  = dedxp{j}{i - 1}*dfdv{1} + ...
                       dedx{i - 1}*dfdvp{j}{1};
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
qu.e   = e;
qu.x   = x;

df.dx  = dfdx;
df.dv  = dfdv;
dE.dx  = spm_cat(dedx);
dE.dp  = spm_cat(dedp);
dE.dv  = spm_cat(dedv(:,1));
dE.dV  = spm_cat(dedv(:,1:d));
dE.dy  = spm_cat(dedy);
dE.dc  = spm_cat(dedc);
for i = 1:np
    dE.dxp{i} = spm_cat(dedxp{i});
    dE.dvp{i} = spm_cat(dedvp{i}(:,1));
    dE.dVp{i} = spm_cat(dedvp{i}(1:n,1:d));
end
if np
    dE.dpx    = spm_cell_swap(dE.dxp);
    dE.dpv    = spm_cell_swap(dE.dvp);
end


