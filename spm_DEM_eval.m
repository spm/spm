function [E DE] = spm_DEM_eval(M,qu,qp)
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

% persistent variables to avoid redundant evaluations
%==========================================================================
persistent Qp dE dgdv dgdx dfdv dfdx dedy dedc dgdvp dgdxp dfdvp dfdxp
du      = 1;
try, du = any(spm_vec(qp.p) - spm_vec(Qp.p)) & M(1).E.linear; end

%==========================================================================
nl   = size(M,2);                       % number of levels
ne   = sum(spm_vec(M.l));               % number of e (errors)
nv   = sum(spm_vec(M.m));               % number of v (casual states)
nx   = sum(spm_vec(M.n));               % number of x (hidden states)
np   = sum(spm_vec(M.p));               % number of p (parameters)
ny   = M(1).l;                          % number of y (inputs)
nc   = M(end).l;                        % number of c (prior causes)

% order parameters (d = n = 1 for static models)
%==========================================================================
d    = M(1).E.d + 1;                        % generalisation order of q(v)
n    = M(1).E.n + 1;                        % embedding order       (n >= d)

% derivatives: dgdp{i,j} = dDi(e)/dp, ...  Di(e) = (d/dt)^i[e], ...
%--------------------------------------------------------------------------
fe        = cell(n,1);
ge        = cell(n,1);
dfdp      = cell(n,1);
dgdp      = cell(n,1);
dfdpi     = cell(nl - 1,nl - 1);
dgdpi     = cell(nl    ,nl - 1);

[fe{:}]   = deal(sparse(nx,1));
[ge{:}]   = deal(sparse(ne,1));
[dfdp{:}] = deal(sparse(nx,np));
[dgdp{:}] = deal(sparse(ne,np));

for i = 1:(nl - 1)
    dgdpi{i + 1,i} = sparse(M(i).m,M(i).p);
    dgdpi{i    ,i} = sparse(M(i).l,M(i).p);
    dfdpi{i    ,i} = sparse(M(i).n,M(i).p);
end

% create deriavtice w.r.t. states if du
%--------------------------------------------------------------------------
if du
    dfdv      = cell(n,d);
    dfdx      = cell(n,n);
    dgdv      = cell(n,d);
    dgdx      = cell(n,n);
    dedy      = cell(n,n);
    dedc      = cell(n,d);
    dfdy      = cell(n,n);
    dfdc      = cell(n,d);
    [dgdv{:}] = deal(sparse(ne,nv));
    [dgdx{:}] = deal(sparse(ne,nx));
    [dfdv{:}] = deal(sparse(nx,nv));
    [dfdx{:}] = deal(sparse(nx,nx));
    [dedy{:}] = deal(sparse(ne,ny));
    [dedc{:}] = deal(sparse(ne,nc));
    [dfdy{:}] = deal(sparse(nx,ny));
    [dfdc{:}] = deal(sparse(nx,nc));

    % derivatives w.r.t. parameters
    %----------------------------------------------------------------------
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
    dgdvi  = cell(nl    ,nl - 1);
    dgdxi  = cell(nl    ,nl - 1);

    % & fill in hierarchical forms
    %--------------------------------------------------------------------------
    for i = 1:(nl - 1)
        dgdvi{i + 1,i} = sparse(M(i).m,M(i).m);
        dgdxi{i + 1,i} = sparse(M(i).m,M(i).n);
        dgdvi{i    ,i} = sparse(M(i).l,M(i).m);
        dgdxi{i    ,i} = sparse(M(i).l,M(i).n);
        dfdvi{i    ,i} = sparse(M(i).n,M(i).m);
        dfdxi{i    ,i} = sparse(M(i).n,M(i).n);
    end

    if np && du
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
end

% un-concatenate states {v,x} into hierarchical form
%--------------------------------------------------------------------------
v     = qu.v;
x     = qu.x;
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
ip    = 1;
for i = 1:(nl - 1)

    % states and parameters for level i
    %----------------------------------------------------------------------
    xvp        = {xi{i},vi{i},qp.p{i},qp.u{i},M(i).pE};

    % g(x,v), f(x,v) and 1st-order partial derivatives (parameters)
    %----------------------------------------------------------------------
    [dfdpj fj] = spm_diff(h,M(i).f,xvp{:},4);
    [dgdpj gj] = spm_diff(h,M(i).g,xvp{:},4);

    % and place in array
    %----------------------------------------------------------------------
    g{i,1}     = gj;
    f{i,1}     = fj;
    dfdpi{i,i} = dfdpj;
    dgdpi{i,i} = dgdpj;

    % if the system is nonlinear or the parameters have changed
    %======================================================================
    if du

        % 1st and 2nd partial derivatives (states)
        %------------------------------------------------------------------
        [dgdxpj dgdxj] = spm_diff(h,M(i).g,xvp{:},[2 4]);
        [dgdvpj dgdvj] = spm_diff(h,M(i).g,xvp{:},[3 4]);
        [dfdxpj dfdxj] = spm_diff(h,M(i).f,xvp{:},[2 4]);
        [dfdvpj dfdvj] = spm_diff(h,M(i).f,xvp{:},[3 4]);

        % place 1st derivatives in array
        %------------------------------------------------------------------
        dgdxi{i,i} = dgdxj;
        dgdvi{i,i} = dgdvj;
        dfdxi{i,i} = dfdxj;
        dfdvi{i,i} = dfdvj;

        % and add constant terms
        %------------------------------------------------------------------
        dgdvi{i + 1,i} = -speye(M(i).m,M(i).m);

        % place 2nd derivatives in array
        %------------------------------------------------------------------
        for j = 1:M(i).p
            dgdxpi{ip}{i,i} = dgdxpj{j};
            dgdvpi{ip}{i,i} = dgdvpj{j};
            dfdxpi{ip}{i,i} = dfdxpj{j};
            dfdvpi{ip}{i,i} = dfdvpj{j};
            ip              = ip + 1;
        end
    end
end

% concatenate hierarchical forms
%--------------------------------------------------------------------------
if du
    dgdv{1} = spm_cat(dgdvi);
    dgdx{1} = spm_cat(dgdxi);
    dfdv{1} = spm_cat(dfdvi);
    dfdx{1} = spm_cat(dfdxi);
    for j = 1:np
        dgdvp{j}{1} = spm_cat(dgdvpi{j});
        dgdxp{j}{1} = spm_cat(dgdxpi{j});
        dfdvp{j}{1} = spm_cat(dfdvpi{j});
        dfdxp{j}{1} = spm_cat(dfdxpi{j});
    end
end

% prediction error (E) - causes
%--------------------------------------------------------------------------
ge{1} = [y{1}; v{1}] - [spm_vec(g); u{1}];
for i = 2:n
    ge{i} = dedy{1}*y{i} + dedc{1}*u{i} ...  % generalised response
          - dgdx{1}*x{i} - dgdv{1}*v{i};     % and prediction
end

% prediction error (E) - states
%--------------------------------------------------------------------------
try
    fe{1} = x{2} - spm_vec(f);
end
for i = 2:n - 1
    fe{i} = x{i + 1} ...                    % generalised motion
          - dfdx{1}*x{i} - dfdv{1}*v{i};    % and prediction
end

% error
%--------------------------------------------------------------------------
E      =  spm_vec({ge, fe});

% dE.dp (parameters)
%--------------------------------------------------------------------------
dfdp{1} = spm_cat(dfdpi);
dgdp{1} = spm_cat(dgdpi);
for p = 1:np
    for i = 2:n
        dgdp{i}(:,p) = dgdxp{p}{1}*x{i} + dgdvp{p}{1}*v{i};
        dfdp{i}(:,p) = dfdxp{p}{1}*x{i} + dfdvp{p}{1}*v{i};
    end
end
dE.dp  = -spm_cat({dgdp; dfdp});

% generalised temporal derivatives: dE.du (states)
%--------------------------------------------------------------------------
if du
    
    % Kronecker forms
    %----------------------------------------------------------------------
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

    % dE.du (states)
    %----------------------------------------------------------------------
    dE.du  = -spm_cat({dgdx, dgdv;
                       dfdx, dfdv});
    dE.dy  =  spm_cat({dedy; dfdy});
    dE.dc  =  spm_cat({dedc; dfdc});

    % 2nd error derivatives
    %----------------------------------------------------------------------
    for i = 1:np
        dE.dup{i} = -spm_cat({dgdxp{i}, dgdvp{i};
                              dfdxp{i}, dfdvp{i}});
    end
    if np
        dE.dpu    =  spm_cell_swap(dE.dup);
    end
end

% remeber parameters and derivatives
%--------------------------------------------------------------------------
Qp   = qp;
DE   = dE;
