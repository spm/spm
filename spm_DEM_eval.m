function [E,dE,f,g] = spm_DEM_eval(M,qu,qp)
% evaluates state equations and derivatives for DEM schemes
% FORMAT [E dE f g] = spm_DEM_eval(M,qu,qp)
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_eval.m 3113 2009-05-11 15:25:13Z karl $

% persistent variables to avoid redundant evaluations
%==========================================================================
persistent Qp dg df

% test for change in parameters and record them
%--------------------------------------------------------------------------
try
    du = any(spm_vec(qp.p) - spm_vec(Qp.p)) & M(1).E.linear;
catch
    du = 1;
end
Qp   = qp;

%==========================================================================
nl   = size(M,2);                       % number of levels
ne   = sum(spm_vec(M.l));               % number of e (errors)
nv   = sum(spm_vec(M.m));               % number of x (causal states)
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

[fe{:}]   = deal(sparse(nx,1));
[ge{:}]   = deal(sparse(ne,1));

df.dp     = cell(nl - 1,nl - 1);
dg.dp     = cell(nl    ,nl - 1);
for i = 1:(nl - 1)
    dg.dp{i + 1,i} = sparse(M(i).m,M(i).p);
    dg.dp{i    ,i} = sparse(M(i).l,M(i).p);
    df.dp{i    ,i} = sparse(M(i).n,M(i).p);
end

% create deriavtice w.r.t. states if du
%--------------------------------------------------------------------------
if du

    % initialise cell arrays for hierarchical structure
    %--------------------------------------------------------------------------
    df.dv  = cell(nl - 1,nl - 1);
    df.dx  = cell(nl - 1,nl - 1);
    dg.dv  = cell(nl    ,nl - 1);
    dg.dx  = cell(nl    ,nl - 1);

    % & fill in hierarchical forms
    %--------------------------------------------------------------------------
    for i = 1:(nl - 1)
        dg.dv{i + 1,i} = sparse(M(i).m,M(i).m);
        dg.dx{i + 1,i} = sparse(M(i).m,M(i).n);
        dg.dv{i    ,i} = sparse(M(i).l,M(i).m);
        dg.dx{i    ,i} = sparse(M(i).l,M(i).n);
        df.dv{i    ,i} = sparse(M(i).n,M(i).m);
        df.dx{i    ,i} = sparse(M(i).n,M(i).n);
    end

    if np
        dg.dvp      = cell(np,1);
        dg.dxp      = cell(np,1);
        df.dvp      = cell(np,1);
        df.dxp      = cell(np,1);
        [dg.dvp{:}] = deal(dg.dv);
        [dg.dxp{:}] = deal(dg.dx);
        [df.dvp{:}] = deal(df.dv);
        [df.dxp{:}] = deal(df.dx);
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
    try
        dfdp   = h(M(i).fp,xvp{:});
        dgdp   = h(M(i).gp,xvp{:});
        fi     = h(M(i).f, xvp{:});
        gi     = h(M(i).g, xvp{:});

    catch
        [dfdp fi]  = spm_diff(h,M(i).f,xvp{:},4);
        [dgdp gi]  = spm_diff(h,M(i).g,xvp{:},4);
    end

    % and place in array
    %----------------------------------------------------------------------
    g{i,1}     = gi;
    f{i,1}     = fi;
    df.dp{i,i} = dfdp;
    dg.dp{i,i} = dgdp;

    % if the system is nonlinear or the parameters have changed
    %======================================================================
    if du

        % 1st and 2nd partial derivatives (states)
        %------------------------------------------------------------------
        try
            [dgdxp dgdx] = spm_diff(h,M(i).gx,xvp{:},4);
            [dgdvp dgdv] = spm_diff(h,M(i).gv,xvp{:},4);
            [dfdxp dfdx] = spm_diff(h,M(i).fx,xvp{:},4);
            [dfdvp dfdv] = spm_diff(h,M(i).fv,xvp{:},4);
        catch
            [dgdxp dgdx] = spm_diff(h,M(i).g,xvp{:},[2 4]);
            [dgdvp dgdv] = spm_diff(h,M(i).g,xvp{:},[3 4]);
            [dfdxp dfdx] = spm_diff(h,M(i).f,xvp{:},[2 4]);
            [dfdvp dfdv] = spm_diff(h,M(i).f,xvp{:},[3 4]);
        end

        % place 1st derivatives in array
        %------------------------------------------------------------------
        dg.dx{i,i}   = dgdx;
        dg.dv{i,i}   = dgdv;
        df.dx{i,i}   = dfdx;
        df.dv{i,i}   = dfdv;

        % and add constant terms
        %------------------------------------------------------------------
        dg.dv{i + 1,i} = -speye(M(i).m,M(i).m);

        % place 2nd derivatives in array
        %------------------------------------------------------------------
        for j = 1:length(dgdxp)
            dg.dxp{ip}{i,i} = dgdxp{j};
            dg.dvp{ip}{i,i} = dgdvp{j};
            df.dxp{ip}{i,i} = dfdxp{j};
            df.dvp{ip}{i,i} = dfdvp{j};
            ip              = ip + 1;
        end
    end
end

% concatenate hierarchical forms
%--------------------------------------------------------------------------
dgdv  =  spm_cat(dg.dv);
dgdx  =  spm_cat(dg.dx);
dfdv  =  spm_cat(df.dv);
dfdx  =  spm_cat(df.dx);
dfdp  = {spm_cat(df.dp)};
dgdp  = {spm_cat(dg.dp)};
for j = 1:np
    dgdvp{j} = spm_cat(dg.dvp{j});
    dgdxp{j} = spm_cat(dg.dxp{j});
    dfdvp{j} = spm_cat(df.dvp{j});
    dfdxp{j} = spm_cat(df.dxp{j});
end

% prediction errors and states
%==========================================================================
dfdy  =  sparse(nx,ny);
dfdc  =  sparse(nx,nc);
dedy  =  spm_speye(ne,ny);
dedc  = -spm_speye(ne,nc,nc - ne);

% prediction error (E) - causes
%--------------------------------------------------------------------------
for i = 1:n
    y{i} = spm_vec(y{i});
end
ge{1} = [y{1}; v{1}] - [spm_vec(g); u{1}];
for i = 2:n
    ge{i} = dedy*y{i} + dedc*u{i} ...  % generalised response
          - dgdx*x{i} - dgdv*v{i};     % and prediction
end

% prediction error (E) - states
%--------------------------------------------------------------------------
try
    fe{1} = x{2} - spm_vec(f);
end
for i = 2:n - 1
    fe{i} = x{i + 1} ...              % generalised motion
          - dfdx*x{i} - dfdv*v{i};    % and prediction
end

% error
%--------------------------------------------------------------------------
E      =  spm_vec({ge, fe});

% Kronecker forms
%==========================================================================

% dE.dp (parameters)
%--------------------------------------------------------------------------
for i = 2:n
    dgdp{i,1} = sparse(ny + nv,np);
    dfdp{i,1} = sparse(nx,np);
    for p = 1:np
        dgdp{i,1}(:,p) = dgdxp{p}*x{i} + dgdvp{p}*v{i};
        dfdp{i,1}(:,p) = dfdxp{p}*x{i} + dfdvp{p}*v{i};
    end
end

% generalised temporal derivatives: dE.du (states)
%--------------------------------------------------------------------------
dedy  = kron(spm_speye(n,n),dedy);
dedc  = kron(spm_speye(n,d),dedc);
dfdy  = kron(spm_speye(n,n),dfdy);
dfdc  = kron(spm_speye(n,d),dfdc);
dgdx  = kron(spm_speye(n,n),dgdx);
dgdv  = kron(spm_speye(n,d),dgdv);
dfdv  = kron(spm_speye(n,d),dfdv);
dfdx  = kron(spm_speye(n,n),dfdx) - kron(spm_speye(n,n,1),speye(nx,nx));

for k = 1:np
    dgdxp{k} = kron(spm_speye(n,n),dgdxp{k});
    dfdxp{k} = kron(spm_speye(n,n),dfdxp{k});
    dgdvp{k} = kron(spm_speye(n,d),dgdvp{k});
    dfdvp{k} = kron(spm_speye(n,d),dfdvp{k});
end

% 1st error derivatives dE.du (states)
%----------------------------------------------------------------------
dE.dy  =  spm_cat({dedy; dfdy});
dE.dc  =  spm_cat({dedc; dfdc});
dE.dp  = -spm_cat({dgdp; dfdp});
dE.du  = -spm_cat({dgdx, dgdv;
                   dfdx, dfdv});

% 2nd error derivatives
%----------------------------------------------------------------------
for i = 1:np
    dE.dup{i} = -spm_cat({dgdxp{i}, dgdvp{i};
                          dfdxp{i}, dfdvp{i}});
end
if np
    dE.dpu    =  spm_cell_swap(dE.dup);
end


