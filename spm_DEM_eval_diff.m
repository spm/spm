function [D] = spm_DEM_eval_diff(x,v,qp,M)
% evaluates derivatives for DEM schemes
% FORMAT [D] = spm_DEM_eval_diff(x,v,qp,M)
%
% v{i} - casual states
% x(i) - hidden states
% qp - conditional density of parameters
%  qp.p{i} - parameter deviates for i-th level
%  qp.u(i) - basis set
%  qp.x(i) - expansion point ( = prior expectation)
% M  - model structure
%
% D.dgdv  - derivatives
% ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_eval_diff.m 3655 2009-12-23 20:15:34Z karl $

% get dimensions
%==========================================================================
nl    = size(M,2);                       % number of levels
ne    = sum(spm_vec(M.l));               % number of e (errors)
nx    = sum(spm_vec(M.n));               % number of x (hidden states)
ny    = M(1).l;                          % number of y (inputs)
nc    = M(end).l;                        % number of c (prior causes)


% initialise cell arrays for hierarchical structure
%--------------------------------------------------------------------------
df.dv = cell(nl - 1,nl - 1);
df.dx = cell(nl - 1,nl - 1);
dg.dv = cell(nl    ,nl - 1);
dg.dx = cell(nl    ,nl - 1);
df.dp = cell(nl - 1,nl - 1);
dg.dp = cell(nl    ,nl - 1);
for i = 1:(nl - 1)
    dg.dp{i + 1,i} = sparse(M(i).m,M(i).p);
    dg.dp{i    ,i} = sparse(M(i).l,M(i).p);
    df.dp{i    ,i} = sparse(M(i).n,M(i).p);
end


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

for i = 1:(nl - 1)
    dg.dvp{i}      = cell(M(i).p,1);
    dg.dxp{i}      = cell(M(i).p,1);
    df.dvp{i}      = cell(M(i).p,1);
    df.dxp{i}      = cell(M(i).p,1);
    [dg.dvp{i}{:}] = deal(dg.dv);
    [dg.dxp{i}{:}] = deal(dg.dx);
    [df.dvp{i}{:}] = deal(df.dv);
    [df.dxp{i}{:}] = deal(df.dx);
end


% Derivatives at each hierarchical level
%==========================================================================

% inline function for evaluating projected parameters
%--------------------------------------------------------------------------
h     = 'feval(f,x,v,spm_unvec(spm_vec(p) + u*q,p))';
h     = inline(h,'f','x','v','q','u','p');
for i = 1:(nl - 1)

    % states level i
    %----------------------------------------------------------------------
    xvp = {x{i},v{i},qp.p{i},qp.u{i},M(i).pE};

    % Constant terms (linking causes over levels)
    %----------------------------------------------------------------------
    dg.dv{i + 1,i} = -speye(M(i).m,M(i).m);

    % 1st and 2nd partial derivatives (states)
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    dg.dx{i,i} = dgdx;
    dg.dv{i,i} = dgdv;
    df.dx{i,i} = dfdx;
    df.dv{i,i} = dfdv;

    % 1st-order partial derivatives (parameters)
    %----------------------------------------------------------------------
    try
        dfdp   = h(M(i).fp,xvp{:});
        dgdp   = h(M(i).gp,xvp{:});
    catch
        dfdp   = spm_diff(h,M(i).f,xvp{:},4);
        dgdp   = spm_diff(h,M(i).g,xvp{:},4);
    end

    % and place in array
    %----------------------------------------------------------------------
    df.dp{i,i} = dfdp;
    dg.dp{i,i} = dgdp;

    % place 2nd derivatives in array
    %----------------------------------------------------------------------
    for j = 1:length(dgdxp)
        dg.dxp{i}{j}{i,i} = dgdxp{j};
        dg.dvp{i}{j}{i,i} = dgdvp{j};
        df.dxp{i}{j}{i,i} = dfdxp{j};
        df.dvp{i}{j}{i,i} = dfdvp{j};
    end
end

% concatenate hierarchical forms
%==========================================================================
D.dgdv  = spm_cat(dg.dv);
D.dgdx  = spm_cat(dg.dx);
D.dfdv  = spm_cat(df.dv);
D.dfdx  = spm_cat(df.dx);
D.dfdp  = spm_cat(df.dp);
D.dgdp  = spm_cat(dg.dp);

D.dgdvp = {};
D.dgdxp = {};
D.dfdvp = {};
D.dfdxp = {};

k     = 1;
for i = 1:length(dg.dvp)
    for j = 1:length(dg.dvp{i})
        D.dgdvp{k} = spm_cat(dg.dvp{i}{j});
        D.dgdxp{k} = spm_cat(dg.dxp{i}{j});
        D.dfdvp{k} = spm_cat(df.dvp{i}{j});
        D.dfdxp{k} = spm_cat(df.dxp{i}{j});
        k = k + 1;
    end
end

% fixed derivatives w.r.t. prediction errors and states
%--------------------------------------------------------------------------
D.dfdy  =  sparse(nx,ny);
D.dfdc  =  sparse(nx,nc);
D.dedy  =  spm_speye(ne,ny);
D.dedc  = -spm_speye(ne,nc,nc - ne);
