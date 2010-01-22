function [D] = spm_DEM_eval_diff(x,v,qp,M,bilinear)
% evaluates derivatives for DEM schemes
% FORMAT [D]       = spm_DEM_eval_diff(x,v,qp,M,bilinear)
% v{i} - casual states
% x(i) - hidden states
% qp - conditional density of parameters
%  qp.p{i} - parameter deviates for i-th level
%  qp.u(i) - basis set
%  qp.x(i) - expansion point ( = prior expectation)
% M        - model structure
% bilinear - optional flag to suppress second-order derivatives
%
% D        - derivatives
%   D.dgdv
%   ...
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_eval_diff.m 3695 2010-01-22 14:18:14Z karl $

% check for evaluation of bilinear terms
%--------------------------------------------------------------------------
try
    bilinear;
catch
    bilinear = 1;
end


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
df.dp = cell(nl - 1,nl - 1);
dg.dv = cell(nl    ,nl - 1);
dg.dx = cell(nl    ,nl - 1);
dg.dp = cell(nl    ,nl - 1);

for i = 1:(nl - 1)
    dg.dv{i + 1,i} = sparse(M(i).m,M(i).m);
    dg.dx{i + 1,i} = sparse(M(i).m,M(i).n);
    dg.dp{i + 1,i} = sparse(M(i).m,M(i).p);
    dg.dv{i    ,i} = sparse(M(i).l,M(i).m);
    dg.dx{i    ,i} = sparse(M(i).l,M(i).n);
    dg.dp{i    ,i} = sparse(M(i).l,M(i).p);
    df.dv{i    ,i} = sparse(M(i).n,M(i).m);
    df.dx{i    ,i} = sparse(M(i).n,M(i).n);
    df.dp{i    ,i} = sparse(M(i).n,M(i).p);
end

if bilinear
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
    
    
    % 1st and 2nd partial derivatives (states)
    %----------------------------------------------------------------------
    if bilinear
        try
            [dgdxp dgdx] = spm_diff(h,M(i).gx,xvp{:},4,'q');
            [dgdvp dgdv] = spm_diff(h,M(i).gv,xvp{:},4,'q');
            [dfdxp dfdx] = spm_diff(h,M(i).fx,xvp{:},4,'q');
            [dfdvp dfdv] = spm_diff(h,M(i).fv,xvp{:},4,'q');
        catch
            [dgdxp dgdx] = spm_diff(h,M(i).g,xvp{:},[2 4]);
            [dgdvp dgdv] = spm_diff(h,M(i).g,xvp{:},[3 4]);
            [dfdxp dfdx] = spm_diff(h,M(i).f,xvp{:},[2 4]);
            [dfdvp dfdv] = spm_diff(h,M(i).f,xvp{:},[3 4]);
        end
    else
        try
            dgdx = h(M(i).gx,xvp{:});
            dgdv = h(M(i).gv,xvp{:});
            dfdx = h(M(i).fx,xvp{:});
            dfdv = h(M(i).fv,xvp{:});
        catch
            dgdx = spm_diff(h,M(i).g,xvp{:},2);
            dgdv = spm_diff(h,M(i).g,xvp{:},3);
            dfdx = spm_diff(h,M(i).f,xvp{:},2);
            dfdv = spm_diff(h,M(i).f,xvp{:},3);
        end
    end


    % 1st-order partial derivatives (parameters)
    %----------------------------------------------------------------------
    try
        dfdp   = h(M(i).fp,xvp{:});
        dgdp   = h(M(i).gp,xvp{:});
    catch
        dfdp   = spm_diff(h,M(i).f,xvp{:},4);
        dgdp   = spm_diff(h,M(i).g,xvp{:},4);
    end
    
%     % check which dervatives need to be evaluated
%     %====================================================================
%     D(i).dgdv = nnz(dgdv) + nnz(spm_vec(dgdvp));
%     D(i).dgdx = nnz(dgdx) + nnz(spm_vec(dgdxp));
%     D(i).dfdv = nnz(dfdv) + nnz(spm_vec(dfdvp));
%     D(i).dfdx = nnz(dfdx) + nnz(spm_vec(dfdxp));
    
    
    % Constant terms (linking causes over levels)
    %----------------------------------------------------------------------
    dg.dv{i + 1,i} = -speye(M(i).m,M(i).m);
    
    % place 1st derivatives in array
    %----------------------------------------------------------------------
    dg.dx{i,i} = dgdx;
    dg.dv{i,i} = dgdv;
    df.dx{i,i} = dfdx;
    df.dv{i,i} = dfdv;
    df.dp{i,i} = dfdp;
    dg.dp{i,i} = dgdp;

    % place 2nd derivatives in array
    %----------------------------------------------------------------------
    if bilinear
        for j = 1:length(dgdxp)
            dg.dxp{i}{j}{i,i} = dgdxp{j};
            dg.dvp{i}{j}{i,i} = dgdvp{j};
            df.dxp{i}{j}{i,i} = dfdxp{j};
            df.dvp{i}{j}{i,i} = dfdvp{j};
        end
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

% fixed derivatives w.r.t. prediction errors and states
%--------------------------------------------------------------------------
D.dfdy  =  sparse(nx,ny);
D.dfdc  =  sparse(nx,nc);
D.dedy  =  spm_speye(ne,ny);
D.dedc  = -spm_speye(ne,nc,nc - ne);

% bilinear terms if required
%--------------------------------------------------------------------------
if bilinear
    D.dgdvp = {};
    D.dgdxp = {};
    D.dfdvp = {};
    D.dfdxp = {};
    for i = 1:length(dg.dvp)
        for j = 1:length(dg.dvp{i})
            D.dgdvp{end + 1} = spm_cat(dg.dvp{i}{j});
            D.dgdxp{end + 1} = spm_cat(dg.dxp{i}{j});
            D.dfdvp{end + 1} = spm_cat(df.dvp{i}{j});
            D.dfdxp{end + 1} = spm_cat(df.dxp{i}{j});
        end
    end
end


