function [u,dg,df] = spm_ADEM_diff(M,u)
% Evaluate an active model given innovations z{i} and w{i}
% FORMAT [u dg df] = spm_ADEM_diff(M,u);
%
% M    - generative model
%
% u.a - active states
% u.v - causal states - updated
% u.x - hidden states - updated
% u.z - innovation (causal state)
% u.w - innovation (hidden states)
%
% dg.dv, ...  components of the Jacobian in generalised coordinates
%
% The system is evaluated at the prior expectation of the parameters.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ADEM_diff.m 5962 2014-04-17 12:47:43Z spm $


% number of states and parameters
%--------------------------------------------------------------------------
nl    = size(M,2);                        % number of levels
 
% order parameters (n = 1 for static models)
%==========================================================================
n     = M(1).E.n + 1;                     % order of embedding
 
% initialise arrays for hierarchical form
%--------------------------------------------------------------------------
dfdvi = cell(nl,nl);
dfdxi = cell(nl,nl);
dfdai = cell(nl,nl);
dgdvi = cell(nl,nl);
dgdxi = cell(nl,nl);
dgdai = cell(nl,nl);
for i = 1:nl
    dgdvi{i,i} = sparse(M(i).l,M(i).l);
    dgdxi{i,i} = sparse(M(i).l,M(i).n);
    dgdai{i,i} = sparse(M(i).l,M(i).k);
    dfdvi{i,i} = sparse(M(i).n,M(i).l);
    dfdxi{i,i} = sparse(M(i).n,M(i).n);
    dfdai{i,i} = sparse(M(i).n,M(i).k);
end
 
% partition states {x,v,z,w} into distinct vector arrays v{i}, ...
%--------------------------------------------------------------------------
vi    = spm_unvec(u.v{1},{M.v});
xi    = spm_unvec(u.x{1},{M.x});
ai    = spm_unvec(u.a{1},{M.a});
zi    = spm_unvec(u.z{1},{M.v});


% Derivatives for Jacobian
%==========================================================================
vi{nl} = zi{nl};
for  i = (nl - 1):-1:1
 
    % evaluate
    %----------------------------------------------------------------------
    [dgdx,g] = spm_diff(M(i).g,xi{i},vi{i + 1},ai{i + 1},M(i).pE,1);
    [dfdx,f] = spm_diff(M(i).f,xi{i},vi{i + 1},ai{i + 1},M(i).pE,1);
    dgdv     = spm_diff(M(i).g,xi{i},vi{i + 1},ai{i + 1},M(i).pE,2);
    dfdv     = spm_diff(M(i).f,xi{i},vi{i + 1},ai{i + 1},M(i).pE,2);
    dgda     = spm_diff(M(i).g,xi{i},vi{i + 1},ai{i + 1},M(i).pE,3);
    dfda     = spm_diff(M(i).f,xi{i},vi{i + 1},ai{i + 1},M(i).pE,3);
    
    % g(x,v) & f(x,v)
    %----------------------------------------------------------------------
    gi{i}    = g;
    fi{i}    = f;
    vi{i}    = spm_vec(gi{i}) + spm_vec(zi{i});
    
    % and partial derivatives
    %----------------------------------------------------------------------
    dgdxi{i,    i} = dgdx;
    dgdvi{i,i + 1} = dgdv;
    dgdai{i,i + 1} = dgda;
    dfdxi{i,    i} = dfdx;
    dfdvi{i,i + 1} = dfdv;
    dfdai{i,i + 1} = dfda;
 
end
 
% concatenate hierarchical arrays
%--------------------------------------------------------------------------
dg.da = spm_cat(dgdai);
dg.dv = spm_cat(dgdvi);
dg.dx = spm_cat(dgdxi);
df.da = spm_cat(dfdai);
df.dv = spm_cat(dfdvi);
df.dx = spm_cat(dfdxi);
 
% update generalised coordinates
%--------------------------------------------------------------------------
u.v{1}  = spm_vec(vi);
u.x{2}  = spm_vec(fi) + u.w{1};
for i = 2:(n - 1)
    u.v{i}     = dg.dv*u.v{i} + dg.dx*u.x{i} + dg.da*u.a{i} + u.z{i};
    u.x{i + 1} = df.dv*u.v{i} + df.dx*u.x{i} + df.da*u.a{i} + u.w{i};
end
u.v{n}  = dg.dv*u.v{n} + dg.dx*u.x{n} + dg.da*u.a{n} + u.z{n};

