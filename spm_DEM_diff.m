function [u, dgdv, dgdx, dfdv, dfdx] = spm_DEM_diff(M,u);
% evaluates a hierarchical model given innovations z{i} and w{i}
% FORMAT [u dgdv dgdx dfdv dfdx] = spm_DEM_diff(M,u);
%
% M    - generative model
%
% u.v - causal states
% u.x - hidden states
% u.z - innovation (causal state)
% u.w - innovation (hidden states)
%
% dgdv, ...  components of the Jacobian in generalised coordinates
%
% The system is evaluated at the prior expectation of the parameters
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id$
 
% number of states and parameters
%--------------------------------------------------------------------------
nl    = size(M,2);                        % number of levels
nv    = sum(spm_vec(M.l));                % number of v (casual states)
nx    = sum(spm_vec(M.n));                % number of x (hidden states)
 
% order parameters (n = 1 for static models)
%==========================================================================
n     = M(1).E.n + 1;                     % order of embedding
 
% initialise arrays for hierarchical form
%--------------------------------------------------------------------------
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
 
% partition states {x,v,z,w} into distinct vector arrays v{i}, ...
%------------------------------------------------------------------
vi    = spm_unvec(u.v{1},{M.v});
xi    = spm_unvec(u.x{1},{M.x});
zi    = spm_unvec(u.z{1},{M.v});
wi    = spm_unvec(u.w{1},{M.x});
 
% Derivatives for Jacobian
%==================================================================
vi{nl} = zi{nl};
for  i = (nl - 1):-1:1
 
    % evaluate
    %--------------------------------------------------------------
    [dgdx g] = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).pE,1);
    [dfdx f] = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).pE,1);
    dgdv     = spm_diff(M(i).g,xi{i},vi{i + 1},M(i).pE,2);
    dfdv     = spm_diff(M(i).f,xi{i},vi{i + 1},M(i).pE,2);
    
    % g(x,v) & f(x,v)
    %--------------------------------------------------------------
    gi{i}    = g;
    fi{i}    = f;
    vi{i}    = gi{i} + zi{i};
    
    % and partial derivatives
    %--------------------------------------------------------------
    dgdxi{i,    i} = dgdx;
    dgdvi{i,i + 1} = dgdv;
    dfdxi{i,    i} = dfdx;
    dfdvi{i,i + 1} = dfdv;
 
end
 
% concatenate hierarchical arrays
%------------------------------------------------------------------
dgdv = spm_cat(dgdvi);
dgdx = spm_cat(dgdxi);
dfdv = spm_cat(dfdvi);
dfdx = spm_cat(dfdxi);
 
% update generalised coordinates
%------------------------------------------------------------------
u.v{1}  = spm_vec(vi);
u.x{2}  = spm_vec(fi) + u.w{1};
for i = 2:(n - 1)
    u.v{i}     = dgdv*u.v{i} + dgdx*u.x{i} + u.z{i};
    u.x{i + 1} = dfdv*u.v{i} + dfdx*u.x{i} + u.w{i};
end

