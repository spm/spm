function [V,X] = spm_DEM_int(M,Z,W)
% Integrates/evaluates a hierarchical model given innovations z{i} and w{i}
% FORMAT [V,X] = spm_DEM_int(M,Z);
%
% M{i}    - model structure
% Z{i}    - innovations (causes)
% W{i}    - innovations (states)
%
% V{i}    - causal states (V{1} = y = response)
% X{i}    - hidden states
%
% The system is evaluated at the prior expectation of the parameters
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id$
 
% set model indices and missing fields
%--------------------------------------------------------------------------
M    = spm_DEM_M_set(M);
 
% innovations
%--------------------------------------------------------------------------
try, Z = spm_cat(Z(:)); end
try, W = spm_cat(W(:)); end
 
% number of states and parameters
%--------------------------------------------------------------------------
nt   = size(Z,2);                        % number of time steps
nl   = size(M,2);                        % number of levels
nv   = sum(spm_vec(M.l));                % number of v (casual states)
nx   = sum(spm_vec(M.n));                % number of x (hidden states)
 
% order parameters (n= 1 for static models)
%==========================================================================
dt   = M(1).E.dt;                        % time step
n    = M(1).E.n + 1;                     % order of embedding
nD   = M(1).E.nD;                        % number of iterations per sample
td   = dt/nD;                            % integration time for D-Step
 
% initialise cell arrays for derivatives z{i} = (d/dt)^i[z], ...
%--------------------------------------------------------------------------
v      = cell(n,1);
x      = cell(n,1);
z      = cell(n,1);
w      = cell(n,1);
[v{:}] = deal(sparse(nv,1));
[x{:}] = deal(sparse(nx,1));
[z{:}] = deal(sparse(nv,1));
[w{:}] = deal(sparse(nx,1));
 
% initialise with starting conditions
%--------------------------------------------------------------------------
v{1}   = spm_vec({M.v});
x{1}   = spm_vec({M.x});
u.v    = v;
u.x    = x;
u.z    = z;
u.w    = w;
 
% derivatives for Jacobian of D-step
%--------------------------------------------------------------------------
Dx    = kron(spm_speye(n,n,1),spm_speye(nx,nx,0));
Dv    = kron(spm_speye(n,n,1),spm_speye(nv,nv,0));
D     = spm_cat(diag({Dv,Dx,Dv,Dx}));
dfdw  = kron(eye(n,n),eye(nx,nx));
 
% initialise conditional estimators of states to be saved (V and X)
%--------------------------------------------------------------------------
for i = 1:nl
    V{i} = sparse(M(i).l,nt);
    X{i} = sparse(M(i).n,nt);
end
 
 
% iterate over sequence (t) and within for static models
%==========================================================================
for t  = 1:nt
    for iD = 1:nD
 
        % sampling time
        %------------------------------------------------------------------
        ts  = (t + (iD - 1)/nD)*dt;
 
        % derivatives of innovations
        %------------------------------------------------------------------
        u.z = spm_DEM_embed(Z,n,ts,dt);
        u.w = spm_DEM_embed(W,n,ts,dt);
 
        % evaluate
        %------------------------------------------------------------------ 
        [u dgdv dgdx dfdv dfdx] = spm_DEM_diff(M,u);
 
        % tensor products for Jabobian
        %------------------------------------------------------------------
        dgdv = kron(spm_speye(n,n,1),dgdv);
        dgdx = kron(spm_speye(n,n,1),dgdx);
        dfdv = kron(spm_speye(n,n,0),dfdv);
        dfdx = kron(spm_speye(n,n,0),dfdx);

        % Save realisation 
        %==================================================================
        vi     = spm_unvec(u.v{1},{M.v});
        xi     = spm_unvec(u.x{1},{M.x});
        if iD == 1
            for i = 1:nl
                V{i}(:,t) = spm_vec(vi{i});
                if M(i).n
                    X{i}(:,t) = spm_vec(xi{i});
                end
            end
        end       
        
        % Jacobian for update
        %------------------------------------------------------------------
        J      = spm_cat({dgdv dgdx Dv  []  ;
                          dfdv dfdx []  dfdw;
                          []   []   Dv  []  ;
                          []   []   []  Dx});
        
        % update states u = {x,v,z,w}
        %------------------------------------------------------------------
        du     = spm_dx(J,D*spm_vec(u),td);
 
        % and unpack
        %------------------------------------------------------------------
        u      = spm_unvec(spm_vec(u) + du,u);
 
    end % iterations over iD
    
end % iterations over t
