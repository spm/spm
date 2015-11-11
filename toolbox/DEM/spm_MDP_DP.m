function [B0,BV] = spm_MDP_DP(MDP)
% dynamic programming using active inference
% FORMAT [B0,BV] = spm_MDP_DP(MDP)
%
% MDP.A(O,N)      - Likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences (prior over future states)
%
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
%
% B0      - optimal state action policy or transition matrix
% BV      - corresponding policy using value iteration
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_DP.m 6598 2015-11-11 19:48:30Z karl $

% set up and preliminaries
%==========================================================================


% generative model and initial states
%--------------------------------------------------------------------------
T     = size(MDP.V,1) + 1;        % number of outcomes
Ns    = size(MDP.B{1},1);         % number of hidden states
Nu    = size(MDP.B,2);            % number of hidden controls
p0    = exp(-8);                  % smallest probability

% likelihood model (for a partially observed MDP implicit in G)
%--------------------------------------------------------------------------
try
    A = MDP.A + p0;
catch
    A = speye(Ns,Ns) + p0;
end
A     = A*diag(1./sum(A));        % normalise
lnA   = log(A);                   % log probabilities
H     = sum(A.*lnA)';             % negentropy of observations

% transition probabilities (priors)
%--------------------------------------------------------------------------
for j = 1:Nu
    B{j}   = MDP.B{j} + p0;
    B{j}   = B{j}*diag(1./sum(B{j}));
    sB{j}  = B{j};
    rB{j}  = spm_softmax(log(B{j})');
    lnB{j} = log(B{j});
end


% terminal probabilities over outcomes (priors)
%--------------------------------------------------------------------------
try
    C = MDP.C;
catch
    C = zeros(No,1);
end

% asume constant preferences over states
%--------------------------------------------------------------------------
if size(C,2) ~= T
    C = C(:,end)*ones(1,T);
end
C = A'*diag(1./sum(A,2))*spm_softmax(C);
C = log(C);

% policies, states and their expectations
%--------------------------------------------------------------------------
V     = MDP.V;
Np    = size(V,2);                % number of allowable policies

% policy iteration
%==========================================================================
for s = 1:Ns
    
    
    % Variational iterations (hidden states)
    %======================================================================
    x     = zeros(Ns,T,Np) + 1/Ns;
    for k = 1:Np
        
        % gradient descent on free energy
        %------------------------------------------------------------------
        for i = 1:16
            
            % hiddens states (x)
            %--------------------------------------------------------------
            x(:,1,k) = 0;
            x(s,1,k) = 1;
            
            for j = 2:T
                
                % current state
                %----------------------------------------------------------
                xj   = x(:,j,k);
                qx   = log(xj);
                v    = 0;
                
                % evaluate free energy and gradients (v = dFdx)
                %----------------------------------------------------------
                if j > 1, v = v + qx - log(sB{V(j - 1,k)}*x(:,j - 1,k)); end
                if j < T, v = v      - log(rB{V(j    ,k)}*x(:,j + 1,k)); end
                
                % update
                %----------------------------------------------------------
                x(:,j,k) = spm_softmax(qx - v/4);
                F(j,k)   = xj'*v;
                
            end
            
            % convergence
            %--------------------------------------------------------------
            if i > 1
                dF = F0 - sum(F(:,k));
                if dF > 1/128, F0 = F0 - dF; else, break, end
            else
                F0 = sum(F(:,k));
            end
            
        end
    end

    % value of policies (Q)
    %======================================================================
    Q     = zeros(Np,1);
    for k = 1:Np
        
        % path integral of expected free energy
        %------------------------------------------------------------------
        for j = 2:T                 
            v    = C(:,j) - log(x(:,j,k)) + H;
            Q(k) = Q(k) + v'*x(:,j,k);
            
        end
    end
    
    % optimal transition from this state
    %======================================================================
    [u,k]   = max(Q);
    B0(:,s) = lnB{V(1,k)}(:,s);
    
end

if nargout < 2, return, end

% value iteration
%==========================================================================
V     = zeros(Ns,1);
C     = C(:,end);
g     = 1 - 1/T;
for i = 1:32
    for s  = 1:Ns
        
        % value of actions (Q)
        %------------------------------------------------------------------
        Q     = zeros(Nu,1);
        for k = 1:Nu
            Q(k)   = B{1,k}(:,s)'*(C + g*V);
        end
        
        % optimal transition from this state
        %------------------------------------------------------------------
        [u,k]   = max(Q);
        BV(:,s) = B{1,k}(:,s);
        
    end
    
    % optimal transition from this state
    %----------------------------------------------------------------------
    dV   = BV'*(C + g*V) - V;
    V    = V + dV;
    
    % convergence
    %----------------------------------------------------------------------
    if norm(dV) < 1e-2, break, end
    
end


