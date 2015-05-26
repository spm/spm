function [B0,BV] = spm_MDP_DP(MDP,OPTION)
% dynamic programming using active inference
% FORMAT [MDP] = spm_MDP_DP(MDP,OPTION,W)
%
% MDP.A(O,N)      - Likelihood of O outcomes given N hidden states
% MDP.B{M}(N,N)   - transition probabilities among hidden states (priors)
% MDP.C(N,1)      - prior preferences (prior over future states)
%
% MDP.V(T - 1,P)  - P allowable policies (control sequences)
%
% OPTION  - {'Free Energy' | 'KL Control' | 'Expected Utility'};
%
% B0      - optimal state action policy or transition matrix
% BV      - corresponding policy using value iteration
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_DP.m 6451 2015-05-26 09:26:03Z karl $

% set up and preliminaries
%==========================================================================

% options
%--------------------------------------------------------------------------
if nargin < 2, OPTION = 'Free Energy'; end

% generative model and initial states
%--------------------------------------------------------------------------
T     = size(MDP.V,1) + 1;        % number of outcomes
Ns    = size(MDP.B{1},1);         % number of hidden states
Nu    = size(MDP.B,2);            % number of hidden controls
p0    = exp(-16);                 % smallest probability

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
    B{j}   = MDP.B{1,j} + p0;
    B{j}   = B{j}*diag(1./sum(B{j}));
    lnB{j} = log(B{j});
end

% terminal probabilities (priors)
%--------------------------------------------------------------------------
try
    C = MDP.C + p0;
    if size(C,2) ~= T
        C = C(:,end)*ones(1,T);
    end
catch
    C = ones(Ns,T);
end
C     = C*diag(1./sum(C));
lnC   = log(C);

% policies, states and their expectations
%--------------------------------------------------------------------------
V     = MDP.V;
Np    = size(V,2);                % number of allowable policies



% policy iteration
%==========================================================================
for s = 1:Ns
    
    
    % Variational iterations (hidden states)
    %======================================================================
    x     = zeros(Ns,T,Np);
    for k = 1:Np
        
        for i = 1:2
            
            % hiddens states (x)
            %--------------------------------------------------------------
            x(s,1,k) = 1;
            
            % future states
            %--------------------------------------------------------------
            for j = 2:(T - 1)
                v = lnB{V(j - 1,k)} *x(:,j - 1,k) + ...
                    lnB{V(j,    k)}'*x(:,j + 1,k);
                x(:,j,k) = spm_softmax(v);
            end
            
            % last state
            %--------------------------------------------------------------
            v = lnB{V(T - 1,k)} *x(:,T - 1,k);
            x(:,T,k) = spm_softmax(v);
            
        end
        
    end
    
    % value of policies (Q)
    %======================================================================
    Q     = zeros(Np,1);
    for k = 1:Np
        
        % path integral of expected free energy
        %------------------------------------------------------------------
        for j = 2:T
            
            switch OPTION
                case{'Free Energy','FE'}
                    v = lnC(:,j) - log(x(:,j,k)) + H;
                    
                case{'KL Control','KL'}
                    v = lnC(:,j) - log(x(:,j,k));
                    
                case{'Expected Utility','EU','RL'}
                    v = lnC(:,j);
                    
                otherwise
                    disp(['unkown option: ' OPTION])
            end
            Q(k)   = Q(k) + v'*x(:,j,k);
            
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
lnC   = lnC(:,end);
g     = 1 - 1/T;
for i = 1:32
    for s  = 1:Ns
        
        % value of actions (Q)
        %------------------------------------------------------------------
        Q     = zeros(Nu,1);
        for k = 1:Nu
            Q(k)   = B{1,k}(:,s)'*(lnC + g*V);
        end
        
        % optimal transition from this state
        %------------------------------------------------------------------
        [u,k]   = max(Q);
        BV(:,s) = B{1,k}(:,s);
        
    end
    
    % optimal transition from this state
    %----------------------------------------------------------------------
    dV   = BV'*(lnC + g*V) - V;
    V    = V + dV;
    
    % convergence
    %----------------------------------------------------------------------
    if norm(dV) < 1e-2, break, end
    
end


