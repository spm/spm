function spm_MDP(X,V,T)
% Demo for Discrete Markov Decision process (planning)
% FORMAT spm_MDP(X,V,T)
% X    - initial and goal position
% V    - initial and goal velocity
% T    - number of time-steps
%
% P.a  - 0th order coefficients of force
% P.b  - 1st order coefficients of force
% P.c  – 2nd order coefficients of force
% P.d  - action coefficient
%
% M    - model structure
%
% f    - flow dx/dt
%
% This routine uses a Markov decisions process formulation of the mountain
% car problem to illustrate prospective free energy minimization under a
% variational Bayesian learning scheme. The key notion here is that the
% agent represents future states and action (in a pullback sense), where it
% has strong prior beliefs about future states. The intervening states and
% actions are optimized with respect to current sensory data to provide
% predictions about the next sensory state, which action fulfills. The
% result is a planned trajectory through state space that realizes prior
% beliefs in a prospective sense.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP.m 4339 2011-06-02 13:04:49Z karl $
 
% set up and preliminaries
%==========================================================================
TOL  = 64;                                  % tolerance for log
u    = [-2 -1 0 1 2];                       % levels of control
Nu   = length(u);                           % number of control/actions
Nx   = 32;                                  % number of cells per dimension
Ns   = Nx*Nx;                               % number of states
dt   = 2;
 
% fisrt and final states
%--------------------------------------------------------------------------
try, T;  catch, T = 16;    end
try, X;  catch, X = [0 1]; end
try, V;  catch, V = [0 0]; end
 
X(1)  = X(1);                               % first position
V(1)  = V(1);                               % first velocity
X(T)  = X(2);                               % final position
V(T)  = V(2);                               % final velocity
 
 
% set-up transition matrix for mountain car problem
%--------------------------------------------------------------------------
x     = linspace(-2,2,Nx);
v     = linspace(-3,3,Nx);
 
K     = toeplitz(sparse(1,[1 2],[1 1/2],1,Nx));
K     = K + K';
K     = K*diag(1./sum(K,1));
 
for k = 1:Nu
    P{k}  = sparse(Ns,Ns);                 % pushforward
    B{k}  = sparse(Ns,Ns);                 % pullback
    for i = 1:length(x)
        for j = 1:length(v)
            
            % change in state
            %--------------------------------------------------------------
            ds = spm_fx_mountaincar([x(i);v(j)],0,u(k),[],16);
            
            
            % transition probabilities - position
            %--------------------------------------------------------------
            dx = x - (x(i) + ds(1));
            ii = find(dx > 0,1);
            if ii == 1,
                px = sparse(1,1,1,Nx,1);
            elseif isempty(ii)
                px = sparse(Nx,1,1,Nx,1);
            else
                ii = [ii - 1,ii];
                px = pinv([x(ii); 1 1])*[(x(i) + ds(1)); 1];
                px = sparse(ii,1,px,Nx,1);
            end
            
            % transition probabilities - velocity
            %--------------------------------------------------------------
            dv = v - (v(j) + ds(2));
            ii = find(dv > 0,1);
            if ii == 1,
                pv = sparse(1,1,1,Nx,1);
            elseif isempty(ii)
                pv = sparse(Nx,1,1,Nx,1);
            else
                ii = [ii - 1,ii];
                pv = pinv([v(ii); 1 1])*[(v(j) + ds(2)); 1];
                pv = sparse(ii,1,pv,Nx,1);
            end
            
            % place in P
            %--------------------------------------------------------------
            px     = K*px;
            pv     = K*pv;
            p      = px*pv';
            P{k}(:,(j - 1)*Nx + i) = p(:);
            
        end
    end
end
 
% Parameters of MDP model
%--------------------------------------------------------------------------
A     = speye(Ns,Ns);
A     = A*diag(1./sum(A,1));
for i = 1:Nu
    B{i}   = P{i}';
    P{i}   = P{i}*diag(1./sum(P{i},1));
    B{i}   = B{i}*diag(1./sum(B{i},1));
end
 
% make into discrete transition probabilities
%--------------------------------------------------------------------------
for i = 1:Nu
    P{i} = P{i}^dt;
    B{i} = B{i}^dt;
end
 
% graphics (transition probabilities)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf
 
[d i] = min(abs(x - X(2)));
[d j] = min(abs(v - V(2)));
j     = (j - 1)*Nx + i;
s     = sparse(j,1,1,Ns,1);
 
k     = [0 4 8];
for i = 1:Nu
    for j = 1:3
        subplot(Nu,3,(i - 1)*3 + j)
        imagesc(x,v,reshape(B{i}^k(j)*s,Nx,Nx)'), axis xy square
        xlabel('position','FontSize',12)
        ylabel('velocity','FontSize',12)
    end
end
drawnow;
 
% initial and final states
%==========================================================================
spm_figure('GetWin','Figure 2');clf
 
 
% first state (s)
%--------------------------------------------------------------------------
[d i] = min(abs(x - X(1)));
[d j] = min(abs(v - V(1)));
X(1)  = x(i);
V(1)  = v(j);
s     = sparse(i,j,1,Nx,Nx);
s     = s(:);
 
% final state (c)
%--------------------------------------------------------------------------
[d i] = min(abs(x - X(T)));
[d j] = min(abs(v - V(T)));
X(T)  = x(i);
V(T)  = v(j);
c     = sparse(i,j,1,Nx,Nx);
c     = c(:);
 
% (uniform) cost over control (d)
%--------------------------------------------------------------------------
d     = ones(Nu,1);
c     = c/sum(c);
d     = d/sum(d);
 
% ensure normalization and get log-transforms
%--------------------------------------------------------------------------
lnA   = log(A + exp(-TOL));
for i = 1:Nu
    lnB{i} = log(B{i} + exp(-TOL));
end
lnc   = log(c);
lnd   = log(d);
 
 
% effector action (e) and sufficient statistics of proposal density
%--------------------------------------------------------------------------
e      = zeros(Nu,T);                       % action
a      = ones(Ns,T);                        % state probability
b      = ones(Nu,T);                        % control probability
a(:,1) = s;                                 % first state
a(:,T) = c;                                 % final state
a      = a*diag(1./sum(a,1));
b      = b*diag(1./sum(b,1));
 
% set up state space graphically
%--------------------------------------------------------------------------
subplot(2,1,1)
for i = 1:length(x);
    plot(x(i),v,'k','color',[1 1 1]/2), hold on
    axis([-2 2 -2 2])
end
 
% solve
%==========================================================================
[i j] = find(reshape(s(:,1),Nx,Nx));
Sx(1) = X(1);
Sv(1) = V(1);
for k = 1:(T - 1)
    
    % forward and backward passes at this time point
    %----------------------------------------------------------------------
    for i = 1:4
        for t = [(T - 1):-1:max(k,2) max(k,2):(T - 1)]
            
            % get data likelihood if available at this time
            %--------------------------------------------------------------
            try
                at = lnA'*s(:,t);
            catch
                at = sparse(Ns,1);
            end
            
            % and accumulate empirical priors
            %--------------------------------------------------------------
            for j = 1:Nu
                at = at + b(j,t + 1) * lnB{j} *a(:,t + 1);
                at = at + b(j,t    ) * lnB{j}'*a(:,t - 1);
            end
            for j = 1:Nu
                bt(j,1) = lnd(j);
                bt(j,1) = bt(j,1) + a(:,t - 1)'* lnB{j} *a(:,t    );
            end
            
            % update sufficient statistics of hidden states and control
            %--------------------------------------------------------------
            at     = exp(at);
            bt     = exp(bt);
            a(:,t) = at/sum(at);
            b(:,t) = bt/sum(bt);
            
            % conditional expectations
            %--------------------------------------------------------------
            Q      = reshape(a(:,t),Nx,Nx);
            Qx     = sum(Q,2);
            Qv     = sum(Q,1);
            X(t)   = x*Qx /sum(Qx);
            V(t)   = v*Qv'/sum(Qv);
            
        end
        
        % plot
        %------------------------------------------------------------------
        subplot(2,1,1)
        plot(X,V,':','color',[1 1 1]*.8), hold on
        axis([-2 2 -2 2])
        
        subplot(2,2,3)
        imagesc(b); axis square
        title('Anticipated control','FontSize',16)
        xlabel('time','FontSize',12)
        drawnow
    end
    
    
    % action-dependent selection of next sensory state
    %----------------------------------------------------------------------
    for i = 1:Nu
        ss(i)  = P{i}(:,find(s(:,k)))'*lnA*a(:,k + 1);
    end
    [m i]      = max(ss);                    % optimum action
    [m j]      = max(P{i}*s(:,k));           % next state index
    s(:,k + 1) = sparse(j,1,1,Ns,1);         % next state
    e(i,k + 1) = 1;
    
    % plot (sensed) position
    %----------------------------------------------------------------------
    [i j]     = find(reshape(s(:,k + 1),Nx,Nx));
    Sx(k + 1) = x(i);
    Sv(k + 1) = v(j);
    
    subplot(2,1,1)
    plot(Sx,Sv,'ok','LineWidth',4)
    plot(Sx,Sv,'k','LineWidth',2)
    title('Anticipated and actual trajectory','FontSize',16)
    xlabel('position','FontSize',12)
    ylabel('velocity','FontSize',12)
    axis([-2 2 -2 2])
    
    subplot(2,2,4)
    imagesc(e); axis square
    title('Actual action','FontSize',16)
    xlabel('time','FontSize',12)
    drawnow
    
end

