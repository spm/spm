function spm_MDP_offer
% Demo for active inference with limited offer game
%__________________________________________________________________________
%
% This routine uses a Markov decisions process formulation of the mountain
% car problem to illustrate prospective free energy minimization under a
% variational Bayesian learning scheme. The key notion here is that the
% agent represents future states and action (in a pullback sense), where it
% has strong prior beliefs about future states. The intervening states and
% actions are optimized with respect to current sensory data to provide
% predictions about the next sensory state, which action fulfils. The
% result is a planned trajectory through state space that realizes prior
% beliefs in a prospective sense.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP.m 4804 2012-07-26 13:14:18Z karl $
 
% set up and preliminaries
%==========================================================================
 
% initial and final states
%==========================================================================
T     = 16;                         % number of offers
Pa    = 1/2;                        % probability of a high offer
Pb    = 1/2;                          % probability of withdrawn offer
   
% transition probabilities (B{1} - decline; B{2} - accept)
%--------------------------------------------------------------------------
for i = 1:T
    
    a       = 1 - (1 - Pa)^(1/T);
    b       = (i*Pb/T)*(1 - a);
    B{i,1}  = [(1 - a - b) 0 0 0 0;
                a          0 0 0 0;
                b          1 1 0 0;
                0          0 0 1 0;
                0          0 0 0 1];
    
    B{i,2}  = [ 0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 1 1;
                1 0 0 0 0;
                0 1 0 0 0];
end
      
 
% initial state
%--------------------------------------------------------------------------
S     = [1 0 0 0 0]';
 
% priors over final state
%--------------------------------------------------------------------------
C     = [0 0 0 1 2]';
 
% prior over control
%--------------------------------------------------------------------------
D     = [1 1]';

% prior over states
%--------------------------------------------------------------------------
E     = [1 1 0 1 1]';
 
 
% solve - an example trial
%==========================================================================
MDP.T = T;                          % process depth (the horizon)
MDP.S = S;                          % initial state
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % control probabilities (priors)
MDP.E = E;                          % state probabilities (priors)
MDP.W = 4;                          % log-precision


MDP.plot = 1;                          % plot convergence



spm_MDP(MDP);

% Generate process (with continuous low offers)
%==========================================================================
 
% True transition probabilities (retaining low offer for all trials)
%--------------------------------------------------------------------------
G{1,1}  = [1 1 0 0 0;
           0 0 0 0 0;
           0 0 1 0 0;
           0 0 0 1 0;
           0 0 0 0 1];
 
G{1,2}   = G{1,1};
MDP.G    = G;
MDP.plot = 0;
 
% probability distribution over time to act
%--------------------------------------------------------------------------
PrY   = @(P)[1 cumprod(P(1,1:end - 1))].*P(2,:);
 
 
% illustrate dependency on latency w.r.t. parameters
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
 
% high offer
%--------------------------------------------------------------------------
p     = linspace(1/16,1 - 1/16,4);
MDPP  = MDP; clear Py str
for i = 1:length(p)
    
    % transition probabilities
    %----------------------------------------------------------------------
    for j = 1:T
        a                  = 1 - (1 - p(i))^(1/T);
        b                  = (j*Pb/T)*(1 - a);
        MDPP.B{j,1}(1:3,1) = [(1 - a - b); a; b];        
    end
    
    [Q,R,S,E,P] = spm_MDP(MDPP);
    Py(i,:)     = PrY(P);
    str{i}      = num2str(p(i));
end
 
subplot(2,2,1)
plot(1:length(P),Py,'.','MarkerSize',16),legend(str), hold on
plot(1:length(P),Py,':'), hold off
xlabel('latency (offers)','FontSize',12)
ylabel('acceptance probability ','FontSize',12)
title('high offer probability ','FontSize',16)
axis square
legend(str)
 
% offer withdrawal
%--------------------------------------------------------------------------
p     = linspace(1/16,1,4);
MDPP  = MDP; clear Py str
for i = 1:length(p)
    
    % transition probabilities
    %----------------------------------------------------------------------
    for j = 1:T
        a                  = 1 - (1 - Pa)^(1/T);
        b                  = (j*p(i)/T)*(1 - a);
        MDPP.B{j,1}(1:3,1) = [(1 - a - b); a; b];        
    end
    
    [Q,R,S,E,P] = spm_MDP(MDPP);
    Py(i,:)     = PrY(P);
    str{i}      = num2str(p(i));
end
 
subplot(2,2,2)
plot(1:length(P),Py,'.','MarkerSize',16),legend(str), hold on
plot(1:length(P),Py,':'), hold off
xlabel('latency (offers)','FontSize',12)
ylabel('acceptance probability','FontSize',12)
title('offer withdrawal','FontSize',16)
axis square
 
 
 
% cost bias C(5)
%--------------------------------------------------------------------------
p     = linspace(1/2,2,4); 
MDPP  = MDP; clear Py str
for i = 1:length(p)
    MDPP.C(5)   = p(i);
    [Q,R,S,E,P] = spm_MDP(MDPP);
    Py(i,:)     = PrY(P);
    str{i}      = num2str(p(i));
end
 
subplot(2,2,3)
plot(1:length(P),Py,'.','MarkerSize',16),legend(str), hold on
plot(1:length(P),Py,':'), hold off
xlabel('latency (offers)','FontSize',12)
ylabel('acceptance probability','FontSize',12)
title('Behaviour and terminal cost','FontSize',16)
axis square
legend(str)
 
% precision
%--------------------------------------------------------------------------
p     = linspace(2,5,4);
MDPP  = MDP; clear Py str
for i = 1:length(p)
    MDPP.W      = p(i);
    [Q,R,S,E,P] = spm_MDP(MDPP);
    Py(i,:)     = PrY(P);
    str{i}      = num2str(p(i));
end
 
subplot(2,2,4)
plot(1:length(P),Py,'.','MarkerSize',16),legend(str), hold on
plot(1:length(P),Py,':'), hold off
xlabel('latency (offers)','FontSize',12)
ylabel('acceptance probability','FontSize',12)
title('Behaviour and precicion','FontSize',16)
axis square
legend(str)
 
    

 
% Changes in uncertainty (Entropy) over successive choices
%==========================================================================
spm_figure('GetWin','Figure 2'); clf
 
[Q,R,S,E,P] = spm_MDP(MDP);
 
 
% uncertainty about current action
%--------------------------------------------------------------------------
H   = sum(-P.*log(P));
 
subplot(2,2,1)
plot(1:length(H),H,'.','MarkerSize',16), hold on
plot(1:length(H),H,':'), hold off
xlabel('latency','FontSize',12)
ylabel('Entropy','FontSize',12)
title('current action','FontSize',16)
axis square
 
 
% uncertainty about future action
%--------------------------------------------------------------------------
for i = 1:(T - 1)
    j    = (i + 1):(T - 1);
    p    = squeeze(R(:,i,j));
    H(i) = mean(sum(-p.*log(p)));  
end
 
subplot(2,2,2)
plot(1:length(H),H,'.','MarkerSize',16), hold on
plot(1:length(H),H,':'), hold off
xlabel('latency','FontSize',12)
ylabel('Entropy','FontSize',12)
title('future action','FontSize',16)
axis square
 
% uncertainty about future states
%--------------------------------------------------------------------------
for i = 1:(T - 1)
    j    = (i + 1):(T - 1);
    p    = squeeze(Q(:,i,j));
    H(i) = mean(sum(-p.*log(p)));  
end
 
subplot(2,2,3)
plot(1:length(H),H,'.','MarkerSize',16), hold on
plot(1:length(H),H,':'), hold off
xlabel('latency','FontSize',12)
ylabel('Entropy','FontSize',12)
title('future states','FontSize',16)
axis square
 
% uncertainty about future states
%--------------------------------------------------------------------------
for i = 1:(T - 1)
    for j = 1:(T - 1)
        p      = squeeze(Q(:,i,j));
        H(i,j) = sum(-p.*log(p)); 
    end
end
 
subplot(2,2,4)
imagesc(H)
xlabel('future state','FontSize',12)
ylabel('current trial','FontSize',12)
title('Over the future','FontSize',16)
axis square
 
 
return



% simulate multiple trials and record when an offer was accepted
%==========================================================================
spm_figure('GetWin','Figure 3'); clf

% trials
%--------------------------------------------------------------------------
for i = 1:32
    [Q,R,S,E,P] = spm_MDP(MDP);
    try
        Y(i)    = find(E(2,:),1);
    end
    fprintf('trial %0.00f\n',i);
end
 
% probability distribution over time to act
%--------------------------------------------------------------------------
Py    = PrY(P);
 
 
% plot
%--------------------------------------------------------------------------
subplot(2,2,1)
hist(Y,1:T);
xlabel('choice latency','FontSize',12)
ylabel('sample frequnecy','FontSize',12)
title('sample distribution of latencies','FontSize',16)
axis square
 
subplot(2,2,2)
bar(Py)
xlabel('choice latency','FontSize',12)
ylabel('probability','FontSize',12)
title('predicted probability','FontSize',16)
axis square
 
 
% Infer prior beliefs from observed responses (meta-modelling)
%==========================================================================
p     = linspace(1/16,1 - 1/16,32);
MDPP  = MDP;
for i = 1:length(p);
    
    % transition probabilities
    %----------------------------------------------------------------------
    for j = 1:T
        a                  = 1 - (1 - p(i))^(1/T);
        b                  = (j*Pb/T)*(1 - a);
        MDPP.B{j,1}(1:3,1) = [(1 - a - b); a; b];        
    end
    
    % get likelihood for this parameter
    %----------------------------------------------------------------------
    [Q,R,S,E,P] = spm_MDP(MDPP);
    Py          = PrY(P);
    L(i)        = sum(log(Py(Y)));
    
end
 
% approximate the MAP with the ML and use the Laplace assumption
%--------------------------------------------------------------------------
[l i] = max(L);
dp    = p(2) - p(1);
dLdpp = (L(i + 1) - L(i) - L(i) + L(i - 1))/(dp^2);
Cp    = inv(-dLdpp);
Ep    = p(i);
 
 
% plot likelihood
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(p,L)
xlabel('latency','FontSize',12)
ylabel('probabaility','FontSize',12)
title('log-likelihood','FontSize',16)
axis square
    
% plot posterior
%--------------------------------------------------------------------------
subplot(2,2,4)
plot(p,spm_Npdf(p,Ep,Cp)), hold on
plot([Pa Pa],[0 8],':'),   hold off
xlabel('latency','FontSize',12)
ylabel('probability','FontSize',12)
title('posterior probability','FontSize',16)
axis square
 
