function spm_MDP_urn
% Demo for active inference with the urn task
%__________________________________________________________________________
%
% This demonstration uses the Urn or Beads Task to illustrate how choice
% behaviour can be simulated using active inference - in the context of
% Markov decision processes. In the urn task, a succession of draws
% from one of two urns are made and the agent has to decide whether the
% balls are being drawn from an urn with predominantly red or green balls.
% We model this in terms of a state-space with four dimensions: number of
% balls drawn (n),  number of green balls drawn (k), choice (undecided, red 
% or green)and the true (hidden) state of the urn (red or green). With
% this relatively simple state-space, the utility of any hidden state is
% simply quantified by the log-odds ratio of making a correct
% decision. From binomial theory this is (2k - n)*log(p/(1 - p)), where p
% is the proportion of red or green balls. Having defined the utility
% function of states, we can then use the MDP formulation of active
% inference using variational Bayes to simulate choice behaviour. 
%
% This routine first provides an illustration of a game in which a decision
% is delayed until the last draw to look at inferences during successive
% draws - with a special focus on precision. The illustration here shows
% a decrease in precision when an unexpected (green ball) is drawn during a
% sequence of red balls.
%
% We then characterise changes in choice probability (and latency to the
% decision) in terms of its dependency on threshold criteria (on the odds
% ratio) and hyperpriors about precision (alpha or the scale parameter of a 
% standard gamma distribution). The routine concludes with an illustration 
% of how to estimate model parameters using the likelihood of observed
% (simulated) choices.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_urn.m 6592 2015-11-06 16:20:48Z guillaume $
 
% set up and preliminaries
%==========================================================================
rng('default')

T     = 8;                         % number of trials in each game
Pa    = .85;                       % proportion of red or green goals
 
% hidden (and initial) states (a red urn with no drawers)
%--------------------------------------------------------------------------
D          = zeros(T,T,3,2);       % #balls - 1 x #green - 1 x #u x #urns
D(1,1,1,1) = 1;
S          = spm_vec(D);
D(1,1,1,2) = 1;
D          = spm_vec(D);

 
% likelihood - everything is seen  apart from the hidden states of the urn
%--------------------------------------------------------------------------
A     = kron([1 1],eye(3*T*T));

 
% transition probabilities (B{1} - wait; B{2} - red; B{3} - green)
%--------------------------------------------------------------------------
for i = 1:6
    Bn{i,i} = kron(spm_speye(T,T),spm_speye(T,T));
end
Br      = spm_speye(T,T,0)*Pa + spm_speye(T,T,-1)*(1 - Pa);
Bg      = spm_speye(T,T,0)*(1 - Pa) + spm_speye(T,T,-1)*Pa;
Bn{1,1} = kron(Br,spm_speye(T,T,-1));
Bn{4,4} = kron(Bg,spm_speye(T,T,-1));
 
B{1}    = spm_cat(Bn);
B{2}    = kron(eye(2),kron([0 0 0;1 1 0;0 0 1],eye(T*T)));
B{3}    = kron(eye(2),kron([0 0 0;0 1 0;1 0 1],eye(T*T)));


% priors over final state (exp(utility))
%--------------------------------------------------------------------------
L     = zeros(T,T,3,2);
W     = zeros(T,T,3,2);
rho   = log(Pa/(1 - Pa));
for n = 0:(T - 1)
    W(n + 1,:,:,:) = - n/8;
    for k = 0:n
        L(n + 1,k + 1,2,1) = (n - 2*k)*rho;
        L(n + 1,k + 1,3,2) = (2*k - n)*rho;
    end
end
C     = spm_softmax(spm_vec((L > 3) + W)*4);

 
% allowable policies (one decision before the game ends)
%--------------------------------------------------------------------------
V     = [ones(T,T) + eye(T,T), ones(T,T) + 2*eye(T,T)];
 
 
% MDP Structure
%==========================================================================
MDP.T = T;                          % process depth (the horizon)
MDP.S = S;                          % initial state
MDP.A = A;                          % likelihood
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % prior over initial states
MDP.V = V;                          % allowable policies
 
% Solve - an example game (with high offer at t = 10)
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
 
% create a sequence of draws (outcomes) - with and oddball on the 4th
%--------------------------------------------------------------------------
k     = cumsum([0 0 0 1 0 0 0 0]);
a     = [1 1 1 1 1 1 2];
o     = zeros(1,length(a) - 1);
for n = 1:length(o)
    o(n) = sub2ind([T,T,3],n,k(n) + 1,1);
end
 
% Active inference (precluding a decision until the last trial)
%--------------------------------------------------------------------------
MDP.o    = o; 
MDP.a    = a;
MDP.plot = gcf;
MDP.N    = 8;
MDP      = spm_MDP_game(MDP);
 
% plot convergence and precision
%--------------------------------------------------------------------------
subplot(4,2,7)
plot(MDP.d)
xlabel('Latency (updates)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Expected precision','FontSize',16)
spm_axis tight
 
% deconvolve to simulate dopamine responses
%--------------------------------------------------------------------------
nd  = length(MDP.d);
K   = tril(toeplitz(exp(-((1:nd) - 1)'/8)));
 
subplot(4,2,8)
plot(pinv(K)*MDP.d(:)), hold on
xlabel('Latency (updates)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Simulated dopamine responses','FontSize',16)
spm_axis tight


% Illustrate dependency on parameters
%==========================================================================
spm_figure('GetWin','Figure 2'); clf

% create a sequence of draws (outcomes) - no oddballs
%--------------------------------------------------------------------------
k     = cumsum([0 0 0 0 0 0 0 0]);
o     = zeros(1,length(k));
for n = 1:length(k)
    o(n) = sub2ind([T,T,3],n,k(n) + 1,1);
end

 
% probability distribution over time: P(1,:) is no action
%--------------------------------------------------------------------------
PrT      = @(P) [1 cumprod(P(1,1:end - 1))].*(1 - P(1,:));
MDP.plot = 0;                        % plot convergence
MDP.N    = 4;                        % number of variational iterations
MDP.a    = ones(1,T);                % and action
MDP.o    = o;                        % and outcomes

% beliefs about final state - decision threshold (log likelihood)
%--------------------------------------------------------------------------
DP    = MDP;
PF    = [];
DF    = [];
p     = linspace(0,8,8);
for i = 1:length(p)
    DP.C    = spm_softmax(spm_vec((L > p(i)) + W)*4);
    DP      = spm_MDP_game(DP);
    PF(i,:) = 1 - DP.P(1,:);
    DF(i,:) = PrT(DP.P);
end
 
% probability of accepting
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(PF')
xlabel('Latency (trials)','FontSize',12)
ylabel('Probability of deciding','FontSize',12)
title('Increasing decision threshold','FontSize',16)
axis square xy
 
% distribution of acceptance latencies
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(DF')
xlabel('Latency (trials)','FontSize',12)
ylabel('Latency of decision','FontSize',12)
title('Latency of decision','FontSize',16)
axis square xy
 
 
% Hyperpriors - prior precision (alpha)
%--------------------------------------------------------------------------
DP    = MDP;
PF    = [];
DF    = [];
p     = linspace(2,16,8);
for i = 1:length(p)
    DP.alpha  = p(i);
    DP      = spm_MDP_game(DP);
    PF(i,:) = 1 - DP.P(1,:);
    DF(i,:) = PrT(DP.P);
end
 
% probability of accepting
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(PF')
xlabel('Latency (trials)','FontSize',12)
ylabel('Probability of deciding','FontSize',12)
title('Increasing prior precision','FontSize',16)
axis square xy
 
% distribution of acceptance latencies
%--------------------------------------------------------------------------
subplot(2,2,4)
plot(DF')
xlabel('Latency (trials)','FontSize',12)
ylabel('Latency of decision','FontSize',12)
title('Latency of decision','FontSize',16)
axis square xy
 
 
% Simulate multiple trials and try to infer prior precision (alpha)
%==========================================================================
spm_figure('GetWin','Figure 3'); clf
 
 
% Simulate multiple trials and record likelihood
%--------------------------------------------------------------------------
MDP.plot = 0;
MDP.N    = 4;
 
p     = linspace(2,16,8);
for t = 1:8
    
    % run game
    %----------------------------------------------------------------------
    MDP.s = [];
    MDP.o = [];
    MDP.a = [];
    MDP   = spm_MDP_game(MDP);
    
    % place outcomes in DP
    %----------------------------------------------------------------------
    DP    = MDP;
    [a j] = find(MDP.U);
    [o j] = find(MDP.O);
    DP.o  = o;
    DP.a  = a;
    y(t)  = find(a > 1);
        
    % get log-likelihood for different parameter values
    %----------------------------------------------------------------------
    for i = 1:length(p);
        DP.alpha = p(i);
        DP       = spm_MDP_game(DP);
        LL(i,t)  = sum(log(DP.P(find(DP.U))));
    end
    
end
 
% approximate the MAP with the ML and use the Laplace assumption
%--------------------------------------------------------------------------
Lp    = sum(LL,2);
dp    = mean(diff(p,1));
dLdpp = diff(Lp,2)/(dp^2);
[l i] = max(Lp(2:end - 1) + (dLdpp < 0)*1024);
Cp    = inv(-dLdpp(i));
Ep    = p(i + 1);
 

% plot responses
%--------------------------------------------------------------------------
subplot(2,2,1)
hist(y,1:T)
xlabel('Latency','FontSize',12)
ylabel('Probability','FontSize',12)
title('Distribution of responses','FontSize',16)
axis square
  
% plot log likelihood over trials
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(p,LL - min(LL(:)))
xlabel('Parameter','FontSize',12)
ylabel('Probability','FontSize',12)
title('Log-likelihood over games','FontSize',16)
axis square
 
% plot likelihood
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(p,Lp - min(Lp))
xlabel('Parameter','FontSize',12)
ylabel('Log-likelihood','FontSize',12)
title('Log-likelihood','FontSize',16)
axis square
    
% plot posterior
%--------------------------------------------------------------------------
subplot(2,2,4)
tp  = 8;
pp  = linspace(p(1),p(end),64);
plot(pp,spm_Npdf(pp,Ep,Cp)), hold on
plot([tp tp],[0 1/8],':'),   hold off
xlabel('Parameter','FontSize',12)
ylabel('Probability','FontSize',12)
title('Posterior probability','FontSize',16)
axis square
 
 
return
