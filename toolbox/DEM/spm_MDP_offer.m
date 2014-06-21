function spm_MDP_offer
% Demo for active inference with limited offer game
%__________________________________________________________________________
%
% This demonstration routine uses variational Bayes to minimise the free
% energy to model decision-making. The particular focus here is on
% decisions that are time-sensitive, requiring an explicit representation
% of future states. The example considered here represents a limited offer
% game, where a low offer can be converted to a high offer, which may or
% may not occur. Furthermore, offers may be withdrawn. The objective is
% to understand model choices about accepting or declining the current
% offer in terms of active inference, under prior beliefs about future
% states. The model is specified in a fairly general way in terms of
% probability transition matrices and beliefs about future states. The
% particular inversion scheme used here is spm_MDP_game, which uses a
% mean-field approximation between hidden control and hidden states. It is
% assumed that the agent believes that it will select a particular action
% (accept or decline) at a particular time.
%
% We run an exemplar game, examine the distribution of time to acceptance
% as a function of different beliefs (encoded by parameters of the
% underlying Markov process) and demonstrate how the model can be used to
% produce trial-specific changes in uncertainty – or how one can use
% behaviour to identify the parameters used by a subject.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_offer.m 6062 2014-06-21 11:00:15Z karl $
 
% set up and preliminaries
%==========================================================================
T     = 16;                           % number of offers
Pa    = 1/2;                          % probability of a high offer
Pb    = 1/16;                         % probability of withdrawn offer
Plos  = @(t,Pb)(1 - (1 - Pb).^t);
Pwin  = @(T,Pa)(1 - (1 - Pa)^(1/T));
 
 
% transition probabilities (B{1} - decline; B{2} - accept)
%--------------------------------------------------------------------------
for t = 1:T
    
    a       = Pwin(T,Pa);
    b       = Plos(t,Pb);
    B{t,1}  = [(1 - a + a*b - b) 0 0 0 0;
                a*(1 - b)        0 0 0 0;
                b                1 1 0 0;
                0                0 0 1 0;
                0                0 0 0 1];
    
    B{t,2}  = [ 0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 1 1;
                1 0 0 0 0;
                0 1 0 0 0];
end
      
 
% initial state
%--------------------------------------------------------------------------
S     = [1 0 0 0 0]';
 
% priors over final state (exp(utility))
%--------------------------------------------------------------------------
C     = spm_softmax([1 1 1 2 4]');

% allowable policies (one shift at different times)
%--------------------------------------------------------------------------
V     = eye(T,T) + 1;
 
 
% MDP Structure
%==========================================================================
MDP.T = T;                          % process depth (the horizon)
MDP.S = S;                          % initial state
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.V = V;                          % allowable policies
 
% Solve - an example game (with high offer at t = 10)
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
 
MDP.s    = [1 1 1 1 1 1 1 1 1 2];
MDP.a    = [1 1 1 1 1 1 1 1 1];
MDP.plot = gcf;
MDP.N    = 8;
MDP      = spm_MDP_game(MDP);

% plot convergence and precision
%--------------------------------------------------------------------------
subplot(4,2,7)
plot(MDP.W)
xlabel('Latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Expected precision','FontSize',16)
spm_axis tight

% deconvolve to simulate dopamine responses
%--------------------------------------------------------------------------
subplot(4,2,8)
plot(MDP.da), hold on
xlabel('Latency (iterations)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Simulated dopamine responses','FontSize',16)
spm_axis tight


% Solve - an example game (with low offer at t = 5)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
MDP.s    = [1 1 1 1 3];
MDP.a    = ones(1,T);
MDP.plot = gcf;
MDP      = spm_MDP_game(MDP);
 
% plot convergence and precision
%--------------------------------------------------------------------------
subplot(4,2,7)
plot(MDP.W)
xlabel('Latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Expected precision','FontSize',16)
spm_axis tight
 
subplot(4,2,8)
plot(MDP.da), hold on
xlabel('Latency (iiterations)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Simulated dopamine responses','FontSize',16)
spm_axis tight


% Illustrate dependency parameters
%==========================================================================
spm_figure('GetWin','Figure 3'); clf
 
% probability distribution over time: P(1,:) is no action
%--------------------------------------------------------------------------
PrT      = @(P) [1 cumprod(P(1,1:end - 1))].*(1 - P(1,:));
MDP.plot = 0;                        % plot convergence
MDP.N    = 4;                        % number of variational iterations
MDP.s    = ones(1,T);                % suppress withdrawal of low offer
MDP.a    = ones(1,T);                % and action
 
% beliefs about final state
%--------------------------------------------------------------------------
DP    = MDP;
p     = linspace(0,8,16);
for i = 1:length(p)
    DP.C    = spm_softmax([1 1 1 p(i) 4]');
    DP      = spm_MDP_game(DP);
    BF      = DP.P;
    DP      = spm_MDP_game(DP,'Expected Utility');
    BE      = DP.P;
    PF(i,:) = BF(2,:);
    PE(i,:) = BE(2,:);
    DF(i,:) = PrT(BF);
    DE(i,:) = PrT(BE);
end
 
% probability of accepting
%--------------------------------------------------------------------------
subplot(2,2,1)
imagesc(1:(T - 1),p,1 - PF)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Conditional divergence','FontSize',16)
axis square xy
 
% compare with expected utility
%--------------------------------------------------------------------------
subplot(2,2,2)
imagesc(PE)
imagesc(1:(T - 1),p,1 - PE)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Expected utility','FontSize',16)
axis square xy
 
% distribution of acceptance latencies
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(1:T,p,1 - DF)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Latency of accepting','FontSize',16)
axis square xy
 
% compare with expected utility
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(1:T,p,1 - DE)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Latency of accepting','FontSize',16)
axis square xy
 
 
% Changes in uncertainty (Entropy) over successive choices
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
 
% uncertainty about current action
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = ones(1,T);
MDP.o = ones(1,T);
 
MDP.plot = 0;
MDP.N    = 4;
 
MDP.C    = spm_softmax([1 1 1 3 8]');
MDP      = spm_MDP_game(MDP);
P        = MDP.P;
W        = MDP.W;
H        = sum(-P.*log(P),1);
 
% action entropy
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(1:length(P),P,'-.'), hold on
plot(1:length(H),H,'.r','MarkerSize',16), hold on
plot(1:length(H),H,':r'), hold off
xlabel('Offer','FontSize',12)
ylabel('Probability and entropy (nats)','FontSize',12)
title('Uncertainty about action','FontSize',16)
spm_axis tight
axis square
 
% precision
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(1:length(W),W,'.','MarkerSize',16), hold on
plot(1:length(W),W,':'), hold off
xlabel('Offer','FontSize',12)
ylabel('Precision','FontSize',12)
title('Precision','FontSize',16)
spm_axis tight
axis square
 
% expected utility and entropy (components of expected divergence)
%--------------------------------------------------------------------------
beta = 1;
D    = 1./W - beta;
MDP  = spm_MDP_game(MDP,'Expected Utility');
W    = MDP.W;
V    = 1./W - beta;
 
subplot(2,2,3)
plot(1:length(D),V - D,'.r','MarkerSize',16), hold on
plot(1:length(V),-V,   '.b','MarkerSize',16), hold off
xlabel('Offer','FontSize',12)
ylabel('Expected values','FontSize',12)
title('Expected utility','FontSize',16)
legend({'Entropy','Expected utility'})
spm_axis tight
axis square

 
% Effective utility (the effect of optimising precision or sensitivity)
%==========================================================================
spm_figure('GetWin','Figure 5'); clf
 
% how does precision depend on beliefs about high offer?
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = ones(1,T);
MDP.o = ones(1,T);
MDP.N = 4;
MDP.K = 1;
 
DP    = MDP;
p     = linspace(2,8,16);
for i = 1:length(p)
 
    % high offer utility (with an EU agent)
    %----------------------------------------------------------------------
    DP.C    = spm_softmax([1 1 1 2 p(i)]');
    UP(i,:) = log(DP.C);
    DP      = spm_MDP_game(DP);
    DW(i,:) = DP.W;
    
end
 
subplot(2,2,1)
imagesc(1:T,p,DW)
title('Precision','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('Utility of high offer','FontSize',12)
axis square xy
 
subplot(2,2,2)
plot(1:T,DW','r')
title('Precision','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('Precision','FontSize',12)
spm_axis tight
axis square
 
% subjective utility
%--------------------------------------------------------------------------
subplot(2,2,3)
UP    = UP(:,[4 5]);
for i = 1:T
    SU = diag(DW(:,i))*UP ;
    plot(UP,SU), hold on
end
hold off, title('Subjective utility','FontSize',16)
xlabel('utility of low and high offers','FontSize',12)
ylabel('Subjective (behavioural) utility','FontSize',12)
spm_axis tight
axis square

% the efect of decreasing (fixed) precision
%--------------------------------------------------------------------------
DP    = MDP;
DP.C  = spm_softmax([1 1 1 2 3]');
p     = linspace(0,1,16);
for i = 1:length(p)
    DP.w    = zeros(1,T) + p(i);
    DP      = spm_MDP_game(DP,'Expected Utility');
    DT(i,:) = PrT(DP.P);
end

subplot(2,2,4)
imagesc(1:T,p,1 - DT)
xlabel('Latency (offers)','FontSize',12)
ylabel('Fixed precision','FontSize',12)
title('Latency of accepting','FontSize',16)
axis square xy

 
% Simulate multiple trials and record when an offer was accepted
%==========================================================================
spm_figure('GetWin','Figure 6'); clf
 
% trials with no higher offer
%--------------------------------------------------------------------------
MDP.C = spm_softmax([1 1 1 2 4]');
MDP.s = ones(1,T);
MDP.a = [];
MDP.o = [];
 
MDP.plot = 0;
MDP.N    = 4;
 
% trials
%--------------------------------------------------------------------------
for i = 1:256
    MDP = spm_MDP_game(MDP);
    try
        Y(i)  = find(MDP.U(2,:),1);
    catch
        Y(i)  = T;
    end
    fprintf('trial %0.00f\n',i);
end
 
% probability distribution over time to act
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = ones(1,T);
MDP   = spm_MDP_game(MDP);
Py    = PrT(MDP.P);
 
% plot
%--------------------------------------------------------------------------
subplot(2,2,1)
hist(Y,1:T);
xlabel('Acceptance latency','FontSize',12)
ylabel('Sample frequency','FontSize',12)
title('Sample distribution of latencies','FontSize',16)
axis square
 
subplot(2,2,2)
bar(Py)
xlabel('Acceptance latency','FontSize',12)
ylabel('Probability','FontSize',12)
title('Predicted probability','FontSize',16)
axis square
 
 
% Infer utility from observed responses (meta-modelling)
%==========================================================================
p     = linspace(1,6,32);
DP    = MDP;
for i = 1:length(p);
    
    % transition probabilities
    %----------------------------------------------------------------------
    DP.C  = spm_softmax([1 1 1 2 p(i)]');
        
    % get likelihood for this parameter
    %----------------------------------------------------------------------
    DP    = spm_MDP_game(DP);
    Py    = PrT(DP.P);
    L(i)  = sum(log(Py(Y) + eps));
    
end
 
% approximate the MAP with the ML and use the Laplace assumption
%--------------------------------------------------------------------------
dp    = mean(diff(p,1));
dLdpp = diff(L,2)/(dp^2);
[l i] = max(L(2:end - 1) + (dLdpp < 0)*1024);
Cp    = inv(-dLdpp(i));
Ep    = p(i + 1);
 
 
% plot likelihood
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(p,L - min(L))
xlabel('Utility','FontSize',12)
ylabel('Probability','FontSize',12)
title('Log-likelihood','FontSize',16)
axis square
    
% plot posterior
%--------------------------------------------------------------------------
subplot(2,2,4)
pp  = spm_Npdf(p,Ep,Cp);                   % posterior probability
tp  = log(MDP.C);                          % true utilities
tp  = tp - min(tp) + 1;
plot(p,pp), hold on
plot([tp(end) tp(end)],[0 1.2*max(pp)],':'),   hold off
xlabel('Utility','FontSize',12)
ylabel('Probability','FontSize',12)
title('Posterior probability','FontSize',16)
axis square
 
 
return
 
 
 
% expected utility
%==========================================================================
function [ED,EU,PT] = PrEU(MDP)
 
% numerical solution
%--------------------------------------------------------------------------
MDP.plot = 0;
MDP.s = [];
MDP.a = [];
MDP.o = [];
 
ST    = 0;
for i = 1:64
    [P,Q,S] = spm_MDP_game(MDP);
    ST = ST + S(:,end);
end
PT    = ST/sum(ST);
i     = find(PT);
EU    = PT'*log(MDP.C/sum(MDP.C));
ED    = EU - PT(i)'*log(PT(i));
 
return


% notes for precision - action value figure 
%==========================================================================
clf; subplot(2,1,1)
Q  = -8:1/32:0;
W  = 8./(1 - Q);
plot(Q,W)
xlabel('Expected action value','FontSize',12)
ylabel('Expected precision','FontSize',12)
title('Precision and action value','FontSize',16)
grid on
axis square




