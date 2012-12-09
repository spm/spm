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
% particular inversion scheme used here is spm_MDP_select, which uses a
% mean-field approximation between hidden control and hidden states. It is
% assumed that the agent believes that it will select a particular action
% (accept or decline) at a particular time.
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
% $Id: spm_MDP_offer.m 5104 2012-12-09 21:47:51Z karl $
 
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
 
 
% MDP Structure
%==========================================================================
MDP.T = T;                          % process depth (the horizon)
MDP.S = S;                          % initial state
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
 
% solve - an example game (with high offer at t = 10)
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
 
MDP.s    = [1 1 1 1 1 1 1 1 1 2];
MDP.a    = [1 1 1 1 1 1 1 1 1 1];
MDP.plot = gcf;
MDP.N    = 8;
 
[P,Q,S,U,W,da] = spm_MDP_select(MDP);
 
% plot convergence and precision
%--------------------------------------------------------------------------
subplot(2,2,2)
spm_axis tight
subplot(2,2,3)
plot(W)
xlabel('Latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Expected precision','FontSize',16)
spm_axis tight
axis square
 
subplot(2,2,4)
plot(da,'b'),   hold on
plot(diff(da),'r'), hold off
xlabel('Latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Precision dynamics','FontSize',16)
spm_axis tight
axis square
 
% solve - an example game (with low offer at t = 5)
%==========================================================================
spm_figure('GetWin','Figure 2'); clf
 
MDP.s    = [1 1 1 1 3];
MDP.a    = ones(1,T);
MDP.plot = gcf;
 
[P,Q,S,U,W,da] = spm_MDP_select(MDP);
 
% plot convergence and precision
%--------------------------------------------------------------------------
subplot(2,2,2)
spm_axis tight
subplot(2,2,3)
plot(W)
xlabel('Latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Expected precision','FontSize',16)
spm_axis tight
axis square
 
subplot(2,2,4)
plot(da,'b'),   hold on
plot(diff(da),'r'), hold off
xlabel('latency (offers)','FontSize',12)
ylabel('Precision of beliefs','FontSize',12)
title('Precision dynamics','FontSize',16)
spm_axis tight
axis square
 
 
% illustrate dependency on latency w.r.t. parameters
%==========================================================================
spm_figure('GetWin','Figure 3'); clf
 
% probability distribution over time to act: P(1,:) is no action
%--------------------------------------------------------------------------
PrT      = @(P)[1 cumprod(P(1,1:end - 1))].*P(2,:);
MDP.plot = 0;                        % plot convergence
MDP.N    = 4;                        % number of variational iterations
MDP.s    = ones(1,T);                % suppress withdrawal of low offer
MDP.a    = ones(1,T);                % and action
 
% beliefs about final state
%--------------------------------------------------------------------------
DP    = MDP;
p     = linspace(0,8,16);
for i = 1:length(p)
    DP.C     = spm_softmax([1 1 1 p(i) 4]');
    BF       = spm_MDP_select(DP);
    BE       = spm_MDP_select(DP,'EU');  
    PF(i,:)  = BF(2,:);
    PE(i,:)  = BE(2,:);
    DF(i,:)  = PrT([BF [0;1]]);
    DE(i,:)  = PrT([BE [0;1]]);
end
 
% probability of accepting
%--------------------------------------------------------------------------
subplot(2,2,1)
imagesc(1:(T - 1),p,1 - PF)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Conditional divergence','FontSize',16)
axis square
 
% compare with expected utility
%--------------------------------------------------------------------------
subplot(2,2,2)
imagesc(PE)
imagesc(1:(T - 1),p,1 - PE)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Expected utility','FontSize',16)
axis square
 
% distribution of acceptance latencies
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(1:(T - 1),p,1 - DF)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Latency of accepting','FontSize',16)
axis square
 
% compare with expected utility
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(PE)
imagesc(1:(T - 1),p,1 - DE)
xlabel('Latency (offers)','FontSize',12)
ylabel('Utility of low offer','FontSize',12)
title('Latency of accepting','FontSize',16)
axis square
 
 
% the effect of sensory precision and its interaction with memory
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
 
MDP.s    = [1 1 1 1 1 1 2];
MDP.a    = [1 1 1 1 1 1 1];
MDP.o    = [1 1 1 4 1 1 2];
MDP.plot = gcf;
 
MDP.N    = 8;
MDP.K    = 1;
 
[P,Q,S,U,W,da] = spm_MDP_select(MDP);
 
% plot true and inferred states
%--------------------------------------------------------------------------
subplot(4,2,6)
imagesc(1 - Q)
title('Inferred states (K > 0)','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('State','FontSize',12)
 
% precisions
%--------------------------------------------------------------------------
subplot(3,1,1)
plot((1:length(da))/MDP.N,da,'r')
title('Expected decision','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('Precision','FontSize',12)
axis square
 
% now repeat but with no memory
%--------------------------------------------------------------------------
MDP.plot = 0;
MDP.K    = 0;
 
[P,Q,S,U,W,da] = spm_MDP_select(MDP);
 
% plot true and inferred states
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4');
 
subplot(4,2,8)
imagesc(1 - Q)
title('Inferred states (K = 0)','FontSize',16)
xlabel('Time','FontSize',12)
ylabel('State','FontSize',12)
 
% precisions
%--------------------------------------------------------------------------
subplot(3,1,1), hold on
plot((1:length(da))/MDP.N,da,'k'), hold off
axis([1 8 1 4])
 
 
% Changes in uncertainty (Entropy) over successive choices
%==========================================================================
spm_figure('GetWin','Figure 5'); clf
 
% uncertainty about current action
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = ones(1,T);
MDP.o = ones(1,T);
 
MDP.plot = 0;
MDP.K    = 1;
MDP.N    = 4;
 
MDP.C       = spm_softmax([1 1 1 3 4]');
[P,Q,S,U,W] = spm_MDP_select(MDP);
H           = sum(-P.*log(P),1);
 
% action entropy
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(1:length(P),P,'-.'), hold on
plot(1:length(H),H,'.r','MarkerSize',16), hold on
plot(1:length(H),H,':r'), hold off
xlabel('Offer','FontSize',12)
ylabel('Probability and entropy (nats)','FontSize',12)
title('Uncertainty about action','FontSize',16)
axis square
 
% precision
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(1:length(W),W,'.','MarkerSize',16), hold on
plot(1:length(W),W,':'), hold off
xlabel('Offer','FontSize',12)
ylabel('Precision','FontSize',12)
title('Precision','FontSize',16)
axis square
 
 
% simulate multiple trials and record when an offer was accepted
%==========================================================================
spm_figure('GetWin','Figure 6'); clf

% trials with no higher offer
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = [];
MDP.o = [];
 
MDP.plot = 0;
MDP.K    = 1;
MDP.N    = 4;

% trials
%--------------------------------------------------------------------------
for i = 1:256
    [P,Q,S,U] = spm_MDP_select(MDP);
    try
        Y(i)  = find(U(2,:),1);
    end
    fprintf('trial %0.00f\n',i);
end
 
% probability distribution over time to act
%--------------------------------------------------------------------------
MDP.s = ones(1,T);
MDP.a = ones(1,T);

[P,Q,S,U] = spm_MDP_select(MDP);
Py        = PrT(P);
 
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
 
 
% Infer prior beliefs from observed responses (meta-modelling)
%==========================================================================
p     = linspace(1/32,1/4,32);
DP    = MDP;
for i = 1:length(p);
    
    % transition probabilities
    %----------------------------------------------------------------------
    for t = 1:T
        a = Pwin(T,Pa);
        b = Plos(t,p(i));
        DP.B{t,1}(:,1) = [1 - a + a*b - b; a*(1 - b); b; 0; 0];       
    end
    
    % get likelihood for this parameter
    %----------------------------------------------------------------------
    P     = spm_MDP_select(DP);
    Py    = PrT(P);
    L(i)  = sum(log(Py(Y)));
    
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
xlabel('Latency','FontSize',12)
ylabel('Probabaility','FontSize',12)
title('Log-likelihood','FontSize',16)
axis square
    
% plot posterior
%--------------------------------------------------------------------------
subplot(2,2,4)
pp    = spm_Npdf(p,Ep,Cp);
plot(p,pp), hold on
plot([Pb Pb],[0 1.2*max(pp)],':'),   hold off
xlabel('latency','FontSize',12)
ylabel('Probability','FontSize',12)
title('Posterior probability','FontSize',16)
axis square
 
 
return
 
 
% expected utility
%==========================================================================
function [EC,PT] = PrE(MDP)
 
% numerical solution
%--------------------------------------------------------------------------
MDP.plot = 0;
 
ST    = 0;
for i = 1:256
    [P,Q,S] = spm_MDP_select(MDP,1);
    ST = ST + S(:,end);
end
PT    = ST/sum(ST);
EC    = PT'*log(MDP.C/sum(MDP.C));
 
return
