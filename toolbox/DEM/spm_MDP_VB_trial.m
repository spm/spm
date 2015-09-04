function spm_MDP_VB_trial(MDP)
% auxiliary plotting routine for spm_MDP_VB - single trial
% FORMAT spm_MDP_VB_trial(MDP)
%
% MDP.P(M,T)      - probability of emitting action 1,...,M at time 1,...,T
% MDP.Q(N,T)      - an array of conditional (posterior) expectations over
%                   N hidden states and time 1,...,T
% MDP.X           - and Bayesian model averages over policies
% MDP.R           - conditional expectations over policies
% MDP.O(O,T)      - a sparse matrix encoding outcomes at time 1,...,T
% MDP.S(N,T)      - a sparse matrix encoding states at time 1,...,T
% MDP.U(M,T)      - a sparse matrix encoding action at time 1,...,T
% MDP.W(1,T)      - posterior expectations of precision
%
% MDP.un  = un;   - simulated neuronal encoding of hidden states
% MDP.xn  = Xn;   - simulated neuronal encoding of policies
% MDP.wn  = wn;   - simulated neuronal encoding of precision
% MDP.da  = dn;   - simulated dopamine responses (deconvolved)
% MDP.rt  = rt;   - simulated dopamine responses (deconvolved)
%
% please see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_VB_trial.m 6539 2015-09-04 08:47:25Z karl $
 
% graphics
%==========================================================================
 
% posterior beliefs about hidden states
%--------------------------------------------------------------------------
subplot(3,2,1)
image(64*(1 - MDP.X)),hold on
if size(MDP.X,1) > 128
    spm_spy(MDP.X,16,1)
end
[s,t] = find(MDP.S); plot(t,s,'.c','MarkerSize',16), hold off
title('Hidden states (and utility)','FontSize',14)
xlabel('trial','FontSize',12)
ylabel('hidden state','FontSize',12)
 
% posterior beliefs about control states
%--------------------------------------------------------------------------
subplot(3,2,2), image(64*(1 - MDP.P)), hold on
[a,t] = find(MDP.U); plot(t,a,'.c','MarkerSize',16), hold off
title('Inferred and selected action','FontSize',14)
xlabel('trial','FontSize',12)
ylabel('action','FontSize',12)
 
% policies
%--------------------------------------------------------------------------
subplot(3,2,3)
imagesc(MDP.V')
title('allowable policies','FontSize',14)
ylabel('policy','FontSize',12)
xlabel('trial','FontSize',12)
 
% expectations over policies
%--------------------------------------------------------------------------
subplot(3,2,4)
image(64*(1 - MDP.un))
title('Posterior probability','FontSize',14)
ylabel('policy','FontSize',12)
xlabel('updates','FontSize',12)
 
% sample (observation)
%--------------------------------------------------------------------------
subplot(3,2,5)
if size(MDP.O,1) > 128
    spm_spy(MDP.O,16,1)
else
    imagesc(1 - MDP.O)
end
title('Observed states','FontSize',14)
xlabel('trial','FontSize',12)
ylabel('outcome','FontSize',12)
 
% expected action
%--------------------------------------------------------------------------
subplot(3,2,6), hold on
bar(4*MDP.dn,'c'), plot(MDP.wn,'k'), hold off
title('Expected precision (dopamine)','FontSize',14)
xlabel('updates','FontSize',12)
ylabel('precision','FontSize',12)
spm_axis tight
drawnow
