function spm_MDP_VB_trial(MDP)
% auxiliary plotting routine for spm_MDP_VB - single trial
% FORMAT spm_MDP_VB_trial(MDP)
%
% MDP.P(M,T)      - probability of emitting action 1,...,M at time 1,...,T
% MDP.Q(N,T)      - an array of conditional (posterior) expectations over
%                   N hidden states and time 1,...,T
% MDP.X           - and Bayesian model averages over policies
% MDP.R           - conditional expectations over policies
% MDP.o           - outcomes at time 1,...,T
% MDP.s           - states at time 1,...,T
% MDP.u           - action at time 1,...,T
%
% MDP.un  = un;   - simulated neuronal encoding of hidden states
% MDP.xn  = Xn;   - simulated neuronal encoding of policies
% MDP.wn  = wn;   - simulated neuronal encoding of precision
% MDP.da  = dn;   - simulated dopamine responses (deconvolved)
% MDP.rt  = rt;   - simulated reaction times
%
% please see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_trial.m 6652 2015-12-21 10:51:54Z karl $

% graphics
%==========================================================================

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
if iscell(MDP.X)
    Nf = numel(MDP.B);                 % number of hidden state factors
    Ng = numel(MDP.A);                 % number of outcome factors
    X  = MDP.X;
    C  = MDP.C;
else
    Nf = 1;
    Ng = 1;
    X  = {MDP.X};
    C  = {MDP.C};
end



% posterior beliefs about hidden states
%--------------------------------------------------------------------------
for f  = 1:Nf
    subplot(3*Nf,2,(f - 1)*2 + 1)
    image(64*(1 - X{f})), hold on
    if size(X{f},1) > 128
        spm_spy(X{f},16,1)
    end
    plot(MDP.s(f,:),'.c','MarkerSize',16), hold off
    try
        title(sprintf('Hidden states - %s',MDP.Bname{f}));
    catch
        if f < 2, title('Hidden states'); end
    end
    xlabel('time')
    ylabel('hidden state')
end

% posterior beliefs about control states
%--------------------------------------------------------------------------
for f  = 1:Nf
    subplot(3*Nf,2,f*2)
    if Nf == 1
        P = MDP.P;
    elseif Nf == 2
        if f == 1
            P = squeeze(mean(MDP.P,2));
        else
            P = squeeze(mean(MDP.P,1))';
        end
    elseif Nf == 3
        if f == 1
            P = squeeze(mean(squeeze(mean(MDP.P,3),2)));
        elseif f == 2
            P = squeeze(mean(squeeze(mean(MDP.P,3),1)));
        else
            P = squeeze(mean(squeeze(mean(MDP.P,1),2)));
        end
    end
    
    image(64*(1 - P)), hold on
    plot(MDP.u(f,:),'.c','MarkerSize',16), hold off
    try
        title(sprintf('Inferred and selected action - %s',MDP.Bname{f}));
    catch
        if f < 2, title('Inferred and selected action'); end
    end
    xlabel('time')
    ylabel('action')
end

% policies
%--------------------------------------------------------------------------
for f  = 1:Nf
    subplot(3*Nf,2,(Nf + f - 1)*2 + 1)
    imagesc(MDP.V(:,:,f)')
    try
        title(sprintf('Allowable policies - %s',MDP.Bname{f}));
    catch
        if f < 2, title('Allowable policies'); end
    end
    ylabel('policy')
    xlabel('time')
end

% expectations over policies
%--------------------------------------------------------------------------
subplot(3,2,4)
image(64*(1 - MDP.un))
title('Posterior probability','FontSize',14)
ylabel('policy','FontSize',12)
xlabel('updates','FontSize',12)

% sample (observation) and preferences
%--------------------------------------------------------------------------
for g  = 1:Ng
    subplot(3*Ng,2,(2*Ng + g - 1)*2 + 1)
    if size(C{g},1) > 128
        spm_spy(C{g},16,1), hold on
    else
        imagesc(1 - C{g}), hold on
    end
    plot(MDP.o(g,:),'.c','MarkerSize',16), hold off
    try
        title(sprintf('Outcomes and preferences - %s',MDP.Aname{g}));
    catch
        if f < 2, title('Outcomes and preferences'); end
    end
    xlabel('time')
    ylabel('outcome')
end

% expected precision
%--------------------------------------------------------------------------
subplot(3,2,6), hold on
if size(MDP.dn,2) > 1
    plot(MDP.dn,'r:'), plot(MDP.wn,'k'), hold off
else
    bar(MDP.dn,1.1,'k'),   plot(MDP.wn,'k'), hold off
end
title('Expected precision (dopamine)','FontSize',14)
xlabel('updates','FontSize',12)
ylabel('precision','FontSize',12)
spm_axis tight
drawnow
