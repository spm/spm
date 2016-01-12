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
% $Id: spm_MDP_VB_trial.m 6672 2016-01-12 12:28:31Z karl $

% graphics
%==========================================================================

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
if iscell(MDP.X)
    Nf = numel(MDP.B);                 % number of hidden state factors
    Ng = numel(MDP.A);                 % number of outcome factors
    X  = MDP.X;
    C  = MDP.C;
    for f = 1:Nf
        Nu(f) = size(MDP.B{f},3) > 1;
    end
else
    Nf = 1;
    Ng = 1;
    Nu = 1;
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
    if f == Nf, xlabel('time'), end
    ylabel('hidden state')
end

% posterior beliefs about control states
%--------------------------------------------------------------------------
Nu     = find(Nu);
Np     = length(Nu);
for f  = 1:Np
    subplot(3*Np,2,f*2)
    P = MDP.P;
    if Nf > 1
        ind     = 1:Nf;
        for dim = 1:Nf
            if dim ~= ind(Nu(f));
                P = sum(P,dim);
            end
        end
        P = squeeze(P);
    end
    
    % display
    %----------------------------------------------------------------------
    image(64*(1 - P)), hold on
    plot(MDP.u(Nu(f),:),'.c','MarkerSize',16), hold off
    try
        title(sprintf('Inferred and selected action - %s',MDP.Bname{Nu(f)}));
    catch
        if f < 2, title('Inferred and selected action'); end
    end
    if f == Np, xlabel('time'), end
    ylabel('action')
end

% policies
%--------------------------------------------------------------------------
for f  = 1:Np
    subplot(3*Np,2,(Np + f - 1)*2 + 1)
    imagesc(MDP.V(:,:,Nu(f))')
    try
        title(sprintf('Allowable policies - %s',MDP.Bname{Nu(f)}));
    catch
        if f < 2, title('Allowable policies'); end
    end
    if Np == 1, xlabel('time'), end
    ylabel('policy')
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
    if g == Ng, xlabel('time'), end
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
