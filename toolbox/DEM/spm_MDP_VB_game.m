function spm_MDP_VB_game(MDP)
% auxiliary plotting routine for spm_MDP_VB - multiple trials
% FORMAT spm_MDP_VB_game(MDP)
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
% $Id: spm_MDP_VB_game.m 6529 2015-08-21 13:27:38Z karl $

% graphics
%==========================================================================
col   = {'r.','b.','g.','c.','m.','k.'};
Nt    = length(MDP);               % number of trials
NT    = size(MDP(1).V,1) + 1;      % number of transitions
Np    = size(MDP(1).V,2) + 1;      % number of policies
for i = 1:Nt
    for j = 1:NT
        for k = 1:NT
            xi{k,j} = MDP(i).xn(:,:,j,k);
        end
    end
    x{i,1} = xi;
    u(:,i) = MDP(i).R(:,end);
    s(:,i) = MDP(i).S(:,1);
    o(:,i) = MDP(i).O(:,end);
    d(:,i) = MDP(i).d/sum(MDP(i).d);
    w(:,i) = MDP(i).dn;
    p(i)   = trace(MDP(i).C'*MDP(i).S)/NT;
    p(i)   = trace(log(MDP(i).A*spm_softmax(MDP(i).C))'*MDP(i).O)/NT;
    q(i)   = sum(MDP(i).rt(2:end));
end

% Initial tates and expected policies (habit in red)
%--------------------------------------------------------------------------
subplot(6,1,1), 
[s,t] = find(s);
if Nt < 64
    MarkerSize = 32;
else
    MarkerSize = 16;
end
image(64*(1 - u)),  hold on
for i = 1:max(s)
    j = find(s == i);
    plot(t(j),j - j + 1,col{rem(i - 1,6)+ 1},'MarkerSize',MarkerSize), hold on
end
plot(Np*(1 - u(end,:)),'r'), hold off
title('Polcy selection and inital state','FontSize',16)
xlabel('Trial','FontSize',12),ylabel('Policy','FontSize',12)

% Performance
%--------------------------------------------------------------------------
q     = q - mean(q);
q     = q/std(q);
[o,t] = find(o);
subplot(6,1,2), bar(p,'k'),   hold on
plot(q,'.c','MarkerSize',16), hold on
plot(q,':c'),                 hold on
for i = 1:max(o)
    j = find(o == i);
    plot(t(j),j - j + 4,col{rem(i - 1,6)+ 1},'MarkerSize',MarkerSize), hold on
end
title('Final outcome, performance and reaction times','FontSize',16)
ylabel('Expected utility','FontSize',12), spm_axis tight

% Initial states (context)
%--------------------------------------------------------------------------
subplot(6,1,3)
plot(spm_cat(x))
title('State estimation (ERPs)','FontSize',16)
ylabel('Response','FontSize',12), spm_axis tight

% Precision (dopamine)
%--------------------------------------------------------------------------
subplot(6,1,4), bar(spm_vec(w),'k')
title('Precision (dopamine)','FontSize',16)
ylabel('Precision','FontSize',12), spm_axis tight

% learning - D
%--------------------------------------------------------------------------
subplot(6,1,5), image(64*(1 - d))
title('Learning (D and H)','FontSize',16)
ylabel('Hidden state','FontSize',12)

% Habit learning
%--------------------------------------------------------------------------
k     = round(linspace(1,Nt,6));
for j = 1:length(k)
    h = MDP(k(j)).c;
    h = h*diag(1./sum(h));
    subplot(6,6,30 + j), image(64*(1 - h))
    axis image
end
