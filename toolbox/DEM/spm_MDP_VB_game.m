function Q = spm_MDP_VB_game(MDP)
% auxiliary plotting routine for spm_MDP_VB - multiple trials
% FORMAT Q = spm_MDP_VB_game(MDP)
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
% MDP.un  = un    - simulated neuronal encoding of hidden states
% MDP.xn  = Xn    - simulated neuronal encoding of policies
% MDP.wn  = wn    - simulated neuronal encoding of precision
% MDP.da  = dn    - simulated dopamine responses (deconvolved)
% MDP.rt  = rt    - simulated dopamine responses (deconvolved)
%
% returns summary of performance:
%
%     Q.X  = x    - expected hidden states
%     Q.R  = u    - final policy expectations
%     Q.S  = s    - initial hidden states
%     Q.O  = o    - final outcomes
%     Q.p  = p    - performance
%     Q.q  = q    - reaction times
%
% please see spm_MDP_VB
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_VB_game.m 7307 2018-05-08 09:44:04Z karl $

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
if iscell(MDP(1).X)
    Nf = numel(MDP(1).B);                 % number of hidden state factors
    Ng = numel(MDP(1).A);                 % number of outcome factors
else
    Nf = 1;
    Ng = 1;
end

% graphics
%==========================================================================
Nt    = length(MDP);               % number of trials
Ne    = size(MDP(1).V,1) + 1;      % number of epochs per trial
Np    = size(MDP(1).V,2) + 1;      % number of policies
for i = 1:Nt
    
    % assemble expectations of hidden states and outcomes
    %----------------------------------------------------------------------
    for j = 1:Ne
        for k = 1:Ne
            for f = 1:Nf
                try
                    x{f}{i,1}{k,j} = gradient(MDP(i).xn{f}(:,:,j,k)')';
                catch
                    x{f}{i,1}{k,j} = gradient(MDP(i).xn(:,:,j,k)')';
                end
            end
        end
    end
    s(:,i) = MDP(i).s(:,1);
    o(:,i) = MDP(i).o(:,end);
    u(:,i) = MDP(i).R(:,end);
    w(:,i) = mean(MDP(i).dn,2);
    
    
    % assemble context learning
    %----------------------------------------------------------------------
    for f = 1:Nf
        try
            try
                D = MDP(i).d{f};
            catch
                D = MDP(i).D{f};
            end
        catch
            try
                D = MDP(i).d;
            catch
                D = MDP(i).D;
            end
        end
        d{f}(:,i) = D/sum(D);
    end
    
    % assemble performance
    %----------------------------------------------------------------------
    p(i)  = 0;
    for g = 1:Ng
        try
            U = spm_softmax(MDP(i).C{g});
        catch
            U = spm_softmax(MDP(i).C);
        end
        for t = 1:Ne
            p(i) = p(i) + log(U(MDP(i).o(g,t),t))/Ne;
        end
    end
    q(i)   = sum(MDP(i).rt(2:end));
    
end

% assemble output structure if required
%--------------------------------------------------------------------------
if nargout
    Q.X  = x;            % expected hidden states
    Q.R  = u;            % final policy expectations
    Q.S  = s;            % inital hidden states
    Q.O  = o;            % final outcomes
    Q.p  = p;            % performance
    Q.q  = q;            % reaction times
    return
end


% Initial states and expected policies (habit in red)
%--------------------------------------------------------------------------
col   = {'r.','g.','b.','c.','m.','k.'};
t     = 1:Nt;
subplot(6,1,1)
if Nt < 64
    MarkerSize = 24;
else
    MarkerSize = 16;
end
image(64*(1 - u)),  hold on
for f = 1:Nf
    for i = 1:max(s(f,:))
        j = find(s(f,:) == i);
        plot(t(j),j - j + f,col{rem(i - 1,6)+ 1},'MarkerSize',MarkerSize)
    end
end
try
    plot(Np*(1 - u(Np,:)),'r')
end
try
    E = spm_softmax(spm_cat({MDP.e}));
    plot(Np*(1 - E(end,:)),'r:')
end
title('Initial state and policy selection')
xlabel('Trial'),ylabel('Policy'), hold off


% Performance
%--------------------------------------------------------------------------
q     = q - mean(q);
q     = q/std(q);
subplot(6,1,2), bar(p,'k'),   hold on
plot(q,'.c','MarkerSize',16), hold on
plot(q,':c')
for g = 1:Ng
    for i = 1:max(o(g,:))
        j = find(o(g,:) == i);
        plot(t(j),j - j + 3 + g,col{rem(i - 1,6)+ 1},'MarkerSize',MarkerSize)
    end
end
title('Final outcome, performance and reaction times')
ylabel('Expected utility'), spm_axis tight, hold off, box off

% Initial states (context)
%--------------------------------------------------------------------------
subplot(6,1,3)
col   = {'r','b','g','c','m','k','r','b','g','c','m','k'};
for f = 1:Nf
    if Nf > 1
        plot(spm_cat(x{f}),col{f}), hold on
    else
        plot(spm_cat(x{f}))
    end
end
title('State estimation (ERPs)'), ylabel('Response'), 
spm_axis tight, hold off, box off

% Precision (dopamine)
%--------------------------------------------------------------------------
subplot(6,1,4)
w   = spm_vec(w);
if Nt > 8
    fill([1 1:length(w) length(w)],[0; w.*(w > 0); 0],'k'), hold on
    fill([1 1:length(w) length(w)],[0; w.*(w < 0); 0],'k'), hold off
else
    bar(w,1.1,'k')
end
title('Precision (dopamine)')
ylabel('Precision','FontSize',12), spm_axis tight, box off
YLim = get(gca,'YLim'); YLim(1) = 0; set(gca,'YLim',YLim);
set(gca,'XTickLabel',{});

% learning - D
%--------------------------------------------------------------------------
for f = 1:Nf
    subplot(6*Nf,1,Nf*4 + f), image(64*(1 - d{f}))
    if f < 2
        title('Context Learning')
    end
    set(gca,'XTick',1:Nt);
    if f < Nf
        set(gca,'XTickLabel',{});
    end
    set(gca,'YTick',1);
    try
        set(gca,'YTickLabel',MDP(1).label.factor{f});
    end
    try
        set(gca,'YTickLabel',MDP(1).Bname{f});
    end
    
    
end
if isfield(MDP(1),'c')
    title('Learning (C and D)')
else
    return
end

% Habit learning
%--------------------------------------------------------------------------
k     = round(linspace(1,Nt,6));
for j = 1:length(k)
    h = MDP(k(j)).c;
    h = h*diag(1./sum(h));
    subplot(6,6,30 + j), image(64*(1 - h))
    axis image
end
