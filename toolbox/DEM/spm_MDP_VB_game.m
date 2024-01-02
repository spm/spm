function Q = spm_MDP_VB_game(MDP)
% auxiliary plotting routine for spm_MDP_VB - multiple trials
% FORMAT Q = spm_MDP_VB_game(MDP)
%
% MDP.P(M,T)      - probability of emitting action 1,...,M at time 1,...,T
% MDP.X           - conditional expectations over hidden states
% MDP.R           - conditional expectations over policies
% MDP.O(O,T)      - a sparse matrix encoding outcomes at time 1,...,T
% MDP.S(N,T)      - a sparse matrix encoding states at time 1,...,T
% MDP.U(M,T)      - a sparse matrix encoding action at time 1,...,T
% MDP.W(1,T)      - posterior expectations of precision
%
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

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

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
Nt    = numel(MDP);                      % number of trials
Ne    = MDP(1).T;                        % number of outcomes per trial
Np    = size(MDP(1).V,2);                % number of policies
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
    
    % fill-in selected policies if missing
    %----------------------------------------------------------------------
    try
        MDP(i).v;
    catch
        MDP(i).v = [];
    end
    
    s(:,i) = MDP(i).s(:,1);
    o(:,i) = MDP(i).o(:,end);
    u{1,i} = MDP(i).R;
    v{1,i} = MDP(i).v;
    w(:,i) = MDP(i).dn;
    
    
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
            try
                p(i) = p(i) + log(U(MDP(i).o(g,t),t))/Ne;
            catch
                p(i) = p(i) + log(U(MDP(i).o(g,t),1))/Ne;
            end
        end
    end
    
    % reaction time
    %----------------------------------------------------------------------
    try
        q(i) = sum(MDP(i).rt(2:end));
    catch
        q(i) = 0;
    end
    
end

% assemble output structure if required
%--------------------------------------------------------------------------
if nargout
    Q.X  = x;            % expected hidden states
    Q.R  = u;            % final policy expectations
    Q.S  = s;            % initial hidden states
    Q.O  = o;            % final outcomes
    Q.p  = p;            % performance
    Q.q  = q;            % reaction times
    return
end


% Initial states and expected policies (habit in red)
%--------------------------------------------------------------------------
col   = {'r.','g.','b.','c.','m.','k.'};
u     = spm_cat(u);
v     = spm_cat(v);
t     = 1:Nt;
subplot(6,1,1)
if Nt < 64
    MarkerSize = 24;
else
    MarkerSize = 16;
end
image(t,1:Np,64*(1 - u)),  hold on
try
    plot(linspace(1,Nt,numel(v)),v,'c.','MarkerSize',16)
end
for f = 1:Nf
    for i = 1:max(s(f,:))
        j = find(s(f,:) == i);
        plot(t(j),j - j + f - Nf,col{rem(i - 1,6)+ 1},'MarkerSize',MarkerSize)
    end
end
set(gca,'YLim',[-Nf,Np] + 1/2), box off
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
