function MDP = DEM_demo_MDP_rule
% Demo of active inference for visual salience
%__________________________________________________________________________
%
% This routine simulates a crude form of consciousness using active
% inference and structure learning to vitiate ignorance or nescience. This
% entails learning the hyper parameters of causal structure generating
% outcomes and then using Bayesian model reduction (during sleep) to
% minimise complexity.
%
% We first set up an abstract problem in which an agent has to respond
% according to rules (identify the correct colour depending upon one of
% three rules that are specified by the colour of a cue in the centre of
% vision). If the rule is centre, the colour is always green; however, if
% the colour of the Centre cue is red, the correct colour is on the left
% (and on the right if the queue is blue). Simulations are provided when
% the agent knows the rules. This is then repeated in the absence
% (nescience) of any knowledge about the rules to see if the agent can
% learn causal structure through Bayesian belief updating of the likelihood
% array (A).
%
% We then consider the improvement in performance (in terms of variational
% free energy, its constituent parts and performance) following Bayesian
% model reduction of the likelihood model (heuristically, like slow wave
% sleep), followed by a restitution of posterior beliefs during fictive
% active inference (as in REM sleep). Finally, we address the communication
% of the implicit  structure learning to a conspecific or child to
% demonstrate the improvement under instruction.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_rule.m 7319 2018-05-29 09:33:01Z karl $

% set up and preliminaries
%==========================================================================

% specify the (rule-based) generative process (environment)
%==========================================================================
rng('default')

% prior beliefs about initial states (in terms of counts_: D
%--------------------------------------------------------------------------
D{1} = [1 1 1]';        % rule:   {'left','centre','right'}
D{2} = [1 1 1]';        % what:   {'red','green','blue'}
D{3} = [0 0 0 1]';      % where:  {'left','centre','right','null'}
D{4} = [0 0 0 1]';      % report: {'red','green','blue','undecided'}

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
Ns    = zeros(1,Nf);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
for f1 = 1:Ns(1)                     % rule
    for f2 = 1:Ns(2)                 % correct colour
        for f3 = 1:Ns(3)             % location of fixation
            for f4 = 1:Ns(3)         % decision
                
                % A{1} what: {'red','green','blue','null'}
                %==========================================================
                if f3 == 4, A{1}(4,f1,f2,f3,f4)  = 1;  end
                if f3 == 2, A{1}(f1,f1,f2,f3,f4) = 1;  end
                if f3 == 1
                    if f1 == 1
                        A{1}(f2,f1,f2,f3,f4) = 1;
                    else
                        A{1}(1:3,f1,f2,f3,f4) = 1/3;
                    end
                end
                if f3 == 3
                    if f1 == 3
                        A{1}(f2,f1,f2,f3,f4) = 1;
                    else
                        A{1}(1:3,f1,f2,f3,f4) = 1/3;
                    end
                end
                
                % A{2} where: {'left','centre','right','null'}
                %----------------------------------------------------------
                A{2}(f3,f1,f2,f3,f4) = 1;
                
                % A{3} feedback: {'null','right','wrong'}
                %----------------------------------------------------------
                if f4 == 4,
                    A{3}(1,f1,f2,f3,f4) = 1;                             % undecided
                else
                    if f1 == 2 && f4 == 2,  A{3}(2,f1,f2,f3,f4) = 1; end % right
                    if f1 == 2 && f4 ~= 2,  A{3}(3,f1,f2,f3,f4) = 1; end % wrong
                    if f1 == 1 && f4 == f2, A{3}(2,f1,f2,f3,f4) = 1; end % right
                    if f1 == 1 && f4 ~= f2, A{3}(3,f1,f2,f3,f4) = 1; end % wrong
                    if f1 == 3 && f4 == f2, A{3}(2,f1,f2,f3,f4) = 1; end % right
                    if f1 == 3 && f4 ~= f2, A{3}(3,f1,f2,f3,f4) = 1; end % wrong
                end
            end
        end
    end
end
Ng    = numel(A);
for g = 1:Ng
    No(g) = size(A{g},1);
    A{g}  = double(A{g});
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% control states B(3): where {'left','centre','right','null'}
%--------------------------------------------------------------------------
for k = 1:Ns(3)
    B{3}(:,:,k) = 0;
    B{3}(k,:,k) = 1;
end

% control states B(4): report {'red','green','blue','undecided'}
%--------------------------------------------------------------------------
for k = 1:Ns(4)
    B{4}(:,:,k) = 0;
    B{4}(k,:,k) = 1;
end

% allowable policies (specified as the next action) U
%--------------------------------------------------------------------------
U(1,1,:)  = [1 1 1 4]';         % sample left
U(1,2,:)  = [1 1 2 4]';         % sample centre
U(1,3,:)  = [1 1 3 4]';         % sample right
U(1,4,:)  = [1 1 4 1]';         % return to centre and report red
U(1,5,:)  = [1 1 4 2]';         % return to centre and report green
U(1,6,:)  = [1 1 4 3]';         % return to centre and report blue


% priors: (utility) C; the agent expects to avoid mistakes
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end
% and expects itself to make a decision after the fifth observation
%--------------------------------------------------------------------------
C{3} = [ 0  0  0  0 -8 -8;
         0  0  0  0  0  0;
        -4 -4 -4 -4 -4 -4];

% MDP Structure
%--------------------------------------------------------------------------
mdp.T = size(C{3},2);           % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states

mdp.Aname = {'what','where','feedback'};
mdp.Bname = {'rule','colour','where','decision'};

% illustrate a single trial
%==========================================================================
mdp  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

% illustrate behaviour
%--------------------------------------------------------------------------
subplot(3,2,3)
spm_MDP_rule_plot(MDP)

% illustrate a sequence of trials
%==========================================================================
clear MDP

% create structure array
%--------------------------------------------------------------------------
N     = 8;
for i = 1:N
    MDP(i) = mdp;
end


% Solve an example sequence under different initial states
%==========================================================================
MDP   = spm_MDP_VB_X(MDP);

% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_game(MDP);

[F,Fu] = spm_MDP_F(MDP);

subplot(6,1,4), plot(1:N,F),  xlabel('trial'), spm_axis tight, title('Free energy','Fontsize',16)
subplot(6,1,5), plot(1:N,Fu), xlabel('trial'), spm_axis tight, title('Confidence','Fontsize',16)
subplot(6,1,6), spm_MDP_plot_moves(MDP),       spm_axis tight
title('Action ','Fontsize',16), legend({'saccades','hits'})


% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP);

% illustrate behaviour in more detail
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf;
for i = 1:min(N,15)
    subplot(5,3,i), spm_MDP_rule_plot(MDP(i));
    axis square
end


% repeat with rule learning
%==========================================================================

% retain knowledge about the cue but remove knowledge about the rules
%--------------------------------------------------------------------------
for f1 = 1:Ns(1)                     % rule
    for f2 = 1:Ns(2)                 % correct colour
        for f3 = 1:Ns(3)             % location of fixation
            for f4 = 1:Ns(3)         % decision
                
                % A{1} what: {'red','green','blue','null'}
                %==========================================================
                if f3 == 4, a{1}(4,f1,f2,f3,f4)  = 128;  end
                if f3 == 2, a{1}(f1,f1,f2,f3,f4) = 128;  end
                if f3 == 1
                    a{1}(1:3,f1,f2,f3,f4) = 1;
                end
                if f3 == 3
                    a{1}(1:3,f1,f2,f3,f4) = 1;
                end
            end
        end
    end
end

a{2}   = A{2}*128;
a{3}   = A{3}*128;
mda    = mdp;
mda.a  = a;
mda.a0 = a;

% create structure array
%--------------------------------------------------------------------------
clear MDP
N     = 32;
for i = 1:N
    MDP(i) = mda;
end


% Solve - an example sequence
%==========================================================================
rng('default')
MDP  = spm_MDP_VB_X(MDP);

% show responses to a difficult (right) trial before and after learning
%--------------------------------------------------------------------------
for i = 1:N
    s(:,i) = MDP(i).s(:,1);
end
i = find(ismember(s(1:2,:)',s(1:2,1)','rows'));

spm_figure('GetWin','Figure 6 - before');
spm_MDP_VB_LFP(MDP(i(1)),[],2);
subplot(3,2,2), set(gca,'YLim',[-.1 1])
subplot(3,2,4), set(gca,'YLim',[-.1 .2])
subplot(3,2,3), spm_MDP_rule_plot(MDP(i(1)))
title(['trial ' num2str(i(1))],'FontSize',16)

spm_figure('GetWin','Figure 6 - after' );
spm_MDP_VB_LFP(MDP(i(end)),[],2);
subplot(3,2,2), set(gca,'YLim',[-.1 1]);
subplot(3,2,4), set(gca,'YLim',[-.1 .2])
subplot(3,2,3), spm_MDP_rule_plot(MDP(i(end)))
title(['trial ' num2str(i(end))],'FontSize',16)

% show learning effects in terms of underlying uncertainty
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7 - before'); clf; spm_MDP_VB_trial(MDP(i(1)));
spm_figure('GetWin','Figure 7 - after' ); clf; spm_MDP_VB_trial(MDP(i(end)));


% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
n  = min(N,4);
spm_figure('GetWin','Figure 4a'); clf
spm_MDP_VB_LFP(MDP(1:n),[],2);
for i = 1:n
    subplot(4,n,i + 3*n), spm_MDP_rule_plot(MDP(i));
    axis square
end


% show trial-by-trial measures in terms of free energy terms
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 8'); clf;
spm_MDP_VB_game(MDP);


[F,Fu] = spm_MDP_F(MDP);

subplot(6,1,4), plot(1:N,F),  xlabel('trial'), spm_axis tight, title('Free energy','Fontsize',16)
subplot(6,1,5), plot(1:N,Fu), xlabel('trial'), spm_axis tight, title('Confidence','Fontsize',16)
subplot(6,1,6), spm_MDP_plot_moves(MDP),       spm_axis tight
title('Action ','Fontsize',16), legend({'saccades','hits'})


% illustrate Bayesian model reduction
%--------------------------------------------------------------------------
OPTIONS.g = 1;
OPTIONS.f = 2;
OPTIONS.T = 3;
OPTIONS.m = @(i,i1,i2,i3,i4) i == i2;

for n = 1:N
    sdp{n} = spm_MDP_VB_sleep(MDP(n),OPTIONS);
    cor(n) = corr(spm_vec(MDP(n).A{1} > 0),spm_vec(sdp{n}.a{1}) > 0);
end
[m,n] = max(cor);

spm_figure('GetWin','Figure 9'); clf; str = sprintf('Sleep (trial %d)',n);
for n = n
    subplot(3,2,1), spm_MDP_A_plot(MDP(n).A);
    subplot(3,2,2), spm_MDP_A_plot(MDP(n).a0), title('Before','Fontsize',16)
    subplot(3,2,3), spm_MDP_A_plot(MDP(n).a),  title('After', 'Fontsize',16)
    subplot(3,2,4), spm_MDP_A_plot(sdp{n}.a),  title(str,     'Fontsize',16)
end

% Bayesian model reduction with dreaming
%--------------------------------------------------------------------------
OPTIONS.o = {MDP.o};
rdp       = spm_MDP_VB_sleep(MDP(n),OPTIONS);
sdp       = sdp{n};


% illustrate benefits of sleep (with and without dreaming)
%==========================================================================
RDP   = MDP(n:end);
RDP   = rmfield(RDP,{'o','u'});
for i = 1:length(RDP)
    RDP(i).s = RDP(i).s(:,1);
end
SDP    = RDP;
SDP(1).a = sdp.a;
RDP(1).a = rdp.a;
SDP    = spm_MDP_VB_X(SDP);
RDP    = spm_MDP_VB_X(RDP);

% plot performance and without sleep
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 10'); clf;

[F,Fu] = spm_MDP_F(MDP);

subplot(4,1,1),plot(1:N,F, ':'),  xlabel('trial'), spm_axis tight, hold on
subplot(4,1,2),plot(1:N,Fu,':'),  xlabel('trial'), spm_axis tight, hold on

% and after SWS sleep
%--------------------------------------------------------------------------
[F,Fu] = spm_MDP_F(SDP);

subplot(4,1,1),plot(n:N,F ,'-.'), xlabel('trial'), spm_axis tight, hold on
subplot(4,1,2),plot(n:N,Fu,'-.'), xlabel('trial'), spm_axis tight, hold on

% and after REM sleep
%--------------------------------------------------------------------------
[F,Fu] = spm_MDP_F(RDP);

subplot(4,1,1),plot(n:N,F ,'-'),  xlabel('trial'), spm_axis tight, hold on
title('Free energy (with sleep)','Fontsize',16)
subplot(4,1,2),plot(n:N,Fu,'-'),  xlabel('trial'), spm_axis tight, hold on
title('Confidence (with sleep)','Fontsize',16)

% and correct responses
%--------------------------------------------------------------------------
hit = @(MDP) any(MDP.o(3,:) == 2) & ~any(MDP.o(3,:) == 3);
k   = min(Fu);
hold on
for i = 1:numel(RDP)
    j = i + n - 1;
    if hit(SDP(i)), plot(j,k,'hr','MarkerSize',8); end
end
for i = 1:(n - 1)
    if hit(MDP(i)), plot(i,k,'hb','MarkerSize',8); end
end
hold off


% illustrate the benefits of instruction (communication of model)
%==========================================================================

% create instructed agent
%--------------------------------------------------------------------------
n     = 1:16;
NDP   = MDP(n);
NDP   = rmfield(NDP,{'o','u'});
for i = 1:length(NDP)
    NDP(i).s = NDP(i).s(:,1);
    NDP(i).a = sdp.a0;
end
NDP   = spm_MDP_VB_X(NDP);

% plot performance and without received wisdom
%--------------------------------------------------------------------------
[F,Fu] = spm_MDP_F(MDP(n));

subplot(4,2,5),plot(n,F ,':'), xlabel('trial'), spm_axis tight, hold on
subplot(4,2,6),plot(n,Fu,':'), xlabel('trial'), spm_axis tight, hold on

% and after instruction
%--------------------------------------------------------------------------
[F,Fu] = spm_MDP_F(NDP(n));
k      = min(Fu);

subplot(4,2,5),plot(n,F),       xlabel('trial'), spm_axis tight
title('Free energy (instruction)','Fontsize',16)
subplot(4,2,6),plot(n,Fu),      xlabel('trial'), spm_axis tight
title('Confidence (instruction)','Fontsize',16)

% and correct responses
%--------------------------------------------------------------------------
for i = n
    hold on
    if hit(NDP(i)); plot(i,k,'hr','MarkerSize',8), end
    if hit(MDP(i)); plot(i,k,'hb','MarkerSize',8), end
    hold off
end


% run multiple subjects
%==========================================================================
BMR.g   = 1;
BMR.f   = 2;
BMR.T   = 3;
BMR.m   = @(i,i1,i2,i3,i4) i == i2;
OPT.BMR = BMR;

N     = 32;
Ns    = 64;
for m = 1:Ns
    
    % create structure array and solve
    %----------------------------------------------------------------------
    clear MDP OPTIONS
    for i = 1:N
        MDP(i) = mda;
    end
    rng(m)
    MDP  = spm_MDP_VB_X(MDP);
    
    % free energy and confidence
    %----------------------------------------------------------------------
    [F,Fu]   = spm_MDP_F(MDP);
    Fm(:,m)  = F(:);
    Fum(:,m) = Fu(:);
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(MDP(i)); h(i,m) = 1; else, h(i,m) = 0;  end
    end

    % repeat with BMR
    %----------------------------------------------------------------------
    RDP = MDP;
    RDP = rmfield(RDP,{'u','o'});
    for i = 1:N
        RDP(i).a  = mda.a;
        RDP(i).a0 = mda.a0;
        RDP(i).s  = RDP(i).s(:,1);
    end
    rng(m)
    RDP      = spm_MDP_VB_X(RDP,OPT);
    [F,Fu]   = spm_MDP_F(RDP);
    Rm(:,m)  = F(:);
    Rum(:,m) = Fu(:);
    
    % look for instances of BMR
    %----------------------------------------------------------------------
    vA    = spm_vec(A{1});
    for i = 1:N
       c(i,1) = corr(vA,(spm_vec(RDP(i).a{1}) > 0));
    end
    bmr(:,m) = diff(c);
    
    % preferred locations
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(RDP(i)); r(i,m) = 1; else, r(i,m) = 0;  end
    end
    
    % find run of correct responses
    %----------------------------------------------------------------------
    for i = 1:N
        if hit(RDP(i)); r(i,m) = 1; else, r(i,m) = 0;  end
    end
    
    % find non-learners
    %----------------------------------------------------------------------
    R = spm_conv(r,2,0);
    H = spm_conv(h,2,0);
    
    % show results
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 11');clf
    subplot(4,1,1)
    b = bar(mean(R(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.0),hold on
    b = bar(mean(H(:,:),2)); set(b,'EdgeColor','w','FaceColor',[1 1 1]*.8),hold off
    xlabel('trial'), ylabel('probability of correct'), axis([1/2 (N + 1/2) 1/3 1]);
    title('Average performance','Fontsize',16)
    
    subplot(4,1,2)
    spm_plot_ci(mean(Rm'),var(Rm')), hold on
    plot(mean(Fm'),'r'),              hold off
    xlabel('trial'), ylabel('free energy'); set(gca,'XLim',[1 N])
    title('Average free energy','Fontsize',16)
    
    subplot(4,1,3)
    spm_plot_ci(mean(Rum'),var(Rum')), hold on
    plot(mean(Fum'),'r'),              hold off
    xlabel('trial'), ylabel('confidence');  set(gca,'XLim',[1 N])
    title('Average confidence','Fontsize',16)
    
    % show individual performance
    %----------------------------------------------------------------------
    subplot(4,1,4), image(32*(r'))
    xlabel('trial'), ylabel('subject')
    title('Aha moments','Fontsize',16)
    
    % plot model updates
    %----------------------------------------------------------------------
    hold on
    for i = 1:m
        j = find(bmr(:,i) > 0) + 1;
        try, plot(j(1),i,'.m','MarkerSize',32), end
        try, plot(j(2),i,'.r','MarkerSize',32), end
        j = find(bmr(:,i) < 0);
        plot(j, (j - j + i),'.b','MarkerSize',32)
    end
    hold off, drawnow
    
    save paper
    
end


return


% confidence – negatively over policies
%--------------------------------------------------------------------------
for i = 1:numel(MDP)
    p     = MDP(i).R;
    Fu(i) = sum(sum(p.*log(p)));
end

return

function spm_MDP_plot_moves(MDP)
for i = 1:numel(MDP)
    m(i) = sum(~~diff(MDP(i).u(3,:)));
end
h = bar(m); set(h,'EdgeColor','w','FaceColor',[1 1 1]*.9), hold on

% and correct responses
%--------------------------------------------------------------------------
hit   = @(MDP) any(MDP.o(3,:) == 2) & ~any(MDP.o(3,:) == 3);
for i = 1:numel(MDP)
    if hit(MDP(i))
        h = plot(i,0,'hr','MarkerSize',8);
        set(h,'MarkerFaceColor',[1 1/2 0])
    end
end
return


function spm_MDP_A_plot(A)
% assemble key parts of the likelihood array
%==========================================================================
for i = 1:3
    for j = 1:3;
        a{i,j} = squeeze(A{1}(:,i,:,j,4));
        a{i,j} = a{i,j}*diag(1./sum(a{i,j}));
    end
end
a = spm_cat(a);
imagesc(a);
title( 'Sample: left - center - right', 'FontSize',16)
ylabel('Rule: left - center - right','FontSize',14)
xlabel('Correct color', 'FontSize',14)
set(gca,'XTick',1:9)
set(gca,'YTick',1:12)
set(gca,'XTicklabel',repmat(['r','g','b'],[1 3])')
set(gca,'YTicklabel',repmat(['r','g','b',' '],[1 3])')
axis image

return




function spm_MDP_rule_plot(MDP)
% illustrates visual search graphically
%==========================================================================

% locations
%--------------------------------------------------------------------------
x{1} = [-1  0; 0  1; 1  0; 0  0];
x{2} = [-1 -1; 0 -1; 1 -1; 0 -2]/2;
col  = {'r','g','b','c'};

% plot cues
%--------------------------------------------------------------------------
if strcmp('replace',get(gca,'Nextplot'))
    
    % plot cues
    %----------------------------------------------------------------------
    s     = MDP.s;hold off
    for i = 1:length(MDP.D{3})
        a = MDP.A{1}(:,s(1),s(2),i,1);
        j = find(rand < cumsum(a),1);
        plot(x{1}(i,1),x{1}(i,2),['.',col{j}],'MarkerSize',32), hold on
    end
    
    % plot choices
    %----------------------------------------------------------------------
    for i = 1:length(MDP.D{4})
        a = find(MDP.A{3}(:,s(1),s(2),4,i));
        if a == 2
            plot(x{2}(i,1),x{2}(i,2),['.m'],'MarkerSize',32), hold on
        end
        plot(x{2}(i,1),x{2}(i,2),['.',col{i}],'MarkerSize',16), hold on
        
    end
    axis([-2 2 -2 2]);
    
end

% Extract and plot eye movements and choice
%--------------------------------------------------------------------------
for i = 1:numel(MDP.o(2,:))
    X(i,:) = x{1}(MDP.o(2,i),:);
end
plot(X(:,1),X(:,2),'k')
for i = 1:numel(MDP.s(4,:))
    X(i,:) = x{2}(MDP.s(4,i),:);
end
plot(X(:,1),X(:,2),'k')
axis([-2 2 -2 2]);

return




