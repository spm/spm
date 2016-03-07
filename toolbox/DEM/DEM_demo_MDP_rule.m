function MDP = DEM_demo_MDP_rule
% Demo of active inference for visual salience
%__________________________________________________________________________
%
% This routine uses active inference for Markov decision processes to
% illustrate epistemic foraging in the context of visual searches. Here,
% the agent has to
%
% This demonstration uses a factorised version of the MDP scheme. In
% other words, we assume a mean field approximation to the posterior over
% different hidden states (context, location, etc) – and over
% multiple modalities (what versus where).  This provides a parsimonious
% representation of posterior beliefs over hidden states – but does induce
% degree of overconfidence associated with approximate Bayesian inference.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_rule.m 6741 2016-03-07 10:32:29Z karl $

% set up and preliminaries
%==========================================================================

% second level (semantic)
%==========================================================================
rng('default')

% prior beliefs about initial states (in terms of counts_: D and d
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
for f1 = 1:Ns(1)                     % location of trgaet colour
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
U(1,2,:)  = [1 1 2 4]';         % sample left
U(1,3,:)  = [1 1 3 4]';         % sample left
U(1,4,:)  = [1 1 4 1]';         % return and report red
U(1,5,:)  = [1 1 4 2]';         % return and report green
U(1,6,:)  = [1 1 4 3]';         % return and report blue


% priors: (utility) C; the agent expects to be right and not wrong
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end
C{3} = [ 0  0  0  0  0 -8 -8;
         0  0  0  0  0  0  0;
        -4 -4 -4 -4 -4 -4 -4];

% MDP Structure
%--------------------------------------------------------------------------
mdp.T = size(C{3},2);           % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.s = [1 1 4 4]';             % initial state

mdp.Aname = {'what','where','feedback'};
mdp.Bname = {'rule','colour','where','decision'};

mdp  = spm_MDP_check(mdp);


% illustrate a single trial
%==========================================================================
MDP  = spm_MDP_VB_X(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

subplot(3,2,3)
spm_MDP_rule_plot(MDP)

% return

% illustrate a sequence of trials
%==========================================================================
clear MDP

% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
N      = 32;
s(1,:) = ceil(rand(1,N)*3);
s(2,:) = ceil(rand(1,N)*3);
s(3,:) = 4;
s(4,:) = 4;

for i = 1:N
    MDP(i)   = mdp;      % create structure array
    MDP(i).s = s(:,i);   % context
end


% Solve - an example sequence
%==========================================================================
MDP  = spm_MDP_VB_X(MDP);

% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
    spm_MDP_VB_game(MDP);
Q = spm_MDP_VB_game(MDP);

% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP(1:8));

% illustrate behaviour in more detail
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf;
for i = 1:min(N,15)
    subplot(5,3,i), spm_MDP_rule_plot(MDP(i));
    axis square
end



% repeat with rule learning
%==========================================================================
for f1 = 1:Ns(1)                     % location of trgaet colour
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
mda      = mdp;
mda.a{1} = a{1};
mda.a{2} = mda.A{2}*128;
mda.a{3} = mda.A{3}*128;


% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
clear MDP
for i = 1:N
    MDP(i)   = mda;      % create structure array
    MDP(i).s = s(:,i);   % context
end


% Solve - an example sequence
%==========================================================================
MDP  = spm_MDP_VB_X(MDP);

spm_figure('GetWin','Figure 6'); clf;
for i = 1:min(N,15)
    subplot(5,3,i), spm_MDP_rule_plot(MDP(i + N - min(N,15)));
    axis square
end

spm_figure('GetWin','Figure 7'); clf;
spm_MDP_VB_trial(MDP(1));


spm_figure('GetWin','Figure 8'); clf;
spm_MDP_A_plot(MDP(end).a);


spm_figure('GetWin','Figure 9'); clf;
P = spm_MDP_VB_game(MDP);
subplot(2,1,1), plot(1:N,Q.p,1:N,P.p)
xlabel('trial'),ylabel('Expected utility')
title('Performance','Fontsize',16), set(gca,'YLim',[-8 -2]); axis square
legend({'optimal','learned'})

save

return

function spm_MDP_A_plot(A)
subplot(4,3,1),imagesc(squeeze(A{1}(:,1,:,1,4)));
ylabel('rule left');title('sample left'), axis image
subplot(4,3,2),imagesc(squeeze(A{1}(:,1,:,2,4)));
ylabel('rule left');title('sample centre'), axis image
subplot(4,3,3),imagesc(squeeze(A{1}(:,1,:,3,4)));
ylabel('rule left');title('sample right'), axis image
subplot(4,3,4),imagesc(squeeze(A{1}(:,2,:,1,4)));
ylabel('rule centre');title('sample left'), axis image
subplot(4,3,5),imagesc(squeeze(A{1}(:,2,:,2,4)));
ylabel('rule centre');title('sample centre'), axis image
subplot(4,3,6),imagesc(squeeze(A{1}(:,2,:,3,4)));
ylabel('rule centre');title('sample right'), axis image
subplot(4,3,7),imagesc(squeeze(A{1}(:,3,:,1,4)));
ylabel('rule right');title('sample left'), axis image
subplot(4,3,8),imagesc(squeeze(A{1}(:,3,:,2,4)));
ylabel('rule right');title('sample centre'), axis image
subplot(4,3,9),imagesc(squeeze(A{1}(:,3,:,3,4)));
ylabel('rule right');title('sample  right'), axis image


return

function spm_MDP_rule_plot(MDP)
% illustrates visual search graphically
%--------------------------------------------------------------------------

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

