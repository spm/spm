function MDP = DEM_reasoning
% Demo of active inference for visual salience
%__________________________________________________________________________
% This routine illustrates active reasoning; namely, choosing policies that
% disambiguate among prior beliefs about hypotheses or model structure. In
% this example, the hypotheses concern the likelihood mapping between
% latent causes and a particular outcome modality (correct or not). The
% example used to demonstrate this form of active model selection revisits
% the three-ball paradigm, in which the correct colour is the colour of a
% ball in a location signified by the colour of another ball. In other
% words, this is a minimal setup or puzzle that features an unknown rule,
% invariance, equivariance or symmetry. In this case, a rule is taken to
% mean the equivalence between a control and stayed that depends upon a
% context state. Here, the context is specified by the colour of a ball in
% one of three locations. This colour specifies the location of the
% criterion colour that has to be matched by the agent's choice.
%
% In this example, we assume that the agent has learned the generative
% model of the game but does not know the likelihood mapping between latent
% states and feedback or reward that tells us whether it's choice was
% correct or incorrect. However, the agent is equipped with priors over the
% likelihood mapping in question, encoding all possible (second-order)
% rules of the kind described above. This enables it to evaluate the
% likelihood of each prior model using Bayesian model reduction (starting
% with prior based upon the Bayesian model average of prior models). These
% marginal likelihoods can then be used to evaluate the expected
% information gain over models, conditioned on each policy. The policies in
% question here depend upon where the agent is looking and the choice it
% commits to.
%
% The expected information gain supplements the salience (information gain
% about states) and novelty (information gain about parameters) terms in
% expected free energy; thereby causing the agent to sample information
% that best resolves uncertainty about which model is in the play. When
% sufficiently confident (as scored by its beliefs about models) it can
% then replace its posterior likelihood mapping with the model selected
% from its model priors and, if correct, proceed as usual to secure as many
% wins as possible.
%
% In what follows, we first demonstrate the optimal behaviour of an agent
% that knows the correct rule, and plans two steps into the future under
% sophisticated inference. We then repeat the simulation removing prior
% preferences (and deep planning) to let the expected information gain over
% models predominate in action selection. When the log odds ratio of the
% most likely model being the only viable model exceeds 16, preferences are
% reinstated and the agent proceeds in an exploitative fashion. The ensuing
% behaviour is characterised in terms of the agents confidence; i.e., the
% expected precision or negentropy of beliefs over policies, performance
% and rule discovery (as scored by the KL divergence between the true
% likelihood mapping and the posterior encoded in Dirichlet counts
% accumulated over successive trials.
%
% see also: spm_MDP_BMR.m and spm_MPD_VB_XXX.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%--------------------------------------------------------------------------
rng(1)

% latent states and outcomes
%==========================================================================
label.factor   = {'left', 'centre', 'right', 'where', 'choice'};
label.name     = { ...
    {'red','green','blue'}, ...
    {'red','green','blue'}, ...
    {'red','green','blue'}, ...
    {'left','centre','right','start'}, ...
    {'red','green','blue','?'}};

label.action   = { ...
    {' '}, ...
    {' '}, ...
    {' '}, ...
    {'left','centre','right','start'}, ...
    {'red','green','blue','?'}};


label.modality = {'what', 'where', 'feedback'};
label.outcome  = { ...
    {'red','green','blue','null'}, ...
    {'left','centre','right','start'}, ...
    {'null','right','wrong'}};

% initialise likelihoods A and priors B
%==========================================================================
Nu = [1,1,1,4,4];
Nf = numel(label.factor);
Ng = numel(label.modality);
Ns = zeros(1,Nf);
for f = 1:Nf
    Ns(f) = numel(label.name{f});
end
No = zeros(1,Ng);
for g = 1:Ng
    No(g) = numel(label.outcome{g});
end
for g = 1:Ng
    A{g}    = zeros([No(g),Ns]);
    C{g}    = zeros([No(g),Ns]);
    id.A{g} = 1:Nf;
    id.C{g} = 1:Nf;
end
for f = 1:Nf
    B{f} = zeros(Ns(f),Ns(f),Nu(f));
    D{f} = zeros(Ns(f),1);
end

% for all combinations of latent states
%--------------------------------------------------------------------------
S     = spm_combinations(Ns);
for i = 1:size(S,1)

    % probabilistic mapping from hidden states to outcomes: A
    %----------------------------------------------------------------------
    f1 = S(i,1);                        % colour on left
    f2 = S(i,2);                        % colour at centre
    f3 = S(i,3);                        % colour on right
    f4 = S(i,4);                        % location of fixation
    f5 = S(i,5);                        % choice

    % A{1} what: {'red','green','blue','start'}
    %======================================================================
    if f4 == 1, A{1}(f1,f1,f2,f3,f4,f5) = 1;  end
    if f4 == 2, A{1}(f2,f1,f2,f3,f4,f5) = 1;  end
    if f4 == 3, A{1}(f3,f1,f2,f3,f4,f5) = 1;  end
    if f4 == 4, A{1}( 4,f1,f2,f3,f4,f5) = 1;  end

    % A{2} where: {'left','centre','right','null'}
    %----------------------------------------------------------------------
    A{2}(f4,f1,f2,f3,f4,f5) = 1;

    % A{3} feedback: {'null','right','wrong'}
    %----------------------------------------------------------------------
    if f5 == 4
        A{3}(1,f1,f2,f3,f4,f5) = 1;     % undecided
    else

        if f2 == 1, col = f1; end
        if f2 == 2, col = f2; end
        if f2 == 3, col = f3; end

        if f5 == col
            A{3}(2,f1,f2,f3,f4,f5) = 1; % right
        else
            A{3}(3,f1,f2,f3,f4,f5) = 1; % wrong
        end

    end
end

% alternative mappings - hypothesis or model space
%==========================================================================
pA{3} = spm_MDP_rules(A{3});

% model average and priors
%--------------------------------------------------------------------------
A0    = 0;
for m = 1:numel(pA{3})
    A0 = A0 + pA{3}{m};
end
A0    = A0/sum(A0,'all');

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% control states B(4): where {'left','centre','right','null'}
%--------------------------------------------------------------------------
for k = 1:Ns(4)
    B{4}(:,:,k) = 0;
    B{4}(k,:,k) = 1;
end

% control states B(5): report {'red','green','blue','?'}
%--------------------------------------------------------------------------
for k = 1:Ns(5)
    B{5}(:,:,k) = 0;
    B{5}(k,:,k) = 1;
end

% allowable policies (specified as the next action) U
%--------------------------------------------------------------------------
U    = [0,0,0,1,1];

% priors: (utility) C: the agent expects to avoid mistakes
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end

% and expects itself to be right
%--------------------------------------------------------------------------
C{3} = spm_softmax([0;2;-6]);
C0   = spm_softmax([0;0;-0]);

% prior beliefs about initial states (in terms of counts_: D
%--------------------------------------------------------------------------
D{1} = [1 1 1]';        % what1:  {'red','green','blue'}
D{2} = [1 1 1]';        % what2:  {'red','green','blue'}
D{3} = [1 1 1]';        % what3:  {'red','green','blue'}

D{4} = [0 0 0 1]';      % where:  {'left','centre','right','start'}
D{5} = [0 0 0 1]';      % report: {'red','green','blue','?'}

% MDP Structure
%--------------------------------------------------------------------------
mdp.T = 6;              % number of moves
mdp.U = U;              % allowable policies
mdp.A = A;              % observation model
mdp.B = B;              % transition probabilities
mdp.C = C;              % preferred outcomes
mdp.D = D;              % prior over initial states
mdp.N = 1;              % depth

mdp.label = label;

% illustrate a single trial
%==========================================================================
MDP  = spm_MDP_VB_XXX(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

% illustrate behaviour
%--------------------------------------------------------------------------
subplot(3,2,3)
spm_MDP_rule_plot(MDP)


% repeat with active reasoning
%==========================================================================

% retain knowledge about the cues but remove knowledge about the rule
%--------------------------------------------------------------------------
MDP      = mdp;
MDP.a{1} = A{1}*512;
MDP.a{2} = A{2}*512;
MDP.a{3} = A0;
MDP.C{3} = C0;
MDP.pA   = pA;
MDP.N    = 0;

% show trial-by-trial evidence accumulation
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf;
N     = 32;                                    % number of trials
RA    = MDP.A{3};                              % likelihood rule
for i = 1:N

    % solve
    %----------------------------------------------------------------------
    PDP = spm_MDP_VB_XXX(MDP);

    % update
    %----------------------------------------------------------------------
    MDP = spm_MDP_VB_update(MDP,PDP);

    % update priors if Occam's window contains a winner
    %======================================================================

    % posterior probabilities over rules
    %----------------------------------------------------------------------
    P(i,:) = PDP.Pa{3};
    L      = max(PDP.Pa{3});
    O(i)   = log(L) - log(1 - L);
    if O(i) > 16

        % Baysian model average
        %------------------------------------------------------------------
        sp   = 0;
        for j = 1:numel(PDP.Pa{3})
            sp = sp + PDP.Pa{3}(j)*pA{3}{j};
        end
        MDP.a{3} = sp*512;
        MDP.C{3} = C{3};
        MDP.N    = 1;

    end

    % record
    %----------------------------------------------------------------------
    F(i) = sum(PDP.F);                         % ELBO (states)
    w(i) = sum(PDP.w);                         % policy precision
    o(i) = sum(PDP.o(3,:) == 2);               % wins
    p(i) = sum(PDP.o(3,:) == 3);               % losses
    K(i) = spm_KL_cat(RA,PDP.a{3});            % KL divergence from rule

    % slices of inferred and true likelihood tensor
    %----------------------------------------------------------------------
    a       = spm_dir_norm(PDP.a{3}(:,:,:,:,1,:));
    a3(i,:) = a(2,:);

    a       = spm_dir_norm(PDP.A{3}(:,:,:,:,1,:));
    A3(i,:) = a(2,:);

    % report
    %----------------------------------------------------------------------
    clf
    subplot(6,1,1), plot(w),  xlabel('trial'), spm_axis tight, title('Confidence', 'Fontsize',16), ylabel('nats')
    subplot(6,1,2), plot(o),  xlabel('trial'), spm_axis tight, title('Outcomes',   'Fontsize',16), ylabel('score'), hold on
    subplot(6,1,2), plot(p),  xlabel('trial'), spm_axis tight, title('Outcomes',   'Fontsize',16), ylabel('score'), hold off
    subplot(6,1,3), plot(K),  xlabel('trial'), spm_axis tight, title('Divergence',         'Fontsize',16), ylabel('nats')
    subplot(4,2,5), imagesc(a3),       title('Posterior likelihoods','Fontsize',16), ylabel('trial'), xlabel('latent state')
    subplot(4,2,7), imagesc(A3),       title('True likelihoods',     'Fontsize',16), ylabel('trial'), xlabel('latent state')
    subplot(4,2,6), imagesc(P.^(1/4)), title('Model comparison',     'Fontsize',16), ylabel('trial'), xlabel('rule or model')
    subplot(4,2,8), plot(O),           title('Occams Razor',         'Fontsize',16), xlabel('trial'), ylabel('nats')
    drawnow

end

return

% Comparative solutions 
%==========================================================================
spm_figure('GetWin','Figure 3'); clf;
N     = 64;                                    % number of trials
for k = 1:64                                   % for K agents

    % initalise MDP for this schemes
    %----------------------------------------------------------------------
    MDP      = mdp;
    MDP.a{1} = A{1}*512;
    MDP.a{2} = A{2}*512;
    MDP.a{3} = A0;
    MDP.C{3} = C0;
    MDP.N    = 0;
    MDP.pA   = pA;

    for i = 1:N                                % for N epochs

        % solve
        %------------------------------------------------------------------
        PDP = spm_MDP_VB_XXX(MDP);

        % update
        %------------------------------------------------------------------
        MDP = spm_MDP_VB_update(MDP,PDP);

        % update priors if Occam's window contains a winner
        %==================================================================
        L     = max(PDP.Pa{3});
        O(i)  = log(L) - log(1 - L);
        if O(i) > 16

            % Baysian model average
            %--------------------------------------------------------------
            sp   = 0;
            for j = 1:numel(PDP.Pa{3})
                sp = sp + PDP.Pa{3}(j)*pA{3}{j};
            end
            MDP.a{3} = sp*512;
            MDP.C{3} = C{3};
            MDP.N    = 1;

        end

        % record
        %------------------------------------------------------------------
        o(i) = sum(PDP.o(3,:) == 2);           % wins
        o(i) = o(i) - sum(PDP.o(3,:) == 3);    % minus losses
        v(i) = sum(log(PDP.C{3}(PDP.o(3,:)))); % expected value
        K(i) = spm_KL_cat(RA,PDP.a{3});        % divergence from rule

        pause(1)

    end

    % report
    %----------------------------------------------------------------------
    subplot(4,1,3), plot(cumsum(o),'k'), hold on
    xlabel('trial'), ylabel('score'), axis([0,N,-40,300])
    title('Cumulative score','Fontsize',16)

    subplot(4,1,1), plot(K,'k'), hold on
    xlabel('trial'), ylabel('nats'),  axis([0,N,0,200])
    title('KL Divergence','Fontsize',16)

    subplot(4,1,2), plot(O,'k'), hold on
    xlabel('trial'), ylabel('nats'),  axis([0,N,-5,25])
    title('Occams Razor', 'Fontsize',16)
    drawnow

    % performance
    %----------------------------------------------------------------------
    WP(k)  = sum(o);                          % cumulative wins
    KT(k)  = sum(K > 1);                      % discovery time

end

% summary
%--------------------------------------------------------------------------
subplot(4,2,7), histogram(WP,-40:20:300)
xlabel('score'), ylabel('frequency'), title('Performance','Fontsize',16)
subplot(4,2,8), histogram(KT,0:4:N)
xlabel('trial'), ylabel('frequency'), title('Discovery time','Fontsize',16)

% cd('C:\Users\Karl\Dropbox\matlab')
% save Reasoning.mat


return

% alternative mappings - hypothesis or model space
%==========================================================================
function pa = spm_MDP_rules(A)
% returns (likelihood) hypothesis or model space
% FORMAT pa = spm_MDP_rules(A)
%__________________________________________________________________________


% preliminaries
%--------------------------------------------------------------------------
Ns    = size(A);
Ns(1) = [];
S     = spm_combinations(Ns);
R     = spm_combinations([3,3,3,3]);
pa    = cell(1,size(R,1));
for j = 1:size(R,1)

    % for this rule
    %----------------------------------------------------------------------
    a     = zeros(size(A));
    for i = 1:size(S,1)

        % probabilistic mapping from hidden states to outcomes: A
        %------------------------------------------------------------------
        s1 = S(i,1);     % colour on left
        s2 = S(i,2);     % colour at centre
        s3 = S(i,3);     % colour on right
        s4 = S(i,4);     % location of fixation
        s5 = S(i,5);     % control

        % A{3} feedback: {'null','right','wrong'}
        %------------------------------------------------------------------
        if s5 == 4
            a(1,s1,s2,s3,s4,s5) = 1;     % undecided
        else
            for k = 1:Ns(R(j,1))
                if S(i,R(j,1)) == k
                    criterion = S(i,R(j,k + 1));
                end
            end

            if s5 == criterion
                a(2,s1,s2,s3,s4,s5) = 1; % right
            else
                a(3,s1,s2,s3,s4,s5) = 1; % wrong
            end

        end
    end

    % append this model or rule
    %----------------------------------------------------------------------
    pa{j} = a;

end


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
    for i = 1:3
        plot(x{1}(i,1),x{1}(i,2),['.',col{s(i)}],'MarkerSize',32), hold on
    end

    % plot choices
    %----------------------------------------------------------------------
    for i = 1:4
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
for i = 1:numel(MDP.s(5,:))
    Y(i,:) = x{2}(MDP.s(5,i),:);
end
plot(Y(:,1),Y(:,2),'k')
axis([-2 2 -2 2]);

return

% NOTES for paper graphics
%__________________________________________________________________________

% get results of various solutions
%--------------------------------------------------------------------------
cd('C:\Users\Karl\Dropbox\matlab')

load Reasoning1
W1 = WP;
T1 = KT;
load Reasoning2
W2 = WP;
T2 = KT;
load Reasoning3
W3 = WP;
T3 = KT;
load Reasoning4
W4 = WP;
T4 = KT;

disp('EIG over parameter and states')
disp('Mean time to discovery:'),disp(mean(T1))
disp('Mean correct-incorrect:'),disp(mean(W1))

disp('EIG over states')
disp('Mean time to discovery:'),disp(mean(T2))
disp('Mean correct-incorrect:'),disp(mean(W2))

disp('No BMS')
disp('Mean time to discovery:'),disp(mean(T3))
disp('Mean correct-incorrect:'),disp(mean(W3))

disp('No EIG')
disp('Mean time to discovery:'),disp(mean(T4))
disp('Mean correct-incorrect:'),disp(mean(W4))

% summary
%==========================================================================
spm_figure('GetWin','Figure 4'); clf;

% compare solutions with and without EIG over parameters
%--------------------------------------------------------------------------
axw = [-64,256];
axt = [1,64];
subplot(3,2,1), plot(T1,T2,'ok',axt,axt,'r'), axis([axt,axt]), axis square
xlabel('trials'), ylabel('trials: no EIG over parameters'), title('Discovery time','Fontsize',16)
subplot(3,2,2), plot(W1,W2,'ok',axw,axw,'r'), axis([axw,axw]), axis square
xlabel('score'),  ylabel('score: no EIG over parameters'),  title('Performance','Fontsize',16)

% compare solutions with and without EIG over parameters and states
%--------------------------------------------------------------------------
subplot(3,2,3), plot(T1,T4,'ok',axt,axt,'r'), axis([axt,axt]), axis square
xlabel('trials'), ylabel('trials: : no EIG over parameters or states'), title('Discovery time','Fontsize',16)
subplot(3,2,4), plot(W1,W4,'ok',axw,axw,'r'), axis([axw,axw]), axis square
xlabel('score'),  ylabel('score: no EIG over parameters or states'),  title('Performance','Fontsize',16)

% compare solutions with and without EIG over parameters and states
%--------------------------------------------------------------------------
subplot(3,1,3), plot(W1,W3,'ok',axw,axw,'r'), axis([axw,axw]), axis square
xlabel('score'),  ylabel('score: no model selection'),  title('Performance','Fontsize',16)







