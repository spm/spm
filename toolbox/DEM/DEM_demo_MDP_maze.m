function MDP = DEM_demo_MDP_maze
% Demo of active inference for epistemic foraging
%__________________________________________________________________________
%
% This routine uses the Markov decision process formulation of active
% inference (with variational Bayes) to model foraging for information in a
% three arm maze.  This demo illustrates variational free energy
% minimisation in the context of Markov decision processes, where the agent
% is equipped with prior beliefs that it will minimise expected free energy
% in the future. This free energy is the free energy of future sensory
% states expected under the posterior predictive distribution. It can be
% regarded as a generalisation of the variational formulation of KL control
% in which information gain or epistemic value is formulated explicitly.
%
% In this example, the agent starts at the centre of a three way maze
% which is baited with a reward in one of the two upper arms. However, the
% rewarded arm changes from trial to trial.  Crucially, the agent can
% identify where the reward (US) is located by accessing a cue (CS) in the
% lower arm. This tells the agent whether the reward is on the left or the
% right upper arm.  This means the optimal policy would first involve
% maximising information gain or epistemic value by moving to the lower arm
% and then claiming the reward this signified. Here, there are eight hidden
% states (four locations times right or left reward), four control states
% (that take the agent to the four locations) and 16 outcomes (four
% locations times two cues times two rewards).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically with an 80% schedule.
%
% A single trial is simulated followed by an examination of dopaminergic
% responses to conditioned and unconditioned stimuli (cues and rewards). A
% hierarchical version is then implemented, in which the mapping between
% locations in the generative model and the generative process is unknown
% and has to be learned.
%
% see also: spm_MPD_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_maze.m 6975 2016-12-18 20:27:00Z karl $

% set up and preliminaries
%==========================================================================
rng('default')

% observation probabilities
%--------------------------------------------------------------------------
a      = .9;
A{1,1} = [.5 .5; .5 .5; 0 0; 0 0];
A{2,2} = [0 0; 0 0; a (1 - a); (1 - a) a];
A{3,3} = [0 0; 0 0; (1 - a) a; a (1 - a)];
A{4,4} = [1 0; 0 1; 0 0; 0 0];
A      = spm_cat(A);

% transition probabilities
%--------------------------------------------------------------------------
for i = 1:4
    B{i} = zeros(4,4);
    B{i}([2 3],[2 3]) = eye(2);
    B{i}(i,[1 4])     = 1;
    B{i} = kron(B{i},eye(2));
end

% priors: softmax(utility)
%--------------------------------------------------------------------------
c  = 2;
C  = spm_softmax(kron(ones(4,1),[0; 0; c; -c]));

% prior beliefs about initial state
%--------------------------------------------------------------------------
D  = kron([1 0 0 0],[1 1]/2)';

% true initial state
%--------------------------------------------------------------------------
S  = kron([1 0 0 0],[1 0])';


% allowable policies (of depth T)
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  2  2  2  3  3  3  3  4  4  4  4
      1  2  3  4  1  2  3  4  1  2  3  4  1  2  3  4
      1  2  3  4  1  2  3  4  1  2  3  4  1  2  3  4];

 
% MDP Structure
%==========================================================================
MDP.N = 8;                          % number of variational iterations
MDP.S = S;                          % true initial state
MDP.A = A;                          % observation model
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % initial state probabilities (priors)
MDP.V = V;                          % allowable policies

MDP.alpha  = 64;                    % gamma hyperparameter
MDP.beta   = 4;                     % gamma hyperparameter
MDP.lambda = 1/4;                   % precision update rate


% Solve - an example game
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
MDP.plot = gcf;
MDP      = spm_MDP_game(MDP,'FE');


% different formulations of optimality as a function of preference
%==========================================================================
spm_figure('GetWin','Figure 2'); clf
MDP.plot = 0;
MDP.N    = 4;

n     = 128;                              % number of simulated trials
c     = linspace(0,1,6);
d     = kron(ones(4,1),[0; 0; 1; 0]);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    MDP.C = spm_softmax(kron(ones(4,1),[0; 0; c(i); -c(i)]));
    
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    for j = 1:n
        
        % randomise reward
        %------------------------------------------------------------------
        s       = rand > 1/2;
        MDP.S   = kron([1 0 0 0],[s ~s])';
        
        
        % inverst under different schemes
        %------------------------------------------------------------------
        MDP     = spm_MDP_game(MDP,'FE');
        FE(j,i) = d'*MDP.O(:,end);
        MDP     = spm_MDP_game(MDP,'KL');
        KL(j,i) = d'*MDP.O(:,end);
        MDP     = spm_MDP_game(MDP,'RL');
        RL(j,i) = d'*MDP.O(:,end);
        MDP     = spm_MDP_game(MDP,'FE',1);
        FP(j,i) = d'*MDP.O(:,end);
        try MDP = rmfield(MDP,'w'); end
        
    end
    
end
MDP.S  = S;
MDP.C  = C;

% plot behavioural results
%--------------------------------------------------------------------------
subplot(3,1,1)
bar(c,[mean(FE); mean(KL); mean(RL); mean(FP)]'*100), hold on
plot(c,c*0 + 100*3/8,'--k'), hold on
plot(c,c*0 + 100*a,  '-.k'), hold off
title('Performance','FontSize',16)
xlabel('Prior preference','FontSize',12)
ylabel('success rate (%)','FontSize',12)
spm_axis tight, axis square
legend({'FE','KL','RL','DA'})


% dopamine responses to US and CS
%==========================================================================
OPT   = 'FE';
MDP.N = 16;
MDP.a = [4 2 2];
MDP.o = [1 ((4 - 1)*4 + 1) ((2 - 1)*4 + 3)];
MDP   = spm_MDP_game(MDP,OPT);

% axis
%--------------------------------------------------------------------------
pst   = (1:length(MDP.d))*100/1000;
ax    = [1 pst(end) 0 4];

subplot(3,2,3)
plot(pst,MDP.d,'k'), hold on
subplot(3,2,5)
plot(pst,MDP.da,'k'), hold on


subplot(3,2,6)
r     = 128;
bar(pst,r*MDP.da + randn(size(MDP.da)).*sqrt(r*MDP.da))
title('Simulated (CS & US) responses','FontSize',16)
xlabel('Peristimulus time (sec)','FontSize',12)
ylabel('Rate','FontSize',12)
axis square, set(gca,'XLim',ax(1:2))


% repeat but with the US only
%--------------------------------------------------------------------------
MDP.a = [1 1 2];
MDP.o = [1 ((1 - 1)*4 + 1) ((2 - 1)*4 + 3)];
MDP   = spm_MDP_game(MDP,OPT);

subplot(3,2,3)
plot(pst,MDP.d,'r'), hold off
title('Precision updates','FontSize',16)
xlabel('Peristimulus time (sec)','FontSize',12)
ylabel('Precision','FontSize',12)
axis square, set(gca,'XLim',ax(1:2))

subplot(3,2,5)
plot(pst,MDP.da,'r'), hold off
title('Dopamine responses','FontSize',16)
xlabel('Peristimulus time (sec)','FontSize',12)
ylabel('Response','FontSize',12)
axis square, axis(ax)

subplot(3,2,4)
r     = 128;
bar(pst,r*MDP.da + randn(size(MDP.da)).*sqrt(r*MDP.da))
title('Simulated (US) responses','FontSize',16)
xlabel('Peristimulus time (sec)','FontSize',12)
ylabel('Rate','FontSize',12)
axis square, set(gca,'XLim',ax(1:2))


% different levels of priors and uncertainty
%==========================================================================
spm_figure('GetWin','Figure 3'); clf; 

MDP.a = [4 2 2];
MDP.o = [1 ((4 - 1)*4 + 1) ((2 - 1)*4 + 3)];
n     = 3;
c     = linspace(0,2,n);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    MDP.C = spm_softmax(kron(ones(4,1),[0; 0; c(i); -c(i)]));
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    MDP  = spm_MDP_game(MDP,'FE');
    col  = [1 1 1]*(length(c) - i)/length(c);
    
    subplot(2,2,1)
    plot(pst,MDP.da,'Color',col), hold on
    title('Preference (utility)','FontSize',16)
    xlabel('Peristimulus time (sec)','FontSize',12)
    ylabel('Response','FontSize',12)
    axis square, axis(ax)
    
    subplot(2*n,2,(i - 1)*2 + 2)
    r     = 128;
    bar(pst,r*MDP.da + randn(size(MDP.da)).*sqrt(r*MDP.da))
    ylabel('Rate','FontSize',12)
    axis square, set(gca,'XLim',ax(1:2))
    
end



% and uncertainty
%--------------------------------------------------------------------------
subplot(2,1,2)

MDP.C = C;
c     = linspace(.5,.9,n);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    a      = c(i);
    U{1,1} = [.5 .5; .5 .5; 0 0; 0 0];
    U{2,2} = [0 0; 0 0; a (1 - a); (1 - a) a];
    U{3,3} = [0 0; 0 0; (1 - a) a; a (1 - a)];
    U{4,4} = [1 0; 0 1; 0 0; 0 0];
    MDP.A  = spm_cat(U);
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    MDP  = spm_MDP_game(MDP,'FE');
    col  = [1 1 1]*(length(c) - i)/length(c);
    
    subplot(2,2,3)
    plot(pst,MDP.da,'Color',col), hold on
    title('Uncertainty','FontSize',16)
    xlabel('Peristimulus time (sec)','FontSize',12)
    ylabel('Response','FontSize',12)
    axis square, axis(ax)
    
    subplot(2*n,2,(i - 1)*2 + 2 + n*2)
    r     = 128;
    bar(pst,r*MDP.da + randn(size(MDP.da)).*sqrt(r*MDP.da))
    ylabel('Rate','FontSize',12)
    axis square, set(gca,'XLim',ax(1:2))
    
end
xlabel('Peristimulus time','FontSize',12)

% learniing the hidden context
%==========================================================================

% mapping from hidden lcoations to the (unkown) world
%--------------------------------------------------------------------------
M     = [1 2 3 4; 1 3 2 4; 1 4 3 2; 1 2 4 3]';

% A
%--------------------------------------------------------------------------
m     = size(M,1);
MDP.A = kron(ones(1,m),A);

% B
%--------------------------------------------------------------------------
for i = 1:length(B)
    MDP.B{i} = spm_cat(spm_diag(B(M(i,:))));
end

% C
%--------------------------------------------------------------------------
MDP.C  = C;


% simulate trials and record outcomes
%==========================================================================
spm_figure('GetWin','Figure 4'); clf
MDP.plot = 0;


MDP.N = 8;
MDP.a = [];
MDP.o = [];
RDP   = MDP;

d     = kron(ones(4,1),[0; 0; 1; 0]);
DD    = kron(eye(m,m) + 1/8,D*ones(1,length(D)));
DD    = DD*diag(1./sum(DD));
for j = 1:128
    
    % initial context
    %----------------------------------------------------------------------
    MDP.D  = kron(ones(m,1)/m,D);
    RDP.D  = kron(ones(m,1)/m,D);
    ss     = find(rand < cumsum(ones(1,m)/m),1);
    ss     = sparse(ss,1,1,m,1);
    
    for i = 1:8
        
        % initial state
        %------------------------------------------------------------------
        s      = rand > 1/2;
        s      = kron([1 0 0 0],[s ~s])';
        MDP.S  = kron(ss,s);
        RDP.S  = kron(ss,s);
        
        % trial
        %------------------------------------------------------------------
        MDP    = spm_MDP_game(MDP,'FE');
        RDP    = spm_MDP_game(MDP,'RL');
        
        R(j,i) = d'*MDP.O(:,end);
        H(j,i) = -log(MDP.Q(:,end))'*MDP.Q(:,end);
        G{i,j} = MDP.da;
        
        r(j,i) = d'*RDP.O(:,end);
        h(j,i) = -log(RDP.Q(:,end))'*RDP.Q(:,end);
        
        % Bayesian update
        %------------------------------------------------------------------
        MDP.D  = spm_softmax(log(DD*MDP.Q(:,end)));
        RDP.D  = spm_softmax(log(DD*RDP.Q(:,end)));
        MDP.gamma = MDP.d(end);
        RDP.gamma = RDP.d(end);
        
    end

end

subplot(2,1,1)
t  = 1:size(R,2);
plot(t,100*mean(R),'-k',t,100*mean(R),'ok'), hold on
plot(t,100*mean(H),'-r',t,100*mean(H),'or')
plot(t,100*mean(r),':k',t,100*mean(r),'ok')
plot(t,100*mean(h),':r',t,100*mean(h),'or')
plot(t,t*0 + 100*3/8,':k'), hold on
plot(t,t*0 + 100*a,  '--k'), hold off
title('Exploration and exploitation','FontSize',16)
xlabel('Number of trials','FontSize',12)
ylabel('Performance and uncertainty','FontSize',12)
axis([1 t(end) 0 100])

% plot simulated dopamine responses
%--------------------------------------------------------------------------
r  = 128;
G  = spm_cat(G);
G  = mean(G,2);

subplot(4,1,3)
plot(G)
title('Average dopaminergic response','FontSize',16)
xlabel('Variational updates','FontSize',12)
ylabel('Precision','FontSize',12)
spm_axis tight

subplot(4,1,4)
bar(r*G + randn(size(G)).*sqrt(r*G))
title('Simulated dopaminergic response','FontSize',16)
xlabel('Variational updates','FontSize',12)
ylabel('Spikes per then','FontSize',12)
spm_axis tight








