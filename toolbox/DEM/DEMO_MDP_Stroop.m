function [mdp, Z] = DEMO_MDP_Stroop(p)
% This demo uses a deep temporal partially observed Markov decision
% process to simulate performance of a Stroop task. This task is used to
% illustrate a formulation of cognitive or mental effort. The synthetic
% participants must overcome a prior belief that the normal action on
% being presented with text is to read that word, and instead state the
% colour the text is presented in. In addition, this routine demonstrates
% the fitting of choice and reaction time data to the model, and the
% recovery of parameters that summarise behaviour.
% 
% see also: spm_MDP_VB_X.m, DEM_demo_MDP_reading.m, DEMO_MDP_questions.m
%__________________________________________________________________________

% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%==========================================================================
rng('default')  % For reproducibility

try e = exp(p.e)*0.85;    catch, e = 0.85;  end % Bias towards reading word
try c = exp(p.c);         catch, c = 1;     end % Preference for being correct

W = 1;          % Name colour [1] or read word [2]
N = 16;         % Number of words

% First level POMDP (HMM in this case)
%==========================================================================

label.factor = {'Written word','Word Colour','Task Sequence','Instruction','Response','Correct?'};

% Priors over initial states P(s{j} = i, t = 1) = D{j}(i)
%--------------------------------------------------------------------------
D{1} = ones(4,1); % Written word:  Red, green, blue, yellow
D{2} = ones(4,1); % Word colour:   Red, green, blue, yellow
D{3} = ones(3,1); % Task sequence: Instruction, null, response
D{4} = ones(2,1); % Instruction:   Colour, written
D{5} = ones(2,1); % Response:      Colour, written
D{6} = ones(2,1); % Correct?:      Correct, incorrect

% Transition probabilities P(s{j} = i, t+1|s{j} = k, u = m, t) = B{j}(i,k,m)
%--------------------------------------------------------------------------
B{1} = eye(4);    % Written word does not change during trial
B{2} = eye(4);    % Word colour does not change during trial
B{3} = [1 0 0;    % Stay in instruction phase if started there
        0 0 0;    % Move to response if in null phase
        0 1 1];   % Stay in response phase if there already
B{4} = eye(2);    % Instructions do not change throughout trial
B{5} = eye(2);    % Response does not change throughout trial
B{5}(:,:,2) = eye(2); % Cheat to get around issue in HMM code
B{6} = eye(2);    % Correctness of response does not change

% Likelihood (P(o{m} = i|s{n1} = j1, s{n2} = j2...) = A{m}(i,j1,j2,...)
%--------------------------------------------------------------------------
for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            for f4 = 1:length(D{4})
                for f5 = 1:length(D{5})
                    for f6 = 1:length(D{6})
                        % Written/colour of word, and instruction:
                        %--------------------------------------------------
                        if f3 == 1
                            A{1}(5,f1,f2,f3,f4,f5,f6)  = 1; % No word during instruction
                            A{2}(5,f1,f2,f3,f4,f5,f6)  = 1;
                            A{3}(f4,f1,f2,f3,f4,f5,f6) = 1; % Instruction given
                        else
                            A{1}(f1,f1,f2,f3,f4,f5,f6) = 1; % Otherwise, depends on written word factor
                            A{2}(f2,f1,f2,f3,f4,f5,f6) = 1; % or word colour factor
                            A{3}(3,f1,f2,f3,f4,f5,f6)  = 1; % No instruction given
                        end
                        
                        % Verbal response
                        %--------------------------------------------------
                        if f3 == 3
                            if f5 == 1
                                A{4}(f2,f1,f2,f3,f4,f5,f6) = 1; % Expect to report colour
                            else
                                A{4}(f1,f1,f2,f3,f4,f5,f6) = 1; % Expect to report written word
                            end
                        else
                            A{4}(5,f1,f2,f3,f4,f5,f6) = 1;      % Expect null word
                        end
                    end
                end
            end
        end
    end
end

% Preferences log(P(o{j} = i|C)) = C{j}(i) + const.
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = zeros(size(A{g},1),1);
end

% Self-generated outcomes (verbal response)
%--------------------------------------------------------------------------
n      = zeros(numel(A),2);
n(4,:) = 1;

% Compile POMDP
%--------------------------------------------------------------------------
MDP.T       = 2;
MDP.A       = A;
MDP.B       = B;
MDP.C       = C;
MDP.D       = D;
MDP.n       = n;
MDP.tau     = 8;
MDP.chi     = -exp(64);
MDP.label   = label;
MDP.lambda  = 1/4;
MDP         = spm_MDP_check(MDP);

clear A B C D n label

% Second level POMDP
%==========================================================================

label.factor = {'Narrative','Instruction','Response'};

% Priors over initial states P(s{j} = i, t = 1) = D{j}(i)
%--------------------------------------------------------------------------
D{1} = [1 0]'; % Narrative: instruction, response
D{2} = [1 1]'; % Instruction: colour, written
D{3} = [1 1]'; % Response: colour, written

% Transition probabilities P(s{j} = i, t+1|s{j} = k, u = m, t) = B{j}(i,k,m)
%--------------------------------------------------------------------------
B{1} = [0 0; 1 1]; % Once instruction has been given, responses
B{2} = eye(2);     % Instruction is constant over time
B{3} = zeros(2,2,2);
B{3}(1,:,2) = 1;   % Choose to respond with colour of word
B{3}(2,:,1) = 1;   % Choose to respond by reading written word

% Likelihood probabilities P(s{j} = i, t+1|s{j} = k, u = m, t) = B{j}(i,k,m)
%--------------------------------------------------------------------------
for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            % Predicted task sequence (factor 3 lower level)
            %--------------------------------------------------------------
            if f1 == 1
                A{1}(1,f1,f2,f3) = 1;
            else
                A{1}(2,f1,f2,f3) = 1;
            end
            A{1}(3,f1,f2,f3) = 0;
            
            % Predicted instruction (factor 4 at lower level)
            %--------------------------------------------------------------
            A{2}(f2,f1,f2,f3) = 1;
            
            % Predicted response (factor 5 at lower level)
            %--------------------------------------------------------------
            A{3}(f3,f1,f2,f3) = 1;
            
            % Predicted correct? (factor 6 at lower level)
            %--------------------------------------------------------------
            if f2 == f3
                A{4}(1,f1,f2,f3) = 1;
            else
                A{4}(2,f1,f2,f3) = 1;
            end
        end
    end
end

% Preferences log(P(o{j} = i|C)) = C{j}(i) + const.
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = zeros(size(A{g},1),1);
end
C{4} = [c;-c]; % Prefers to be correct

% Policies
%--------------------------------------------------------------------------
E    = spm_softmax([e;-e]);

% Compile POMDP structure
%--------------------------------------------------------------------------
mdp.MDP = MDP; clear MDP

mdp.A = A;
mdp.B = B;
mdp.C = C;
mdp.D = D;
mdp.E = E;
mdp.T = N;
mdp.s = [1;W;1];
mdp.link = zeros(numel(mdp.MDP.D),numel(mdp.A));
mdp.link(3,1) = 1;
mdp.link(4,2) = 1;
mdp.link(5,3) = 1;
mdp.link(6,4) = 1;
mdp.label     = label;
mdp.tau       = 8;
mdp           = spm_MDP_check(mdp);
OPTIONS.gamma = 1;

if exist('p'), return, end

% Solve POMDP
%--------------------------------------------------------------------------
MDP = spm_MDP_VB_X(mdp,OPTIONS);

% Synthetic electrophysiology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
MDP_Stroop_beliefs(MDP);

% Plot example sequence and beliefs
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
MDP_Stroop_animation(MDP);

% Congruency x task condition (ERPs, errors and reaction times)
%--------------------------------------------------------------------------
mdp.T = 64;
mdp.s = [1;1;1];
MDP1  = spm_MDP_VB_X(mdp);

% Plot effort
%--------------------------------------------------------------------------
% spm_figure('GetWin','Cognitive effort'); clf
% subplot(1,2,1)
% MDP_Stroop_effort(MDP1)

mdp.s = [1;2;1];
MDP2  = spm_MDP_VB_X(mdp);

% subplot(1,2,2)
% MDP_Stroop_effort(MDP2)

spm_figure('GetWin','Figure 3'); clf
[U, ~, ~] = MDP_Stroop_beliefs(MDP1); clf

% Determine behavioural measures from stimuli and responses
%--------------------------------------------------------------------------
[stim, resp] = MDP_Stroop_SR(MDP1);
H1 = zeros(3,numel(resp));
for i = 1:numel(resp)
    H1(1,i) = strcmp(stim.color{i},stim.word{i});      % Congruent?
    H1(2,i) = strcmp(['"' stim.color{i} '"'],resp{i}); % Correct?
end
H1(3,:) = MDP_Stroop_RT(MDP1);                         % Reaction time

[stim, resp] = MDP_Stroop_SR(MDP2);
H2 = zeros(3,numel(resp));
for i = 1:numel(resp)
    H2(1,i) = strcmp(stim.color{i},stim.word{i});      % Congruent?
    H2(2,i) = strcmp(['"' stim.word{i} '"'],resp{i});  % Correct?
end
H2(3,:) = MDP_Stroop_RT(MDP2);                         % Reaction time


% Plot performance summaries in each condition
%--------------------------------------------------------------------------
subplot(3,1,1)
B = [sum(H1(2,H1(1,:)==1))/sum(H1(1,:));
     sum(H1(2,H1(1,:)==0))/sum(1-H1(1,:));
     sum(H2(2,H1(1,:)==1))/sum(H1(1,:));
     sum(H2(2,H1(1,:)==0))/sum(1-H1(1,:))];
X = categorical({'Colour (cong.)','Colour (incong.)','Word (cong.)', 'Word (incong.)'});
X = reordercats(X,{'Colour (cong.)','Colour (incong.)','Word (cong.)', 'Word (incong.)'});
bar(X,B*100,'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
axis square
title('Percentage correct')

% Plot reaction times in each condition
%--------------------------------------------------------------------------
subplot(3,1,2)
histogram(exp(H1(3,H1(1,:)==1)+randn(size(H1(3,H1(1,:)==1)))/16)/2,0.4:0.02:1.2), hold on
histogram(exp(H1(3,H1(1,:)==0)+randn(size(H1(3,H1(1,:)==0)))/16)/2,0.4:0.02:1.2), hold on
histogram(exp(H2(3,H2(1,:)==1)+randn(size(H2(3,H2(1,:)==1)))/16)/2,0.4:0.02:1.2), hold on
histogram(exp(H2(3,H2(1,:)==0)+randn(size(H2(3,H2(1,:)==0)))/16)/2,0.4:0.02:1.2)
axis square
xlabel('Reaction Time (s)')
ylabel('Counts')
legend('Colour (cong.)','Colour (incong.)','Word (cong.)', 'Word (incong.)')
title('Reaction time distribution')

% Plot condition-specific electophysiological responses
%--------------------------------------------------------------------------
subplot(3,1,3)
n1 = 0;
n2 = 0;
u1 = zeros(17,1);
u2 = zeros(17,1);
for i = 1:length(H1(2,:))
    if H1(2,i) == 1
        if H1(1,i) == 1
            u1 = u1 + U{2}(i*16-8:(i+1)*16-8,3)+U{2}(i*16-8:(i+1)*16-8,5); % LFPs associated with response factor at second level
            n1 = n1 + 1;
        else
            u2 = u2 + U{2}(i*16-8:(i+1)*16-8,3)+U{2}(i*16-8:(i+1)*16-8,5); % LFPs associated with response factor at second level
            n2 = n2 + 1;
        end
    end
end
u1 = u1/n1; u1 = u1 - u1(1);
u2 = u2/n2; u2 = u2 - u2(1);
plot((0:16)*500/16, u1), hold on, plot((0:16)*500/16, u2)
ax = gca;
ax.XAxisLocation = 'origin';
axis square
box off
legend('Congruent','Incongruent')
ylabel('a.u.')
xlabel('Time (ms)')
title('Evoked responses')

% Relationship between errors, reaction times, and priors
%--------------------------------------------------------------------------
clear MDP X resp stim

[x, y] = meshgrid((-1:1/2:1)/4,(-1:1/2:1)/4);
X = [x(:),y(:)]; % Design matrix - columns are log(c), log(e)

% Generate data under alternative priors:
for i = 1:size(X,1)
    p.c   = X(i,1);
    p.e   = X(i,2);
    mdp    = DEMO_MDP_Stroop(p);
    mdp.T  = 32;
    rng default % Standardise stimulus stream
    mdp.s(2) = 1;
    MDP(i,1) = spm_MDP_VB_X(mdp,OPTIONS);                                  
    mdp.s(2) = 2;                                                          
    MDP(i,2) = spm_MDP_VB_X(mdp,OPTIONS);                                  
end

% Performance summaries:
B = zeros(1,size(X,1));
R = zeros(2*(MDP(1).T-1),size(X,1));
for k = 1:size(MDP,1)
    [stim{k}, resp{k}] = MDP_Stroop_SR(MDP(k,1));
    H = zeros(2,numel(resp{k}));
    for i = 1:numel(resp{k})
        H(1,i) = strcmp(stim{k}.color{i},stim{k}.word{i});       % Congruent?
        H(2,i) = strcmp(['"' stim{k}.color{i} '"'],resp{k}{i});  % Correct?
    end
    B(k)     = sum(H(2,:))/length(H(2,:));
    R(1:MDP(1).T-1,k)   = MDP_Stroop_RT(MDP(k,1));               % Reaction time 
    R(1:MDP(1).T-1,k)   = exp(R(1:MDP(1).T-1,k) + randn(size(R(1:MDP(1).T-1,k)))/16)/2; 

    [stim{k}, resp{k}] = MDP_Stroop_SR(MDP(k,2));                          
    H = zeros(2,numel(resp{k}));
    for i = 1:numel(resp{k})
        H(1,i) = strcmp(stim{k}.color{i},stim{k}.word{i});       % Congruent?
        H(2,i) = strcmp(['"' stim{k}.word{i} '"'],resp{k}{i});   % Correct? 
    end
    B(k)     = (B(k) + sum(H(2,:))/length(H(2,:)))/2;
    R(MDP(1).T:end,k)   = MDP_Stroop_RT(MDP(k,2));               % Reaction time
    R(MDP(1).T:end,k)   = exp(R(MDP(1).T:end,k) + randn(size(R(MDP(1).T:end,k)))/16)/2;
end

% Fit linear model to synthetic data
%--------------------------------------------------------------------------
U = [ones(size(X,1),1), X, X(:,1).*X(:,2), X(:,1).^2, X(:,2).^2];
M.L       = @(P,M,U,Y) sum(log( spm_Npdf(Y, exp(U*P.beta)./(exp(U*P.beta)+exp(1)), exp(-P.pi)) ));
M.pE.beta = zeros(size(U,2),1);                   % prior means (parameters)
M.pE.pi   = 4;
M.pC      = eye(spm_length(M.pE));                % prior variance (parameters)
[Ep,Cp,F] = spm_nlsi_Newton(M,U,B');
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.F  = F;
DCM.M  = M;
[PCM,~,BMA] = spm_dcm_bmr_all(DCM,'beta');
Bp  = BMA.Ep.beta;
BPp = PCM.Pp.beta;

M.L       = @(P,M,U,Y) sum(log( spm_Npdf(log(Y), U*P.beta, exp(-P.pi)) ));
[Ep,Cp,F] = spm_nlsi_Newton(M,U,mean(R)');
DCM.Ep = Ep;
DCM.Cp = Cp;
DCM.F  = F;
DCM.M  = M;
[PCM,~,BMA] = spm_dcm_bmr_all(DCM,'beta');
Rp  = BMA.Ep.beta;
RPp = PCM.Pp.beta;

spm_figure('GetWin','Figure 4'); clf
subplot(3,2,1)
imagesc(X(:,1),X(:,2),100*reshape(B,[sqrt(length(B)),sqrt(length(B))]))
title('Percentage correct')
colorbar
axis square, axis xy
xlabel('c')
ylabel('e')

subplot(3,2,2)
imagesc(X(:,1),X(:,2),reshape(mean(R),[sqrt(length(mean(R))),sqrt(length(mean(R)))]))
title('Mean reaction time')
colorbar
axis square, axis xy
xlabel('c')
ylabel('e')

subplot(3,2,3)
imagesc(X(:,1),X(:,2),100*reshape(exp(U*Bp)./(exp(U*Bp)+exp(1)),[sqrt(length(B)),sqrt(length(B))]))
title('Percentage correct (fit)')
colorbar
axis square, axis xy
xlabel('c')
ylabel('e')
caxis([min(B*100),max(B*100)])

subplot(3,2,4)
imagesc(X(:,1),X(:,2),reshape(exp(U*Rp),[sqrt(size(R,2)),sqrt(size(R,2))]))
title('Mean reaction time (fit)')
colorbar
axis square, axis xy
xlabel('c')
ylabel('e')
caxis([min(mean(R)),max(mean(R))])

subplot(3,2,5)
bar(BPp,'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
axis square
title('Posterior probabilities (% correct)')
xlabel('Parameter')
ylabel('Probability')

subplot(3,2,6)
bar(RPp,'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
axis square
title('Posterior probabilities (reaction time)')
xlabel('Parameter')
ylabel('Probability')


return

% Model fits
%--------------------------------------------------------------------------
% This section demonstrates the fitting of this model to synthetic data,
% and shows that it is possible to recover (combinations of) the prior
% belief parameters that determine behaviour.

clear M
U  = [];
Ep = [];
Cp = [];
M.L       = @MDP_Stroop_L;
M.G       = @(p) MDP_Stroop_Gen(p, MDP(1));
M.pE.c    = 0;                         % prior means (parameters)
M.pE.e    = 0;
M.pC      = eye(spm_length(M.pE))/256; % prior variance (parameters)
M.ch      = 1; % Include choice data
M.rt      = 1; % Include reaction time data

for i = 1:size(MDP,1)
    Y.o = {};
    for j = 1:length(MDP(i,1).mdp)                                         
        Y.o{j} = [MDP(i,1).mdp(j).o MDP(i,2).mdp(j).o];                    
    end
    Y.r = R(:,i);
    [EP,CP,~] = spm_nlsi_Newton(M,U,Y);
    display(['Inverted model ' num2str(i) '/' num2str(length(MDP))])
    Ep(i,:) = spm_vec(EP)';
    Cp(i,:) = [diag(CP)' CP(1,2)];
    spm_figure('GetWin','Figure 5'); clf
    subplot(3,1,1)
    bar(1:size(Ep,1),Ep(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Ep,1),Ep(:,1),1.65*sqrt(Cp(:,1)),-1.65*sqrt(Cp(:,1)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Ep,1),X(1:i,1),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Demand')
    subplot(3,1,2)
    bar(1:size(Ep,1),Ep(:,2),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Ep,1),Ep(:,2),1.65*sqrt(Cp(:,2)),-1.65*sqrt(Cp(:,2)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Ep,1),X(1:i,2),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Effort')
    subplot(3,1,3)
    bar(1:size(Ep,1),Ep(:,2)-Ep(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
    hold on
    errorbar(1:size(Ep,1),Ep(:,2)-Ep(:,1),1.65*sqrt(Cp(:,1)+Cp(:,2)-2*Cp(:,3)),-1.65*sqrt(Cp(:,1)+Cp(:,2)-2*Cp(:,3)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Ep,1),X(1:i,2)-X(1:i,1),'.r','MarkerSize',20)
    xlabel('Dataset')
    ylabel('Parameters')
    title('Differences')
    hold off
end

Z.Ep = Ep;
Z.Cp = Cp;
save('Z','Z')

% Repeat the above with each data-source omitted (to caculate information
% gained from each source).
%--------------------------------------------------------------------------
M.ch      = 0; % Omit choice data
for i = 1:size(MDP,1)                                                      
    Y.o = {};
    for j = 1:length(MDP(i,1).mdp)                                         
        Y.o{j} = [MDP(i,1).mdp(j).o MDP(i,2).mdp(j).o];                    
    end
    Y.r = R(:,i);
    [EP,CP,~] = spm_nlsi_Newton(M,U,Y);
    display(['Inverted model ' num2str(i) '/' num2str(length(MDP))])
    Eprt(i,:) = spm_vec(EP)';
    Cprt(i,:) = [diag(CP)' CP(1,2)];
    spm_figure('GetWin','Figure 6'); clf
    subplot(3,1,1)
    bar(1:size(Eprt,1),Eprt(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Eprt,1),Eprt(:,1),1.65*sqrt(Cprt(:,1)),-1.65*sqrt(Cprt(:,1)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Eprt,1),X(1:i,1),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Demand')
    subplot(3,1,2)
    bar(1:size(Eprt,1),Eprt(:,2),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Eprt,1),Eprt(:,2),1.65*sqrt(Cprt(:,2)),-1.65*sqrt(Cprt(:,2)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Eprt,1),X(1:i,2),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Effort')
    subplot(3,1,3)
    bar(1:size(Eprt,1),Eprt(:,2)-Eprt(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
    hold on
    errorbar(1:size(Eprt,1),Eprt(:,2)-Eprt(:,1),1.65*sqrt(Cprt(:,1)+Cprt(:,2)-2*Cprt(:,3)),-1.65*sqrt(Cprt(:,1)+Cprt(:,2)-2*Cprt(:,3)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Eprt,1),X(1:i,2)-X(1:i,1),'.r','MarkerSize',20)
    xlabel('Dataset')
    ylabel('Parameters')
    title('Differences')
    hold off
end

Z.Eprt = Eprt;
Z.Cprt = Cprt;
save('Z','Z')

M.ch      = 1; % Include choice data
M.rt      = 0; % Omit reaction time data

for i = 1:size(MDP,1)                                                      
    Y.o = {};
    for j = 1:length(MDP(i,1).mdp)                                         
        Y.o{j} = [MDP(i,1).mdp(j).o MDP(i,2).mdp(j).o];                    
    end
    Y.r = R(:,i);
    [EP,CP,~] = spm_nlsi_Newton(M,U,Y);
    display(['Inverted model ' num2str(i) '/' num2str(length(MDP))])
    Epch(i,:) = spm_vec(EP)';
    Cpch(i,:) = [diag(CP)' CP(1,2)];
    spm_figure('GetWin','Figure 6'); clf
    subplot(3,1,1)
    bar(1:size(Epch,1),Epch(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Epch,1),Epch(:,1),1.65*sqrt(Cpch(:,1)),-1.65*sqrt(Cpch(:,1)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Epch,1),X(1:i,1),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Demand')
    subplot(3,1,2)
    bar(1:size(Epch,1),Epch(:,2),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1]), hold on, errorbar(1:size(Epch,1),Epch(:,2),1.65*sqrt(Cpch(:,2)),-1.65*sqrt(Cpch(:,2)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Epch,1),X(1:i,2),'.r','MarkerSize',20)
    hold off
    xlabel('Dataset')
    ylabel('Parameter estimate')
    title('Effort')
    subplot(3,1,3)
    bar(1:size(Epch,1),Epch(:,2)-Epch(:,1),'FaceColor',[.7 .7 .9],'EdgeColor',[1 1 1])
    hold on
    errorbar(1:size(Epch,1),Epch(:,2)-Epch(:,1),1.65*sqrt(Cpch(:,1)+Cpch(:,2)-2*Cpch(:,3)),-1.65*sqrt(Cpch(:,1)+Cpch(:,2)-2*Cpch(:,3)),'CapSize',0,'LineWidth',2,'LineStyle','None','Color',[.3 .3 .6])
    plot(1:size(Epch,1),X(1:i,2)-X(1:i,1),'.r','MarkerSize',20)
    xlabel('Dataset')
    ylabel('Parameters')
    title('Differences')
    hold off
end

Z.Epch = Epch;
Z.Cpch = Cpch;
save('Z','Z')

% Compute information gain in each condition
%--------------------------------------------------------------------------
IG = zeros(25,3);
for i = 1:25
    IG(i,1) = spm_kl_normal(Ep(i,:),[Cp(i,1),Cp(i,3);Cp(i,3),Cp(i,2)],spm_vec(M.pE),M.pC);
    IG(i,2) = spm_kl_normal(Eprt(i,:),[Cprt(i,1),Cprt(i,3);Cprt(i,3),Cprt(i,2)],spm_vec(M.pE),M.pC);
    IG(i,3) = spm_kl_normal(Epch(i,:),[Cpch(i,1),Cpch(i,3);Cpch(i,3),Cpch(i,2)],spm_vec(M.pE),M.pC);
end

Z.IG = IG;
save('Z','Z')

function [stim, resp] = MDP_Stroop_animation(MDP)
% This routine returns the stimuli and response sequence from a solved 
% POMDP problem, along with an animation of that sequence.

str = {'red','green','blue','yellow',' '};
col = {[1,0,0],[0 1 0],[0 0 1],[1 1 0]};

if isfield(MDP,'mdp')
    stim.word  = [];
    stim.color = [];
    resp       = [];
    for i = 1:length(MDP.mdp)
        subplot(6,2,1), cla, axis off, title('Stimulus')
        subplot(6,2,2), cla, axis off, title('Response')
        pause(1/64)
        [s, r] = MDP_Stroop_animation(MDP.mdp(i));
        stim.word  = [stim.word, s.word];
        stim.color = [stim.color, s.color];
        resp       = [resp, r];
        if i>1
            subplot(6,2,3)
            text(0.5,0.5-i/4, stim.word(i-1),'Color',stim.color{i-1},'HorizontalAlignment','center','units','normalized','FontSize',10), hold on
            axis off
            subplot(6,2,4)
            text(0.5,0.5-i/4, resp(i-1),'Color','k','HorizontalAlignment','center','units','normalized','FontSize',10), hold on
            axis off
        end
    end
else
    stim.color = [];
    stim.word  = [];
    resp       = [];
    for t = 1:MDP.T
        subplot(6,2,1), cla
        if MDP.o(3,t) == 3
            text(0.5,0.5,str(MDP.o(1,t)),'Color',col{MDP.o(2,t)},'HorizontalAlignment','center','units','normalized','FontSize',12);
            if t == 2, stim.color = [stim.color, {col{MDP.o(2,t)}}]; stim.word = [stim.word, str(MDP.o(1,t))]; end
        else
            if MDP.o(3,t) == 1
                text(0.5,0.5,'Give the colours of the font of the following words','HorizontalAlignment','center','units','normalized','FontSize',10)
            else
                text(0.5,0.5,'Read the words that follow','HorizontalAlignment','center','units','normalized','FontSize',10)
            end
        end
        axis off
        title('Stimulus')
        subplot(6,2,2), cla
        if MDP.o(4,t) ~= 5
            text(0.5,0.5,['"' str{MDP.o(4,t)} '"'],'HorizontalAlignment','center','units','normalized','FontSize',12)
        end
        if t==2 && ~isempty(stim.word), resp  = [resp, {['"' str{MDP.o(4,t)} '"']}]; end
        axis off
        title('Response')
        pause(1/4)
    end
end

function [U, V, P] = MDP_Stroop_beliefs(MDP)
% This function reports the beliefs of a synthetic participant, presented
% as if they were electrophyisological data.

% Preliminaries
%--------------------------------------------------------------------------
r   = 16;
s   = 0.2;
h   = 64;                % Suppression of high frequencies
T   = size(MDP.mdp(1).xn{1},3);
XX  = [];
uu  = [];
C   = zeros(2,MDP.T*T);  % Congruency (first row), Correct (second row)

% Beliefs at second level
%--------------------------------------------------------------------------
for i = 1:numel(MDP.xn)
    X = [];
    u = [];
    for j = 1:size(MDP.xn{i},3)
        X = [X;MDP.xn{i}(:,:,j,j)];
        u = [u;gradient(MDP.xn{i}(:,:,j,j)')'];
    end
    XX = [XX X];
    uu = [uu u];
end

% Convert (second level) belief states into synthetic electrophysiology
%--------------------------------------------------------------------------
P  = MDP.R;
XX = kron(XX,ones(round(r*T),r));
R  = rand(size(XX));
XX(R>XX-s) = 0;
XX(R<XX+s) = 1;
duu   = spm_dctmtx(size(uu,1))*uu;
k     = repmat(exp((1-(1:size(duu,1)))/h)',1,size(duu,2));
uu    = spm_dctmtx(size(uu,1))'*(k.*duu);
U{1} = XX;
U{2} = uu;
U{3} = 250+(1:size(uu,1))*250*T/16; % Time

P = kron(P,ones(r,round(r*T*size(MDP.xn{1},1))));
R = rand(size(P));
P(R>P-s) = 0;
P(R<P+s) = 1;

% Plotting of second level belief states
%--------------------------------------------------------------------------
subplot(3,2,3)
imagesc(1 - XX'), colormap gray
title('Beliefs - second level')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('States')

subplot(3,2,4)
plot(U{3},uu)
xlim([0 size(uu,1)*250*T/16])
xlabel('Time (ms)')
ylabel('a.u.')
title('Local field potentials')

% First level beliefs
%--------------------------------------------------------------------------
XXX = [];
uuu = [];
for t = 1:length(MDP.mdp)
    XX = [];
    uu = [];
    for i = 1:numel(MDP.mdp(t).xn)
        X = [];
        u = [];
        for j = 1:size(MDP.mdp(t).xn{i},3)
            if i == 1
                C(1,(t-1)*T + j) = MDP.mdp(t).o(1,j) == MDP.mdp(t).o(2,j);
                if MDP.mdp(t).o(4,j) == 5, C(2,(t-1)*T + j) = 0.5;
                elseif MDP.s(2,1) == 1
                    C(2,(t-1)*T + j) = MDP.mdp(t).o(2,j) == MDP.mdp(t).o(4,j);
                else
                    C(2,(t-1)*T + j) = MDP.mdp(t).o(1,j) == MDP.mdp(t).o(4,j);
                end
            end
            X = [X;MDP.mdp(t).xn{i}(:,:,j,j)];
            u = [u;gradient(MDP.mdp(t).xn{i}(:,:,j,j)')'];
        end
        XX = [XX X];
        uu = [uu u];
    end
    XXX = [XXX; XX];
    uuu = [uuu; uu];
end

% Convert first level beliefs into synthetic electrophysiology
%--------------------------------------------------------------------------
XXX = kron(XXX,ones(r,r));
R   = rand(size(XXX));
XXX(R>XXX-s) = 0;
XXX(R<XXX+s) = 1;
V{1} = XXX;
duu  = spm_dctmtx(size(uuu,1))*uuu;
k    = repmat(exp((1-(1:size(duu,1)))/h)',1,size(duu,2));
k    = k./sum(k,2);
uuu  = spm_dctmtx(size(uuu,1))'*(k.*duu);
V{2} = uuu;
V{3} = (1:size(uuu,1))*250/16;

% Plot first level beliefs
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(1 - P), colormap gray
title('Beliefs - policy')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('Policies')

subplot(3,2,5)
imagesc(1 - XXX'), colormap gray
title('Beliefs - first level')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('States')

subplot(3,2,6)
plot(V{3},uuu)
xlabel('Time (ms)')
ylabel('a.u.')
title('Local field potentials')
subplot(3,2,2)
imagesc(C)
set(gca,'YTick',1:2)
set(gca,'YTickLabel',{'Congruency','Correctness'})
xlabel('Time steps (discrete)')

function rt = MDP_Stroop_RT(MDP)
% Simulate reaction times (based upon predictive entropy)

for i = 2:length(MDP.mdp)
    x = [];
    for k = 1:size(MDP.mdp(i).xn{1},1)
        for j = 1:numel(MDP.mdp(i).xn)
            xn{j} = MDP.mdp(i).xn{j}(k,:,2,1);
        end
        x(:,end+1) = spm_dot(MDP.mdp(i).A{4},xn);
    end
    v       = -diag(x'*spm_log(x));
    rt(i-1) = find(v<8/10,1);
    rt(i-1) = v(end);
end

function [stim, resp] = MDP_Stroop_SR(MDP)
% This routine reports simuli and response sequences without generating an
% animation (see above)

str = {'red','green','blue','yellow',' '};
col = {'red','green','blue','yellow'};

if isfield(MDP,'mdp')
    stim.word  = [];
    stim.color = [];
    resp       = [];
    for i = 1:length(MDP.mdp)
        [s, r] = MDP_Stroop_SR(MDP.mdp(i));
        stim.word  = [stim.word, s.word];
        stim.color = [stim.color, s.color];
        resp       = [resp, r];
    end
else
    stim.color = [];
    stim.word  = [];
    resp       = [];
    for t = 1:MDP.T
        if MDP.o(3,t) == 3
            if t == 2, stim.color = [stim.color, {col{MDP.o(2,t)}}]; stim.word = [stim.word, str(MDP.o(1,t))]; end
        end
        if t==2 && ~isempty(stim.word), resp  = [resp, {['"' str{MDP.o(4,t)} '"']}]; end
    end
end


function L = MDP_Stroop_L(P,M,~,Y)
% Likelihood function for model fitting.

L = 0;

mdp(1)   = M.G(P);
mdp(2)   = M.G(P); 

% Assign outcomes
%--------------------------------------------------------------------------
for i = 1:mdp(1).T                                                         
    mdp(1).mdp(i).o = Y.o{i}(:,1:2);  
    mdp(2).mdp(i).o = Y.o{i}(:,3:4);    
end

% Invert model (forcing it to take same actions)
%--------------------------------------------------------------------------
OPTIONS.gamma = 1;
MDP(1) = spm_MDP_VB_X(mdp(1),OPTIONS);                                     
MDP(2) = spm_MDP_VB_X(mdp(2),OPTIONS); 

if M.ch == 1
    
% Evaluate likelihood of actions
%--------------------------------------------------------------------------
    for i = 2:length(MDP(1).mdp)
        x1 = [];
        x2 = [];
        for k = 1:size(MDP(1).mdp(i).xn{1},1)                              
            for j = 1:numel(MDP(1).mdp(i).xn)                              
                xn{j} = MDP(1).mdp(i).xn{j}(end,:,2,1);                    
            end
            x1(:,end+1) = spm_softmax(MDP(1).mdp(i).lambda*spm_log(spm_dot(MDP(1).mdp(i).A{4},xn)));

            for j = 1:numel(MDP(2).mdp(i).xn)                              
                xn{j} = MDP(2).mdp(i).xn{j}(end,:,2,1);                    
            end
            x2(:,end+1) = spm_softmax(MDP(2).mdp(i).lambda*spm_log(spm_dot(MDP(2).mdp(i).A{4},xn))); 

            L = L + spm_log(x1(MDP(1).mdp(i).o(4,2),end)) + spm_log(x2(MDP(2).mdp(i).o(4,2),end));
        end
    end
    
end

if M.rt == 1
% Evaluate likelihood of reaction times
%--------------------------------------------------------------------------
    R1 = MDP_Stroop_RT(MDP(1));
    R2 = MDP_Stroop_RT(MDP(2));
    R  = [R1 R2];
    L = L + sum(log( spm_Npdf(log(Y.r), R' - log(2), 1/256 )));
end

function MDP = MDP_Stroop_Gen(p, MDP)
% POMDP for generative model
e = exp(p.e)*0.85;
c = exp(p.c);
MDP.E    = spm_softmax([e;-e]);
MDP.C{4} = [c;-c];


function MDP_Stroop_effort(MDP)
% Plot effort over time
p  = spm_softmax(log(MDP.E) + MDP.G);
xi = sum(p.*(MDP.G) - log(MDP.E'*exp(MDP.G)));
barh(xi)

return

% Illustrate concept of cognitive effort
%------------------------------------------------------------------------
o{1} = [0.9;0.1];  % Policy 1 leads to outcome 1 with probability 0.9
o{2} = [0.1;0.9];  % Policy 2 leads to outcome 2 with probability 0.9
effort = zeros(100,100);
prob   = zeros(100,100);
X      = zeros(100,100);
Y      = zeros(100,100);
for i  = 1:100
  for j = 1:100
      C = spm_softmax([i;-i]/10);
      E = spm_softmax([-j;j]/10);
      G = [o{1}'*(log(o{1}) - log(C));o{2}'*(log(o{2}) - log(C))]; 
      p = spm_softmax(log(E) - G);
      effort(i,j) = p'*(-G) - log(E'*exp(-G));
      prob(i,j) = p(1);
      X(i,j) = E(2);
      Y(i,j) = C(1);
  end
end
spm_figure('GetWin','Effort'); clf
subplot(2,1,1), imagesc(Y(:,1),X(1,:),effort'), axis square, axis xy, title('Effort'), ylabel('Cognitive demand (E)'), xlabel('Preference (C)'), colorbar
subplot(2,1,2), imagesc(Y(:,1),X(1,:),prob'), axis square, axis xy, title('Probability of overcoming cognitive demand'), ylabel('Cognitive demand (E)'), xlabel('Preference (C)'), colorbar
