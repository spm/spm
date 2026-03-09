function MDP = DEM_thermostat
% Demo of multiple constraint satisfaction (smart thermostat)
%__________________________________________________________________________
% 
% This demonstration routine illustrates the application of active
% inference to a multiple constraint satisfaction problem. In this example,
% there are two controllers (i.e., agents) the first can boost the charge
% of a domestic battery. The second is a heating, ventilation and air
% conditioning system (HVAC). These two controllers have to operate in
% tandem to satisfy multiple constraints; namely, keep electricity costs
% down by boosting during periods of low time of use (ToU) rates, while
% maintaining room temperature at a comfortable level, when and only when
% the rule is in use (i.e., occupied). The problem is solved using
% constraints (i.e., prior preferences) on outcomes. These include
% consumption (i.e., time of use rates multiplied by the state of the
% booster), room temperature and battery charge (state of charge; which
% cannot be at the highest or lowest level of charge).
%
% This example illustrates the discretisation or course graining of latent
% states and outcomes at several timescales. First, it illustrates home to
% use the solutions of underlying state space model (e.g., differential
% equations) to specify outcomes. Second, it uses a model to illustrate how
% slow fluctuations (day-to-day) contextualise fast fluctuations (hour two
% hour) in hidden states generating outcomes. For example, the time of use
% rate depends upon a daily average, with two hourly peaks. Similarly, the
% outside (air) temperature fluctuates relatively smoothly, with an initial
% (early morning) temperature that changes itinerantly from day to day.%
%
% The results of active inference are shown using a bespoke displays scheme
% that picks out the outcomes that are equipped with constraints and
% evaluates constraint satisfaction in terms of the cost; namely, the
% negative log probability of the realised outcome under prior preferences
% encoded in (context sensitive) constraints (C). Context sensitivity in
% this instance is endowed by conditioning prior constraints on one or more
% hidden factors (here, constraints over room temperature depend upon
% whether the room is in use or not).
%
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging
 
% set up and preliminaries
%==========================================================================
rng('default')

%--------------------------------------------------------------------------
% First, we will specify the POMDP in terms of factors and dependencies.
% This involves specifying the nodes of a factor graph that comprise the
% likelihood and transition priors at the first and second levels,
% respectively. The specification rests upon specifying the requisite (ABC)
% tensors and the edges specifying dependencies (iA,iC,iD,iE): see
% spm_parents.
%--------------------------------------------------------------------------
  
% latent causes: likelihoods A and contraints C
%==========================================================================
n = 8;                                % number of temperature levels
m = 5;                                % number of HVAC levels
N = 1;                                % depth of planning
D = 7;                                % number of days
U = [1,1];                            % enable charge and HVAC


% ToU Rate
%--------------------------------------------------------------------------
A = zeros(22,2,4,12);
s = size(A);
c = spm_combinations(s(2:end));
for i = 1:size(c,1)
    s = c(i,:);

    %   Charge x Daily_Rate x hourly fluctuations
    %----------------------------------------------------------------------
    charge     = s(1);
    hour       = (s(3) - 1)/12;
    Daily_Rate = s(2) + 4*(1 + sin(2*pi*2*hour));

    o  = fix(charge*Daily_Rate);
    A(o,s(1),s(2),s(3)) = 1;

end
level(1).outcomes(1).name  = 'ToU_Rate';
level(1).outcomes(1).state = {};
level(1).outcomes(1).A     = A;
level(1).outcomes(1).iA    = {'Charge','Daily_Rate','Hour'};
level(1).outcomes(1).C     = spm_softmax(-(1:size(A,1))');
level(1).outcomes(1).iC    = {};


% State of charge
%--------------------------------------------------------------------------
A = zeros(m + 2,2,2,12,m);
s = size(A);
c = spm_combinations(s(2:end));
for i = 1:size(c,1)
    s = c(i,:);

    %  Charge + Solar - HVAC: Solar = Cloud cover x daylight hours
    %----------------------------------------------------------------------
    charge = s(1);
    cloud  = (s(2) - 1);
    hour   = (s(3) - 1)/12;
    solar  = cloud*(1 + cos(2*pi*hour));
    hvac   = abs(s(4) - 2);

    o      = fix(charge + solar - hvac) + m - 2;
    A(o,s(1),s(2),s(3),s(4)) = 1;

end
level(1).outcomes(2).name = 'SoC';
level(1).outcomes(2).state = {};
level(1).outcomes(2).A    = A;
level(1).outcomes(2).iA   = {'Charge','Cover','Hour','HVAC'};
level(1).outcomes(2).C    = spm_softmax(sparse(2:(m + 1),1,32,m + 2,1));
level(1).outcomes(2).iC   = {};

% Solar power
%--------------------------------------------------------------------------
A = zeros(9,2,12);
s = size(A);
c = spm_combinations(s(2:end));
for i = 1:size(c,1)
    s = c(i,:);

    % Solar = Cloud cover x daylight hours
    %----------------------------------------------------------------------
    cloud  = (s(1) - 1);
    hour   = (s(2) - 1)/12;
    solar  = cloud*(1 + cos(2*pi*hour));

    o      = fix(4*solar);
    A(o + 1,s(1),s(2)) = 1;

end
level(1).outcomes(3).name  = 'Solar';
level(1).outcomes(3).state = {};
level(1).outcomes(3).A     = A;
level(1).outcomes(3).iA    = {'Cover','Hour'};
level(1).outcomes(3).C     = [];
level(1).outcomes(3).iC    = {};

% Room_Temperature temperature
%--------------------------------------------------------------------------
A = zeros(m + n - 1,m,n);
s = size(A);
c = spm_combinations(s(2:end));
for i = 1:size(c,1)
    s = c(i,:);

    % HVAC - (Room temperature - (outside) Temperature)
    %----------------------------------------------------------------------
    o = m - s(1) + s(2);
    o = fix(o);
    A(o,s(1),s(2)) = 1;
end
level(1).outcomes(4).name   = 'Room_Temperature';
level(1).outcomes(4).state  = {};
level(1).outcomes(4).A      = A;
level(1).outcomes(4).iA     = {'HVAC','Temperature'};
level(1).outcomes(4).C(:,1) = spm_Npdf((1:(m + n - 1))',4,8);
level(1).outcomes(4).C(:,2) = spm_Npdf((1:(m + n - 1))',4,1);
level(1).outcomes(4).iC     = {'Occupied'};

% Air temperature
%--------------------------------------------------------------------------
level(1).outcomes(5).name  = 'Temperature';
level(1).outcomes(5).state = {};
level(1).outcomes(5).A     = eye(n,n);
level(1).outcomes(5).iA    = {'Temperature'};
level(1).outcomes(5).C     = [];
level(1).outcomes(5).iC    = {};

% Motion detection
%--------------------------------------------------------------------------
level(1).outcomes(6).name  = 'Motion';
level(1).outcomes(5).state = {'no','yes'};
level(1).outcomes(6).A     = eye(2,2);
level(1).outcomes(6).iA    = {'Occupied'};
level(1).outcomes(6).C     = [];
level(1).outcomes(6).iC    = {};


% latent states: prior transitions B and contraints H
%==========================================================================

% Clock time
%--------------------------------------------------------------------------
level(1).states(1).name  = 'Hour';
level(1).states(1).state = num2cell(2:2:24);
level(1).states(1).B     = spm_speye(12,12,-1,1);
level(1).states(1).iD    = {'Hour'};
level(1).states(1).iE    = {};
level(1).states(1).H     = [];
level(1).states(1).U     = 0;


% Charge
%--------------------------------------------------------------------------
level(1).states(2).name     = 'Charge';
level(1).states(2).state    = {'off','on'};
level(1).states(2).B(:,:,1) = kron([1;0],[1 1]);
level(1).states(2).B(:,:,2) = kron([0;1],[1 1]);
level(1).states(2).iD   = {};
level(1).states(2).iE   = {};
level(1).states(2).H    = [];
level(1).states(2).U    = U(1);

% HVAC
%--------------------------------------------------------------------------
for u = 1:m
    B(:,:,u)  = zeros(m,m);
    B(u,:,u)  = 1;
end
level(1).states(3).name  = 'HVAC';
level(1).states(3).state = {'hot','warm','off','chill','cold'};
level(1).states(3).B     = B;
level(1).states(3).iD    = {};
level(1).states(3).iE    = {};
level(1).states(3).H     = [];
level(1).states(3).U     = U(2);

% Daily rate
%--------------------------------------------------------------------------
level(1).states(4).name  = 'Daily_Rate';
level(1).states(4).state = {'hot','warm','off','chill','cold'};
level(1).states(4).B(:,:,1) = eye(4,4);
level(1).states(4).iD   = {'Daily_Rate'};
level(1).states(4).iE   = {};
level(1).states(4).H    = [];
level(1).states(4).U    = 0;

% Cloud cover
%--------------------------------------------------------------------------
level(1).states(5).name  = 'Cover';
level(1).states(5).state = {'yes','clear'};
level(1).states(5).B(:,:,1) = [4, 8; 1, 1];
level(1).states(5).B(:,:,2) = [1, 1; 1, 2];
level(1).states(5).iD   = {};
level(1).states(5).iE   = {'Cover'};
level(1).states(5).H    = [];
level(1).states(5).U    = 0;

% External temperature
%--------------------------------------------------------------------------
level(1).states(6).name  = 'Temperature';
level(1).states(6).state = num2cell(1:n);
level(1).states(6).B(:,:,1) = spm_conv(eye(n,n),2);
level(1).states(6).iD   = {'Temperature'};
level(1).states(6).iE   = {};
level(1).states(6).H    = [];
level(1).states(6).U    = 0;


% Room occupancy
%--------------------------------------------------------------------------
level(1).states(7).name  = 'Occupied';
level(1).states(7).state = {'no','yes'};
level(1).states(7).B(:,:,1) = diag([8,1]) + 1;
level(1).states(7).B(:,:,2) = diag([1,8]) + 1;
level(1).states(7).iD   = {};
level(1).states(7).iE   = {'Occupied'};
level(1).states(7).H    = [];
level(1).states(7).U    = 0;


% second (day) level
%==========================================================================

% Daily Rate
%--------------------------------------------------------------------------
level(2).outcomes(1).name  = 'Daily_Rate';
level(2).outcomes(1).state = num2cell(1:4);
level(2).outcomes(1).A     = spm_softmax(rand(4,7)*8);
level(2).outcomes(1).iA    = {'Day'};
level(2).outcomes(1).C     = [];
level(2).outcomes(1).iC    = {};

% Hours
%--------------------------------------------------------------------------
level(2).outcomes(2).name  = 'Hour';
level(2).outcomes(2).state = num2cell(2:2:24);
level(2).outcomes(2).A     = sparse(1,1:7,1,12,7);
level(2).outcomes(2).iA    = {'Day'};
level(2).outcomes(2).C     = [];
level(2).outcomes(2).iC    = {};

% External (x3) temperature
%--------------------------------------------------------------------------
level(2).outcomes(3).name  = 'Temperature';
level(2).outcomes(3).state = num2cell(1:8);
level(2).outcomes(3).A     = eye(n,n);
level(2).outcomes(3).iA    = {'Weather'};
level(2).outcomes(3).C     = [];
level(2).outcomes(3).iC    = {};

% Room occupancy
%--------------------------------------------------------------------------
level(2).outcomes(4).name   = 'Occupied';
level(2).outcomes(4).state  = {'no','yes'};
level(2).outcomes(4).A(1,:) = 1:n;
level(2).outcomes(4).A(2,:) = flip(1:n);
level(2).outcomes(4).iA     = {'Weather'};
level(2).outcomes(4).C      = [];
level(2).outcomes(4).iC     = {};

% Cloud cover
%--------------------------------------------------------------------------
level(2).outcomes(5).name   = 'Cover';
level(2).outcomes(5).state  = {'yes','no'};
level(2).outcomes(5).A(1,:) = flip(1:n);
level(2).outcomes(5).A(2,:) = 1:n;
level(2).outcomes(5).iA     = {'Weather'};
level(2).outcomes(5).C      = [];
level(2).outcomes(5).iC     = {};


% Calender time
%--------------------------------------------------------------------------
level(2).states(1).name  = 'Day';
level(2).states(1).state = {'Mo','Tu','We','Th','Fr','Sa','Su'};
level(2).states(1).B     = spm_speye(7,7,-1,1);
level(2).states(1).iD    = {};
level(2).states(1).iE    = {};
level(2).states(1).H     = [];
level(2).states(1).U     = 0;


% Weather
%--------------------------------------------------------------------------
level(2).states(2).name  = 'Weather';
level(2).states(2).state = num2cell(1:n);
level(2).states(2).B     = spm_conv(eye(n,n),n/2);
level(2).states(2).iD    = {};
level(2).states(2).iE    = {};
level(2).states(2).H     = [];
level(2).states(2).U     = 0;

% Create MDP structure and initialise day and hour
%--------------------------------------------------------------------------
MDP      = spm_make_MDP(level,[12,D],[N,0]);
MDP.D{1} = sparse(1,1,1,7,1);
MDP.D{2} = sparse(1,1,1,8,1);


% Illustrate with an example over D days
%==========================================================================
MDP = spm_MDP_VB_XXX(MDP);
 
% illustrate active inference (second level)
%--------------------------------------------------------------------------
spm_figure('GetWin','Weeky'); clf
spm_MDP_VB_trial(MDP);
 
% illustrate active inference (first level: last epoch)
%--------------------------------------------------------------------------
spm_figure('GetWin','Last week'); clf
spm_MDP_VB_trial(MDP.MDP,1:7,1:6);

% illustrate both by retrieving recorded posteriors and outcomes in MDP.Q
%--------------------------------------------------------------------------
spm_show_Q(MDP)

return


function [MDP,RDP] = spm_make_MDP(level,T,N)
% create deep (Recursive) MDP from level structures
% FORMAT [MDP,RDP] = spm_make_MDP(level,T,N)
%__________________________________________________________________________


% defaults
%--------------------------------------------------------------------------
Nm  = numel(level);
if nargin < 2, T = zeros(1,Nm) + 2; end
if nargin < 3, N = zeros(1,Nm) + 0; end

% for each hierarchical level
%--------------------------------------------------------------------------
for n = 1:Nm

    % states at this level
    %======================================================================
    for f = 1:numel(level(n).states)

        % for each factor
        %------------------------------------------------------------------
        MDP.B{f} = full(level(n).states(f).B);
        MDP.U(f) = full(level(n).states(f).U);
        MDP.H{f} = full(level(n).states(f).H);

        % labels
        %------------------------------------------------------------------
        MDP.label.factor{f} = level(n).states(f).name;
        if numel(level(n).states(f).state)
           MDP.label.name{f} = level(n).states(f).state;
        end

        % parents
        %------------------------------------------------------------------
        if n < Nm
            
            % inital states
            %--------------------------------------------------------------
            state = level(n).states(f).iD;
            [~,j] = ismember(state,{level(n + 1).outcomes.name});
            MDP.id.D{f} = j;

            % inital paths
            %--------------------------------------------------------------
            paths = level(n).states(f).iE;
            [~,j] = ismember(paths,{level(n + 1).outcomes.name});
            MDP.id.E{f} = j;

        end
    end

    % outcomes at this level
    %======================================================================
    for g = 1:numel(level(n).outcomes)

        % for each factor
        %------------------------------------------------------------------
        MDP.A{g} = level(n).outcomes(g).A;
        MDP.C{g} = level(n).outcomes(g).C;

        % labels
        %------------------------------------------------------------------
        MDP.label.modality{g} = level(n).outcomes(g).name;
        if numel(level(n).outcomes(g).state)
           MDP.label.outcome{g} = level(n).outcomes(g).state;
        end

        % parents of likelihoods
        %------------------------------------------------------------------
        state = level(n).outcomes(g).iA;
        [~,j] = ismember(state,{level(n).states.name});
        MDP.id.A{g} = j;

        % parents of contraints
        %------------------------------------------------------------------
        state = level(n).outcomes(g).iC;
        [~,j] = ismember(state,{level(n).states.name});
        MDP.id.C{g} = j;

    end

    % MDP at this level
    %======================================================================
    MDP.N   = N(n);
    MDP.T   = T(n);
    MDP.L   = n;
    if n < Nm
        PDP.MDP = MDP;
        MDP     = PDP;
    end

end

% hierachical form and checking
%--------------------------------------------------------------------------
RDP = spm_check_edges(MDP);

return


function spm_show_Q(MDP)
% plot latent states and outcomes of subordinate states in MDP.Q
% FORMAT [MDP,RDP] = spm_make_MDP(level,T,N)
% MDP - recursive MDP structure after inversion
%__________________________________________________________________________


% Get cell array form and work through levels
%--------------------------------------------------------------------------
PDP   = spm_check_edges(MDP);
NL    = MDP.L - 1;
for L = 1:NL
    
    % figure and sizes
    %----------------------------------------------------------------------
    spm_figure('GetWin',sprintf('Level %i',L)); clf
    Nf  = size(MDP.Q.X{L},1);
    Ng  = size(MDP.Q.Y{L},1);

    % latent states and posteriors
    %----------------------------------------------------------------------
    for f = 1:Nf
        subplot(Nf + Ng + 2,1,f)
        imagesc(1 - spm_cat(MDP.Q.X{L}(f,:)))
        hold on, plot(MDP.Q.s{L}(f,:),'.r','MarkerSize',16)
        title(PDP{L}.label.factor{f})
    end

    subplot(Nf + Ng + 2,1,Nf + 2)
    plot(MDP.Q.E{L},'-r','LineWidth',4),spm_axis tight
    title('ELBO','FontSize',16)

    % latent states and posteriors
    %----------------------------------------------------------------------
    for g = 1:Ng
        subplot(Nf + Ng + 2,1,Nf + 2 + g)
        imagesc(1 - spm_cat(MDP.Q.Y{L}(g,:)))
        hold on, plot(MDP.Q.o{L}(g,:),'.r','MarkerSize',16)
        title(PDP{L}.label.modality{g})
    end

    % cost analyis
    %======================================================================
    spm_figure('GetWin',sprintf('Cost %i',L)); clf

    % Get modalities with constraints
    %----------------------------------------------------------------------
    C     = [];
    for g = 1:numel(PDP{L}.C)
        if numel(PDP{L}.C{g})
            if any(diff(PDP{L}.C{g}),'all')
                C(g) = 1;
            end
        end
    end
    C     = find(C);
    Ng    = numel(C);

    % Predictive posteriors and outcomes
    %----------------------------------------------------------------------
    o     = [];
    X     = {};
    P     = {};
    for i = 1:Ng
        g = C(i);

        % Predictive posteriors and outcomes
        %------------------------------------------------------------------
        o(i,:) = MDP.Q.o{L}(g,:);
        subplot(Ng + Ng,1,i)
        imagesc(1 - spm_cat(MDP.Q.Y{L}(g,:)))
        hold on, plot(o(i,:),'.r','MarkerSize',16)
        title(PDP{L}.label.modality{g})

        % Preference distributions
        %------------------------------------------------------------------
        if numel(PDP{L}.id.C{g})
            f      = PDP{L}.id.C{g};
            X{i}   = spm_cat(MDP.Q.X{L}(f,:));
        else
            X{i}   = {};
        end
        P{i} = spm_dir_norm(PDP{L}.C{g});
    end

    % Evaluate costs
    %----------------------------------------------------------------------
    for t = 1:size(o,2)
        for i = 1:size(o,1)
            p = PDP{L}.C{C(i)};
            if numel(X{i})
                U(i,t) = 0;
                for j = 1:size(X{i},1)
                    U(i,t) = U(i,t) - spm_log(p(o(i,t),j))*X{i}(j,t);
                end
            else
                U(i,t) = -spm_log(p(o(i,t)));
            end

        end
    end

    % Plot costs and free energy
    %----------------------------------------------------------------------
    subplot(6,1,4), plot(-MDP.Q.E{L},'-k','LineWidth',2), spm_axis tight
    title('Variational free energy','FontSize',16), spm_axis tight
    hold on, plot([0,t],[3,3],':r'), plot([0,t],[5,5],'-.r'), hold off

    subplot(6,1,5), plot(U','LineWidth',2), legend(PDP{L}.label.modality(C))
    title('Cost'), ylabel('Nats'), xlabel('time'), spm_axis tight
    hold on, plot([0,t],[3,3],':r'), plot([0,t],[5,5],'-.r'), hold off

    subplot(6,1,6), plot(sum(U),'r','LineWidth', 2)
    title('Cost (total)'), ylabel('Nats'), xlabel('time'), spm_axis tight
    hold on, plot([0,t],[3,3],':r'), plot([0,t],[5,5],'-.r'), hold off

end

return

