function DCM = DEM_Diabetes
% dynamic causal model of Indian-type diabetes
%__________________________________________________________________________
% 
% This demonstration routine illustrates a dynamic causal model of
% Indian-type diabetes, with a focus on genetic and epigenetic mechanisms.
% Technically, it illustrates the use of population density dynamics
% modelled with a Master Equation; namely, a probability transition matrix
% that determines the probability of moving from one state to another.
% Crucially, these states are discrete, assuming that a person sampled at
% random is in one of a finite number of states. These states can be
% factored into a number of relevant dimensions (e.g., age, genotype, et
% cetera), creating a state space of cells in a high dimensional array. The
% population dynamics are determined by the (sparse) probabilities of
% moving from one cell to another cell, generally expressed as rate
% constants (i.e., the inverse of the time a particular cell is occupied).
% Equipped with the transition matrix one can then solve the master
% equation over time where, some parameters can change in response to
% exogenous causes. This provides the basis for a generative model that is
% completed by a likelihood mapping from the probability distribution over
% states to observations over time. Generally, these observations are
% nonnegative prevalences that can be modelled with a Gaussian likelihood,
% based upon a binomial approximation. In this example, the ensuing
% generative model is inverted using variational Laplace: a generic
% variational scheme that uses a gradient descent on variational free
% energy under Gaussian assumptions about free parameters and
% hyperparameters; namely, transition rates and the precision of random
% sampling effects. Priors over the model parameters are specified in terms
% of their means and variances.
%
% This demo sets up the priors over parameters that are used to generate
% timeseries to which appropriate noise is added. These synthetic data are
% inverted to illustrate the recovery of the parameters from sparse data
% that are sampled asynchronously; here, every few years over a century.
%
% The generative model includes all the states that are necessary to
% generate population demographics, clinical data and metrics of economic
% or agricultural activity. In turn, fluctuations in the state of the
% population are caused by discrete events (i.e., famines) that change the
% transition from a food secure to an insecure state. The effects of food
% insecurity are modelled through an epigenetic mechanism, in which
% interactions between the environment (food security) and genetics (high
% versus low genotypic risk) are mediated at the point of birth. This
% mediation is modelled by a simple form of replicator dynamics, in which
% the rate of moving from a gene pool to a particular phenotype (high,
% typical or low BMI) depends on both the environment (food security) and
% genotype. This replication rate can be regarded as reproductive fitness,
% where having a low risk (i.e., normal) genotype reduces reproductive
% fitness when, and only when food insecure. The replicator dynamics are
% situated in transitions among age groups, from child to adults and then
% old. Mortality entails moving to the gene pool that slowly equilibria
% rates over time. This means that replicator dynamics are mediated by
% replication rates that themselves depend upon genetic-environment
% interactions. The clinical phenotype inherits from the BMI phenotype, in
% the sense that the probability of moving from healthy to diabetic is
% increased for a low BMI phenotype (when at high genetic risk). This model
% also features direct environment-phenotype interactions; for example,
% food security determines the transitions between BMI phenotypes, with the
% potential risk of obesity in food secure states. In addition, under food
% insecurity, there is a probability of dying and moving to the pool state.
% The transition rates among the age groups determine the number of people
% in each group and the composition of phenotypes and genotypes that enter
% the genetic pool.
%
% Although a gross simplification, this model features the minimum number
% of states necessary to model interactions between genotype, phenotype and
% environment in a way that incorporates age-dependent effects at the
% population level. With the probability of being in any particular
% combination of such states, one can then generate population data that
% are, in principle, measurable; thereby affording the opportunity to
% estimate the underlying genetic and epigenetic contributions to disease
% burden, where the genetic and epigenetic legacies are expressed within
% and between generations, respectively.
%__________________________________________________________________________
 

% latent states and outcomes
%==========================================================================
M.factor    = {'Age', 'Phenotype', 'Environment', 'Genotype', 'Clinical'};
%--------------------------------------------------------------------------
M.states{1} = {'Child',     'Adult',    'Old', 'Pool'};
M.states{2} = {'Low BMI',   'Typical', 'High BMI'};
M.states{3} = {'Secure',    'Insecure'};
M.states{4} = {'High Risk', 'Low Risk'};
M.states{5} = {'Healthy',   'Diabetic'};
%--------------------------------------------------------------------------
M.data = {'Child','Adult','Old','Economic','Production','BMI','Diabetes'};

% famine onsets (years)
%--------------------------------------------------------------------------
U.T  = [16, 32];

% metadata: sample times in years for each data modality
%--------------------------------------------------------------------------
M.T{1}  = 1:10:121;            % number of children (1000)
M.T{2}  = 1:10:121;            % number of adults
M.T{3}  = 1:10:121;            % number of elderly
M.T{4}  = 1:20:121;            % economic activity (percent)
M.T{5}  = 1:20:121;            % food production (percent)
M.T{6}  = 100:2:121;           % average BMI
M.T{7}  = 100:2:121;           % Diabetes (percent)

for g = 1:numel(M.T)
    N.T{g} = 1:128;
end

% likelihood and prior parameters
%==========================================================================
[pE,pC] = spm_diabetes_priors(U);

% illustrate steady state
%--------------------------------------------------------------------------
spm_diabetes_T(pE);

% generate data
%==========================================================================
[Y,X] = spm_diabetes_gen(pE,M,U);
Ng    = numel(Y);
Nf    = numel(X);

% add random effects
%--------------------------------------------------------------------------
for g = 1:Ng
    V(g) = mean(sqrt(Y{g}));
    Y{g} = Y{g} + sqrt(V(g))*randn(size(Y{g}));
end

% plot
%--------------------------------------------------------------------------
spm_diabetes_plot(Y,X,M)


% data structure with precision components
%==========================================================================
xY.y  = spm_vec(Y);
C     = cell(1,Ng);
hE    = zeros(1,Ng);
for g = 1:Ng
    C{g} = zeros(1,numel(Y{g}));
end
for g = 1:Ng

    % precision components
    %----------------------------------------------------------------------
    Q       = C;
    Q{g}    = Q{g} + 1;
    xY.Q{g} = diag(spm_cat(Q));

    % log precisions based upon total counts
    %----------------------------------------------------------------------
    h     = mean(sqrt(Y{g} + 1));
    hE(g) = 4 - log(h);
end

% feature selection
%--------------------------------------------------------------------------
FS     = @(Y)real(sqrt(spm_vec(Y)));

% model specification
%==========================================================================
M.Nmax = 4;                    % maximum number of iterations
M.G    = @spm_diabetes_gen;    % generative function
M.FS   = FS;                   % feature selection  (link function)
M.pE   = pE;                   % prior expectations (parameters)
M.pC   = pC;                   % prior covariances  (parameters)
M.hE   = hE;                   % prior expectation  (log-precision)
M.hC   = 1/8;                  % prior covariances  (log-precision)

% model inversion with Variational Laplace (Gauss Newton)
%==========================================================================
[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,xY);

% save in DCM structure
%--------------------------------------------------------------------------
DCM.M  = M;
DCM.Y  = Y;
DCM.U  = U;
DCM.Ep = Ep;
DCM.Eh = Eh;
DCM.Cp = Cp;
DCM.F  = F;

% Illustrate estimated states and measurements
%--------------------------------------------------------------------------
M.T   = N.T;
[Y,X] = spm_diabetes_gen(Ep,M,U);
spm_diabetes_plot(Y,X,M)

% show parameter recovery
%--------------------------------------------------------------------------
spm_figure('GetWin','Latent states'), clf

Vp      = spm_unvec(diag(Cp),pE);
qP.P{1} = Ep;
qP.V{1} = Vp;
pP.P{1} = pE;

spm_DEM_qP(qP,pP)


return


function [P,C,str] = spm_diabetes_priors(U)
% Generate log transformed parameters
% FORMAT [P,C,str] = spm_diabetes_priors
% P    - prior expectations
% C    - prior variances
%__________________________________________________________________________

% expectations
%--------------------------------------------------------------------------
P.size      = 1;               % population size
P.diabetes  = 1/128;           % risk of diabetes | not low BMI
P.excess    = 1/4;             % risk of diabetes | low BMI
P.famine    = 1/32;            % P(food insecure | food secure)
P.recovery  = 1/8;             % P(food secure | food insecure)
P.nutrition = 1/4;             % P(normalisation of BMI | food secure)
P.obesity   = 1/8;             % P(obesity | food secure)
P.mortality = 1/16;            % P(morbidity | food insecure)
P.fecundity = 1/2;             % P(replication | insecure and Low risk) 
P.drift     = 1/128;           % genetic exchange between gene pools
P.prod      = 1/8;             % sensitivity of production of food insecurity
P.severity  = ones(size(U.T)); % severity of successive famines

% Variances (on log transformed parameters
%--------------------------------------------------------------------------
V     = exp(-2);               % uninformative priors
W     = exp(-4);               % informative priors (weak)

C.size      = V;
C.diabetes  = V;
C.excess    = W;
C.famine    = W;
C.recovery  = W;
C.nutrition = W;
C.obesity   = W;
C.mortality = W;
C.fecundity = W;
C.drift     = W;
C.prod      = W;
C.severity  = ones(size(U.T))*V;


% log transform
%==========================================================================
P   = spm_vecfun(P,@log);

% field names of random effects
%--------------------------------------------------------------------------
str = fieldnames(P);

return

function [B,Ns] = spm_diabetes_T(P)
% Generate parameterised probability transitions for Master equation
% FORMAT [B,Ns] = spm_diabetes_T(P)
% P    - model parameters
%__________________________________________________________________________

% factors and their states
%--------------------------------------------------------------------------
factor    = {'Age', 'Phenotype', 'Environment', 'Genotype', 'Clinical'};
states{1} = {'Child',     'Adult',   'Old', 'Pool'};
states{2} = {'Low BMI',   'Typical', 'High BMI'};
states{3} = {'Secure',    'Insecure'};
states{4} = {'High Risk', 'Low Risk'};
states{5} = {'Healthy',   'Diabetic'};
%__________________________________________________________________________

% exponentiate
%--------------------------------------------------------------------------
P     = spm_vecfun(P,@exp);               % ensure all P > 0
P     = spm_vecfun(P,@erf);               % ensure all P < 1

% transition priors
%--------------------------------------------------------------------------
Ns    = [4 3 2 2 2];                      % number of states per factor
NS    = prod(Ns);                         % number of joint states
B     = zeros(NS,NS);                     % transition matrix
D     = 88;                               % Gene pool occupancy

% specify generative model under combinations of hidden states
%--------------------------------------------------------------------------
S     = spm_combinations(Ns);
for s = 1:NS

    % Priors (B)
    %======================================================================
    
    % risk of developing Diabetes
    %----------------------------------------------------------------------
    if S(s,5) == 1 && ismember(S(s,1),[2 3])
        
        % if healthy & adult or old, move to diabetic
        %------------------------------------------------------------------
        T      = S(s,:);
        T(5)   = 2;
        t      = ismember(S,T,'rows');

        % low BMI and high genetic risk has excess risk
        %------------------------------------------------------------------
        if S(s,2) == 1 && S(s,4) == 1
            B(t,s) = P.excess;
        else
            B(t,s) = P.diabetes;
        end
    end

    % risk of food insecurity
    %----------------------------------------------------------------------
    if S(s,3) == 1
        
        % if food secure, move to insecure
        %------------------------------------------------------------------
        T      = S(s,:);
        T(3)   = 2;
        t      = ismember(S,T,'rows');
        B(t,s) = P.famine;

    else

        % if food insecure, move to secure
        %------------------------------------------------------------------
        T      = S(s,:);
        T(3)   = 1;
        t      = ismember(S,T,'rows');
        B(t,s) = P.recovery;

    end

    % risk of malnutrition
    %----------------------------------------------------------------------
    if S(s,2) == 1 && S(s,3) == 1
        
        % if low BMI and food secure, move to normal
        %------------------------------------------------------------------
        T      = S(s,:);
        T(2)   = 2;
        t      = ismember(S,T,'rows');
        B(t,s) = P.nutrition;

     elseif S(s,2) == 2 && S(s,3) == 1

        % if Typical BMI and food secure, move to High BMI
        %------------------------------------------------------------------
        T      = S(s,:);
        T(2)   = 3;
        t      = ismember(S,T,'rows');
        B(t,s) = P.obesity;


    elseif S(s,2) == 3 && S(s,3) == 2

        % if high BMI and food insecure, move to normal
        %------------------------------------------------------------------
        T      = S(s,:);
        T(2)   = 2;
        t      = ismember(S,T,'rows');
        B(t,s) = P.nutrition;

    end

    % and pool if food insecure
    %------------------------------------------------------------------
    if ismember(S(s,3),[1 2 3])
        T      = S(s,:);
        T(1)   = 4;
        t      = ismember(S,T,'rows');
        B(t,s) = P.mortality;
    end

    % population dynamics
    %======================================================================
    if S(s,1) == 1
        
        % if born move to adult (8+ years)
        %------------------------------------------------------------------
        T      = S(s,:);
        T(1)   = 2;
        t      = ismember(S,T,'rows');
        B(t,s) = 1/8;


    elseif S(s,1) == 2

        % if adult move to old (38+ years)
        %------------------------------------------------------------------
        T      = S(s,:);
        T(1)   = 3;
        t      = ismember(S,T,'rows');
        B(t,s) = 1/30;
        

     elseif S(s,1) == 3

        % if old move to gene pool (D+ years)
        %------------------------------------------------------------------
        T      = S(s,:);
        T(1)   = 4;
        t      = ismember(S,T,'rows');
        B(t,s) = 1/50;

    elseif S(s,1) == 4

        % replicator dynamics with epigenetic interactions
        %==================================================================
        if S(s,3) == 1                        % food secure

            T      = S(s,:);
            T(1)   = 1;                       % Child
            T(2)   = 2;                       % Typical
            T(5)   = 1;                       % healthy
            t      = ismember(S,T,'rows');
            B(t,s) = 1/D;

        elseif S(s,3) == 2                        % food insecure

            if S(s,4) == 1

                % High genetic risk, move to Low BMI program
                %----------------------------------------------------------
                T      = S(s,:);
                T(1)   = 1;                       % Child
                T(2)   = 1;                       % Low BMI
                T(5)   = 1;                       % healthy
                t      = ismember(S,T,'rows');
                B(t,s) = 1/D;

            else

                % move to typical at reduced replication rate
                %----------------------------------------------------------
                T      = S(s,:);
                T(1)   = 1;                       % Child
                T(2)   = 2;                       % Typical
                T(5)   = 1;                       % healthy
                t      = ismember(S,T,'rows');
                B(t,s) = P.fecundity/D;

            end

        end

        % genetic equilibration
        %----------------------------------------------------------------------
        if S(s,4) == 1

            % if high risk move to low risk
            %------------------------------------------------------------------
            T      = S(s,:);
            T(4)   = 2;
            t      = ismember(S,T,'rows');
            B(t,s) = P.drift;

        elseif S(s,4) == 2

            % if high risk move to low risk
            %------------------------------------------------------------------
            T      = S(s,:);
            T(4)   = 1;
            t      = ismember(S,T,'rows');
            B(t,s) = P.drift;

        end

    end

end


% create probability transition matrix (c.f., Laplacian)
%==========================================================================
B = B + diag(1 - sum(B));

if nargout, return, end

% illustrate
%--------------------------------------------------------------------------
spm_figure('GetWin','Diabetes DAG'), clf

[e,s] = eig(B,'nobalance');
[~,j] = max(diag(s));
p     = spm_dir_norm(e(:,j));
p     = reshape(p,Ns);
m     = spm_marginal(p);
for f = 1:numel(m)
    subplot(4,3,f), bar(m{f}), 
    title(factor{f}),set(gca,'XtickLabel',(states{f}))
    axis square
end

subplot(4,3,6)
imagesc(B), axis square, title('Dynamics')

% directed graph with transitive closure
%--------------------------------------------------------------------------
subplot(4,2,7)
plot(transclosure(digraph(B > 1/16)))

subplot(4,2,8)
T  = -1./log(abs(diag(s)));
T  = sort(T,'descend');
bar(T), set(gca,'YLim',[0,128])
title('Time constants'), xlabel('eigenmode'), ylabel('years'),

% display marginal dynamics
%--------------------------------------------------------------------------
b     = reshape(B,[Ns,Ns]);
Nf    = numel(Ns);
for f = 1:Nf
    subplot(4,Nf,2*Nf + f)
    i    = true(1,Nf);
    i(f) = false;
    i    = [i,i];
    BF   = squeeze(spm_sum(b,find(i)));
    BF   = log(BF);
    imagesc(spm_dir_norm(BF)), axis square
    title(factor{f})
end

return


function [Y,X] = spm_diabetes_gen(P,M,U)
% Generate predictions and hidden states of a COVID model
% FORMAT [Y,X] = spm_SARS_gen(P,M,U,NPI,age)
% P    - model parameters
% M    - model structure
% U    - output times
%__________________________________________________________________________


% duration of time series
%--------------------------------------------------------------------------
if isfield(M,'T')
    tmax = max(spm_vec(M.T));
else
    tmax = 128;
end

% exponentiate
%--------------------------------------------------------------------------
Q      = spm_vecfun(P,@exp);                 % ensure all Q > 0
Q      = spm_vecfun(Q,@erf);                 % ensure all Q < 1

% inital steady state distribution
%--------------------------------------------------------------------------
[T,Ns] = spm_diabetes_T(P);
[e,s]  = eig(T,'nobalance');
[~,j]  = max(diag(s));
p      = spm_dir_norm(e(:,j));

P0    = P.famine;
X     = {};
for t = 1:tmax

    % time dependent parameters
    %----------------------------------------------------------------------
    u  = find(ismember(U.T,t));
    if numel(u)
        P.famine = P.severity(u);
    else
        P.famine = P0;
    end

    % update ensemble density (x)
    %======================================================================
    T   = spm_diabetes_T(P);
    p   = T*p;

    % marginals for population
    %----------------------------------------------------------------------
    n     = reshape(p,Ns);
    n     = spm_marginal(n);

    % marginal densities (p): removing inactive pool
    %----------------------------------------------------------------------
    m     = reshape(p,Ns);
    m(end,:,:,:,:,:,:) = 0;
    m     = spm_marginal(m);
    for j = 1:numel(m)
        m{j}      = spm_dir_norm(m{j});
        X{j}(t,:) = m{j};
    end
    
    % outcomes
    %======================================================================

    % age and other demographics
    %----------------------------------------------------------------------
    Y{1}(t,1)  = 1000*Q.size*n{1}(1);          % number of children (1000)
    Y{2}(t,1)  = 1000*Q.size*n{1}(2);          % number of adults
    Y{3}(t,1)  = 1000*Q.size*n{1}(3);          % number of elderly
    Y{4}(t,1)  = 100*m{3}(1)*m{1}(2);          % economic activity (percent)
    Y{5}(t,1)  = 100*[1,Q.prod]*m{3};          % food production (percent)
    Y{6}(t,1)  = [18,24,32]*m{2};              % average BMI
    Y{7}(t,1)  = 100*m{5}(2);                  % Diabetes (percent)

end

% subsample times in years for each data modality
%--------------------------------------------------------------------------
if isfield(M,'T')
    for i = 1:numel(M.T)
        Y{i} = Y{i}(M.T{i}(:));
    end
end

return

% plot
%==========================================================================
function spm_diabetes_plot(Y,X,M)
% plotting subroutine
%__________________________________________________________________________

% plot data
%--------------------------------------------------------------------------
spm_figure('GetWin','Diabetes data'); clf
Ng    = numel(Y);
Nf    = numel(X);
for g = 1:Ng
    subplot(Ng,1,g)
    plot(M.T{g},Y{g},'.','MarkerSize',16),
    title(M.data{g}), ylabel('data'), set(gca,'XLim',[1,128])
end
xlabel('years')

% plot latent states
%--------------------------------------------------------------------------
spm_figure('GetWin','Latent states'); clf

for f = 1:Nf
    subplot(Nf,1,f)
    plot(X{f}), title(M.factor{f})
    ylabel('proportion'), legend(M.states{f})
end
xlabel('years')
