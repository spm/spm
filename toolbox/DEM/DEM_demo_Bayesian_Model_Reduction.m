function RCM = DEM_demo_Bayesian_Model_Reduction
% This demonstration code illustrates the application of post hoc model 
% optimisation or Bayesian model reduction (BMR) in identifying gene and 
% gene-gene interaction effects in behavioural or physiological variables.  
% The basic idea is to replace conventional heuristics based on the 
% assumption that the contribution of any gene is sampled from a sparse 
% distribution (such that a small number contribute and a large number do 
% not) with an explicit search over a model space that includes all sparse
% and non-sparse models.  This exhaustive search rests upon recent 
% advances in model optimisation based upon variational Bayesian model 
% inversion.  In short, it is possible to estimate the posterior 
% distribution of model parameters under a reduced model, given the 
% posterior and prior distributions under a full model.  
% In this context, a reduced model corresponds to a model in which some 
% parameters are removed (by shrinking their prior variance to zero).  
% This means that it is only necessary to invert the full model and then 
% perform automatic BMR (over all possible combinations of parameters) 
% using a greedy search based upon the free energy approximation to log 
% model evidence.  With sufficient signal to noise, this scheme can 
% recover the small number of effects, even in under determined or 
% ill-posed problems (where the number of potential effects can vastly 
% exceed the number of samples). 
%
% The illustration below uses 128 subjects who have been measured three 
% times (say in three brain regions) and we want to model these 
% measurements in terms of first and second order genetic contributions 
% given 8 (binary) genetic variables and all (unique) pair wise 
% interactions.  This means that there are 36 unknown parameters 
% (excluding a constant and, say, age confounds over subjects).  In the 
% scheme below, each measurement is inverted separately under a simple 
% (polynomial) model with uninformative priors on the parameters and 
% (precision) hyper-parameters describing beliefs about signal to noise.  
% A fixed effects Bayesian model averaging (BMA) scheme is used in 
% combination with BMR to identify the best model out of all possible 
% combinations of first and second order effects.  With the signal to 
% noise and number of samples used in this simulation, the recovery is 
% generally perfect.  This scheme also illustrates inference over a 
% partition of model space (or families of models).

%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_Bayesian_Model_Reduction.m 6306 2015-01-18 20:50:38Z karl $

rng('default')

% Genomic data (G)
%--------------------------------------------------------------------------
Ns = 256;                                   % unmber of subjects
Ng = 8;                                    % unmber of genes
Nr = 3;                                    % unmber of (fMRI) regions

% Genotype
%--------------------------------------------------------------------------
G     = rand(Ns,Ng) > 1/2;

% design matix (U) with second-order (gene-gene interaction) effects
%--------------------------------------------------------------------------
GG    = [];
for i = 1:Ns
    u     = [];
    for j = 1:Ng
        for k = (j + 1):Ng
          u(end + 1) = G(i,j)*G(i,k);
        end
    end
    GG = [GG; u];
end
U  = [G GG];
Nb = size(U,2);
R  = ones(1,Nr);
A  = zeros(Nb,1);

% first-order effects (g1 abd g4) and one interaction (g1 x g2)
%--------------------------------------------------------------------------
A([1 4 (Ng + 3)],:) = [1/2 1 1];             % effect sizes

% Simulate data with a SNR of about 4:1
%--------------------------------------------------------------------------
y  = U*A*R + randn(Ns,Nr);                 % SD(noise) = 1


% Bayesian model inversion
%==========================================================================

% Model specification
%--------------------------------------------------------------------------
M.IS   = @(P,M,U) U*P.A;
M.pE.A = zeros(Nb,1);
M.pC.A = zeros(Nb,1) + 128;
M.hE   = 0;
M.hC   = 1/128;

% confounds (at the between subject level)
%--------------------------------------------------------------------------
age   = rand(Ns,1);
Y.X0  = [age ones(Ns,1)];

for i = 1:Nr
    
    % model inversion for this region
    %----------------------------------------------------------------------
    Y.y       = y(:,i);
    [Ep,Cp]   = spm_nlsi_GN(M,U,Y);
    
    % save model for subsequent BMS
    %----------------------------------------------------------------------
    DCM{i}.M  = M;
    DCM{i}.Ep = Ep;
    DCM{i}.Cp = Cp;
    
end

% Bayesian model reduction: see below for family fun
%==========================================================================
RCM  = spm_dcm_post_hoc(DCM,@fun);
    

% show results
%--------------------------------------------------------------------------
spm_figure('Getwin','Model posterior (over families)'); clf

Qf         = zeros(1,8);
Qf(fun(A)) = 1;                           % true family
Pf         = zeros(1,8);
Pf(1:length(RCM.Pf)) = RCM.Pf;            % psoteror family probailities

subplot(2,1,1)
bar([Pf; Qf]')
xlabel('familiy')
title('Model posterior (over families)','FontSize',16)
axis square
legend({'esimated','true'})

subplot(2,1,2)
bar([RCM.Ep.A A])
title('Recovered and true effects','FontSize',16)
xlabel('parameter')
axis square


% Partition of mdoel space into famllies
%==========================================================================
function k = fun(A)

f1  = any(A(1:4));
f2  = any(A(5:8));
f3  = any(A(9:12));

if ~f1 && ~f2 && ~f3, k = 1; end
if ~f1 && ~f2 &&  f3, k = 2; end
if ~f1 &&  f2 && ~f3, k = 3; end
if ~f1 &&  f2 &&  f3, k = 4; end
if  f1 && ~f2 && ~f3, k = 5; end
if  f1 && ~f2 &&  f3, k = 6; end
if  f1 &&  f2 && ~f3, k = 7; end
if  f1 &&  f2 &&  f3, k = 8; end

