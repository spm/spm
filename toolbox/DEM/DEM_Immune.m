function [DCM, PEB, BMA, BMR] = DEM_Immune
% This demo builds upon the results of the COVID modelling demos, which
% found that the epidemic data could best be explained by appealing to the
% idea that much of the population may not be susceptible and that of those
% who are, some may be resistant and only experience a mild illness.
% This means measures of immunity based upon antibody tests may
% underestimate the effective herd immunity. This demo formalises several
% alternative hypotheses as to the mechanisms of resistance. It
% demonstrates the way in which the underlying model may be inverted to 
% test these hypotheses. 
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
 
% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% For reproducibility
%--------------------------------------------------------------------------
rng default

% Specify priors and simulate 'normal' immune response
%==========================================================================
[pE,pC] = spm_immune_priors;
P{1} = pE;

spm_figure('GetWin','Immune response'); clf
y{1} = spm_immune_plot(P{1},spm_unvec(spm_vec(pC)/16,pC));

% Reduce viral entry to cells
%==========================================================================
P{2} = pE;
P{2}.int = -2.2;

spm_figure('GetWin','Reduced cell entry'); clf
y{2} = spm_immune_plot(P{2},spm_unvec(spm_vec(pC)/16,pC));

% Pre-existing immunity
%==========================================================================
P{3}   = pE;
P{3}.n = 2;

spm_figure('GetWin','Cross-reactive immunity'); clf
y{3} = spm_immune_plot(P{3},spm_unvec(spm_vec(pC)/16,pC));

% T-cell mediated immunity
%==========================================================================
P{4} = pE;
P{4}.TCP = -1/8;

spm_figure('GetWin','Cell-mediated immunity'); clf
y{4} = spm_immune_plot(P{4},spm_unvec(spm_vec(pC)/16,pC));

% Assume measurements every day for one month and add measurement noise
%--------------------------------------------------------------------------
U = 1:30;
U = U*24;

for i = 1:numel(y)
    Y{i} = (sqrt(y{i}(U,:)) + randn(size(y{i}(U,:)))/8).^2;
end

save('Y.mat','Y');
save('U.mat','U');
save('P.mat','P');

% Fit model to synthetic data
%==========================================================================
for i = 1:numel(Y)
    [F,Ep,Cp,pE,pC] = spm_immune(Y{i},U,pE,pC);
    DCM(i).F  = F;
    DCM(i).Ep = Ep;
    DCM(i).Cp = Cp;
    DCM(i).M.pE = pE;
    DCM(i).M.pC = pC;
end

% Use Bayesian model reduction to construct a confusion matrix
%--------------------------------------------------------------------------
for i = 1:numel(Y)
    for j = 1:numel(Y)
        L(i,j) = spm_log_evidence(DCM(j).Ep,DCM(j).Cp,pE,pC,P{i},spm_unvec(spm_vec(pC)/16,pC));
    end
end

% Convert to posterior probabilities and plot matrix
%--------------------------------------------------------------------------
p = spm_softmax(L);

spm_figure('GetWin','Confusion matrix'); clf
imagesc(p), axis square
ylabel('Model')
xlabel('Synthetic dataset')

clear DCM

% Design matrix
%--------------------------------------------------------------------------
X = [ones(32,1),rand(32,1)>0.5];

% Assume BCG increases T-cell response
%--------------------------------------------------------------------------
B = [pE.TCP 2];

for i = 1:size(X,1)
    R = pE;
    R.TCP = X(i,:)*B' + randn/16;
    y{i} = spm_immune_gen(R);
    y{i} = (sqrt(y{i}(U,:)) + randn(size(y{i}(U,:)))/8).^2;
    [F,Ep,Cp,pE,pC] = spm_immune(y{i},U,pE,pC);
    DCM{i,1}.F  = F;
    DCM{i,1}.Ep = Ep;
    DCM{i,1}.Cp = Cp;
    DCM{i,1}.M.pE = pE;
    DCM{i,1}.M.pC = pC;
end

GLM.X  = X;
GLM.Xn = {'const.','BCG'};
[PEB,DCM] = spm_dcm_peb(DCM,GLM,{'TCP'});

% Bayesian model averaging (over reduced models), testing for GLM effects
%--------------------------------------------------------------------------
[BMA,BMR] = spm_dcm_bmr_all(PEB,{'TCP'});
