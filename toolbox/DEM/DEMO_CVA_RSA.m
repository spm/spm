function RSA = DEMO_CVA_RSA
% Canonical Variate Analysis and representational similarity analysis
% FORMAT RSA = DEMO_CVA_RSA
%
% output structure
%--------------------------------------------------------------------------
% RSA.C    - hypotheses matrices
% RSA.c    - (orthogonalised) contrasts
% RSA.W    - (second-order) canonical weights
% RSA.X    - design matrix
% RSA.Y    - data
% RSA.X0   - confounds
% RSA.F    - (BIC) log evidence
%__________________________________________________________________________
%
% This demonstration routine starts with a canonical covariates analysis in
% which hypotheses are specified in terms of second-order matrices (of the
% sort used in representational similarity analysis). This part
% illustrates the inversion of a multivariate linear model over multiple
% subjects, testing for the expression of multivariate responses under each
% of three hypotheses. Furthermore, it illustrates the (Bayesian) model
% comparison under the assumption that only one hypothesisis true.
%
% The three hypotheses correspond to a main effect of a parametric variable
% (e.g., the degree to which something is judged valuable), the main
% effect of a categorical variable (e.g., big or small) and their
% interaction. Note that this requires a specification in terms of
% second-order hypothesis matrices that are not expressed in terms of
% similarities per se. In other words, the second-order hypotheses are
% assumed to be in the form of covariance matrices; as opposed to
% correlation matrices.
%
% This routine demonstrates the testing of hypothesis matrices with a rank
% of one (corresponding to a T-contrast). However, the code has been
% written to handle arbitrary hypothesis matrices (corresponding to F-
% contrasts) that test a subspace of two or more dimensions.
%
% To the extent that this reproduces the hypothesis testing of
% representational similarity analysis, there is an important observation:
% this analysis works for a single voxel. In other words, representational
% similarity analysis is not an inherently multivariate approach.
%
% This illustration deliberately mixes two (main) effects in equal measure,
% within the same region of interest. This is to highlight the
% inappropriate application of hypothesis selection; here demonstrated via
% Bayesian model comparison using the Bayesian information criteria. In
% other words, several hypotheses about a particular region could be true
% at the same time.
%
% We then revisit exactly the same problem (i.e., Bayesian model comparison
% of covariance components of second-order responses) using variational
% Laplace to estimate the contributions of each component of pattern
% explicitly. This has the advantage of enabling parametric empirical Bayes
% at the between subject level - and subsequent Bayesian model reduction.
%
% References:
%
% Characterizing dynamic brain responses with fMRI: a multivariate
% approach. Friston KJ, Frith CD, Frackowiak RS, Turner R. NeuroImage. 1995
% Jun;2(2):166-72.
%
% A multivariate analysis of evoked responses in EEG and MEG data. Friston
% KJ, Stephan KM, Heather JD, Frith CD, Ioannides AA, Liu LC, Rugg MD,
% Vieth J, Keber H, Hunter K, Frackowiak RS. NeuroImage. 1996 Jun;
% 3(3):167-174.
%
% Population level inference for multivariate MEG analysis. Jafarpour A,
% Barnes G, Fuentemilla Lluis, Duzel E, Penny WD. PLoS One. 2013.
% 8(8): e71305
%__________________________________________________________________________
% Copyright (C) 2006-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEMO_CVA_RSA.m 7679 2019-10-24 15:54:07Z spm $


% preliminaries
%--------------------------------------------------------------------------
rng('default')

Nn   = 8;       % number of subjects
Nv   = 32;      % number of voxels in the volume of interest
Ns   = 16;      % number of stimuli (i.e., objects)
Np   = 24;      % number of presentations (per subject)

CV   = [1 1 0]; % true effects (main effects and interaction)
s    = 1/2;     % standard deviation of noise
k    = 1/8;     % spatial correlations of noise

% special convolution kernel
%--------------------------------------------------------------------------
K    = toeplitz(exp(-(0:(Nv - 1)).^2/k^2/2));
K    = K/mean(diag(K*K'));

% canonical contrasts
%==========================================================================
% Imagine Ns objects that have been designed or rated along two attributes,
% say a parametric attribute (e.g., brightness) and a categorical attribute
% that, here, stands in for context (e.g., big or small). we will now
% categorise or classify each object along both attributes and evaluate
% the interaction:
%--------------------------------------------------------------------------
c(:,1) = spm_detrend(randn(Ns,1));          % parametric attribute
c(:,2) = kron([1; - 1],ones(Ns/2,1));       % categorical attribute
c(:,1) = c(:,1) - c(:,2)*(c(:,2)\c(:,1));   % orthogonalise
c(:,3) = c(:,1).*c(:,2);                    % interaction
c      = spm_orth(c,'norm');                % orthonormalise
Nc     = size(c,2);                         % number of contrasts

%--------------------------------------------------------------------------
% These canonical effects can be expressed in terms of contrasts or in
% terms of their outer products that correspond to a second-order contrast
% or a hypothesis matrix
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf
for i = 1:Nc
    
    % second-order contrast (i.e., similarity) matrix
    %----------------------------------------------------------------------
    C{i} = c(:,i)*c(:,i)';
    
    subplot(6,3,i), bar(c(:,i))
    title('Contrast'), xlabel('Stimulus'), box off, spm_axis tight
    
    subplot(6,3,i + 3), imagesc(1 - C{i}')
    title('Hypothesis matrix'), xlabel('Stimulus'), ylabel('Stimulus')
    box off, axis square
end

% synthetic data
%==========================================================================
% Let us now create synthetic data from a region that encodes the main
% effect of the parametric attribute and the context (but not the
% interaction). This functional specialisation can be specified in terms of
% canonical values as follows (assuming each stimulus is presented Np times
% to Ns subjects):
%--------------------------------------------------------------------------
X     = kron(ones(Np,1),eye(Ns,Ns));
X0    = ones(size(X,1),1);
for i = 1:Nn
    
    % functionally specialised responses, randomly distributed over voxels
    %----------------------------------------------------------------------
    B    = c*diag(CV)*randn(Nc,Nv)*K;
    
    % observation error
    %----------------------------------------------------------------------
    e    = randn(size(X,1),Nv)*K*s;
    
    % known confounds
    %----------------------------------------------------------------------
    B0   = randn(size(X0,2),Nv)/16;
    
    % response variable
    %----------------------------------------------------------------------
    Y{i} = X*B + X0*B0 + e;
    
end

% canonical variates analysis (i.e. representational similarity analysis)
%==========================================================================
% We can now recover the canonical effects using CVA and accumulate the
% evidence for different contrasts or hypothesis matrices over subjects.
% Crucially, this analysis can either be specified directly in terms of the
% first-order contrasts - or the second order contrast matrices; i.e., the
% hypothesis matrices. Here the implicit contrasts are recovered using SVD:
%--------------------------------------------------------------------------

% get contrasts from hypothesis matrices
%--------------------------------------------------------------------------
clear c
for i = 1:Nc
    c{i} = full(spm_svd(C{i}));
    
    % ensure (multivariate) contrasts are orthogonal
    %----------------------------------------------------------------------
    if i > 1
        c{i}  = c{i} - cc*(cc\c{i});
    end
    
    % accumulate contrast space
    %----------------------------------------------------------------------
    cc   = full(spm_cat(c));
end

% canonical variates analysis
%--------------------------------------------------------------------------
for i = 1:Nn
    
    % canonical variance analysis (all contrasts)
    %----------------------------------------------------------------------
    CVA       = spm_cva(Y{i},X,X0,cc);
    
    % accumulate canonical vectors (second order statistics)
    %----------------------------------------------------------------------
    W(:,:,i)  = cc*CVA.W*CVA.W'*cc';
    
    % and accumulate the log evidence for each contrast
    %----------------------------------------------------------------------
    for j = 1:Nc
        CVA    = spm_cva(Y{i},X,X0,c{j});
        F(j,i) = CVA.bic;
    end
    
end

% output structure
%--------------------------------------------------------------------------
RSA.C  = C;                   % hypotheses matrices
RSA.c  = c;                   % (orthogonalised) contrasts
RSA.W  = W;                   % (second-order) canonical weights
RSA.X  = X;                   % design matrix
RSA.Y  = Y;                   % data
RSA.X0 = X0;                  % confounds
RSA.F  = F;                   % (BIC) log evidence


% Results
%==========================================================================
% now report results in terms of the average (second-order matrix of)
% canonical vectors
%--------------------------------------------------------------------------
subplot(3,1,2), imagesc(1 - sum(RSA.W,3))
title('Canonical similarity matrix','FontSize',16)
xlabel('Stimulus'), ylabel('Stimulus'), box off, axis square

% as inference about each contrast
%--------------------------------------------------------------------------
subplot(3,3,7), bar(RSA.F), xlabel('contrast'), ylabel('log evidence')
title('Subject-specific effects'), axis square

subplot(3,3,8), bar(sum(RSA.F,2)), xlabel('Contrast'), ylabel('Log-evidence')
title('Pooled'), axis square

% and in terms of a model comparison (i.e., the best hypothesis)
%--------------------------------------------------------------------------
subplot(3,3,9), bar(spm_softmax(sum(RSA.F,2))), xlabel('Contrast'),
ylabel('Posterior probability'), title('Model comparison'), axis square


% now start again using variational Laplace
%--------------------------------------------------------------------------

% variational Bayes (i.e. representational similarity analysis)
%==========================================================================
% We will now repeat the analysis using a covariance component approach
% (c.f., pattern component modelling).
%--------------------------------------------------------------------------

% get contrasts from hypothesis matrices
%--------------------------------------------------------------------------
clear RSA
for i = 1:Nc
    c{i} = full(spm_svd(C{i}));
end

% create covariance components (in the space of pinv(X))
%--------------------------------------------------------------------------
for i = 1:Nc
    Q{i}  = C{i};                                  % components of interest
end
iX        = pinv(X);                               % projector
Q{Nc + 1} = iX*iX';                                % i.i.d. error
R         = eye(size(X,1)) - [X X0]*pinv([X X0]);  % residual projector
X0        = iX*X0;                                 % projected confounds

% component analysis
%--------------------------------------------------------------------------
for i = 1:Nn
    
    % estimate spatial degrees of freedom (Nv)
    %----------------------------------------------------------------------
    e  = R*Y{i};
    e  = e'*e;
    Nv = trace(e)^2/trace(e*e);
    
    
    % second-order response
    %----------------------------------------------------------------------
    YY = iX*Y{i};
    YY = YY*YY';
    
    % variational covariance components analysis
    %----------------------------------------------------------------------
    [Cy,h,Ph,F,Fa,Fc,Eh,Ch,hE,hC,Qh] = spm_reml_sc(YY,X0,Q,Nv,-16,128);
    
    % accumulate
    %----------------------------------------------------------------------
    RSA{i}.M.pE = hE;         % prior expectation of hyperparameters
    RSA{i}.M.pC = hC;         % prior covariances of hyperparameters
    RSA{i}.Ep   = Eh;         % posterior expectations
    RSA{i}.Cp   = Ch;         % posterior covariance
    RSA{i}.Q    = Qh;         % scaled covariance components
    RSA{i}.F    = F;          % free energy
    
end


% parametric empirical Bayes
%==========================================================================
% now use parametric empirical Bayes to compute group average covariance
% hyperparameters
%--------------------------------------------------------------------------
M.X       = ones(Nn,1);
[PEB,RSA] = spm_dcm_peb(RSA(:),M);

% recompute second-order response
%--------------------------------------------------------------------------
G     = 0;
for i = 1:Nc
    G = G + exp(PEB.Ep(i))*RSA{1}.Q{i};
end

% and use Bayesian model reduction for model comparison
%==========================================================================
pE    = PEB.M.pE;
pC    = PEB.M.pC;
qE    = PEB.Ep;
qC    = PEB.Cp;
for i = 1:Nc
    
    % Place precise shrinkage priors on each component
    %----------------------------------------------------------------------
    rC      = pC;
    rC(i,i) = 1/128;
    Fc(i,1) = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    
    
    % and assess evidence for just this component
    %----------------------------------------------------------------------
    rC      = pC;
    j       = 1:Nc; j(i) = [];
    rC(j,j) = 1/128;
    Fs(i,1) = spm_log_evidence(qE,qC,pE,pC,pE,rC);
    
end

PEB.Fc = max(Fc) - Fc;
PEB.Fs = Fs - min(Fs);


% Results
%==========================================================================
spm_figure('GetWin','RSA');clf
for i = 1:Nc
    
    % first-order contrast (i.e., similarity) matrix
    %----------------------------------------------------------------------
    subplot(6,3,i), bar(c{i})
    title('Contrast'), xlabel('Stimulus'), spm_axis tight, axis off
    
    % second-order contrast (i.e., similarity) matrix
    %----------------------------------------------------------------------
    subplot(6,3,i + 3), imagesc(1 - C{i}')
    title('Hypothesis matrix'), xlabel('Stimulus'), ylabel('Stimulus')
    box off, axis square
end

% now report results in terms of average similarity matrix
%--------------------------------------------------------------------------
subplot(3,1,2), imagesc(1 - G)
title('Variational similarity matrix','FontSize',16)
xlabel('Stimulus'), ylabel('Stimulus'), box off, axis square

% and subject specific posterior over hyperparameters
%--------------------------------------------------------------------------
for i = 1:Nn
    EP(:,i) = RSA{i}.Ep;
    CP(:,i) = diag(RSA{i}.Cp);
end

subplot(3,3,7), spm_plot_ci(EP',CP',[],[],'exp'), xlabel('Component')
ylabel('Contribution'), title('Subject-specific effects','Fontsize',14)
axis square tight, set(gca,'YLim',[0 6])

% show the results of Bayesian model comparison
%--------------------------------------------------------------------------
subplot(3,3,8), bar(PEB.Fc),   hold on
plot([0,(Nc + 1)],[3,3],':r'), hold off
xlabel('Component'), ylabel('Log-evidence')
title('Model comparison','Fontsize',14), axis square

% and selection (i.e., the best hypothesis)
%--------------------------------------------------------------------------
subplot(3,3,9), bar(spm_softmax(PEB.Fs(:))), xlabel('component'),
ylabel('Posterior probability'), title('Model selection','Fontsize',14)
axis square





