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
% This demonstration routine illustrates a canonical covariates analysis in
% which hypotheses are specified in terms of second-order matrices (of the
% sort used in representational similarity analysis). This routine
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
% of one(corresponding to a key contrast). However, the code has been
% written to handle arbitrary hypothesis matrices (corresponding to F
% contrasts) that test a subspace of two or more dimensions.
%
% To the extent that this reproduces the hypothesis testing of
% representational similarity analysis, there is an important observation:
% this analysis works for a single voxel. In other words, representational
% similarity analysis is not an inherently multivariate approach.
%
% This routine deliberately mixes two (main) effects in equal measure,
% within the same region of interest. This is to highlight the
% inappropriate application of hypothesis selection; here demonstrated via
% Bayesian model comparison using the Bayesian information criteria. In
% other words, several hypotheses about a particular region could be true
% at the same time.
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
% $Id: DEMO_CVA_RSA.m 7303 2018-04-28 16:00:22Z karl $


% preliminaries
%--------------------------------------------------------------------------
Nn   = 8;       % number of subjects
Nv   = 32;      % number of voxels in the volume of interest
Ns   = 16;      % number of stimuli objects
Np   = 24;      % number of presentations (per subject)


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
    title('contrast'), xlabel('stimulus'), box off
    
    subplot(6,3,i + 3), imagesc(1 - C{i}')
    title('hypothesis matrix'), xlabel('stimulus'), ylabel('stimulus')
    box off, axis square
end

% synthetic data
%==========================================================================
% Let us now creates an synthetic data from a region that encodes the main
% effect of the parametric attribute and the context (but not the
% interaction). This functional specialisation can be specified in terms of
% canonical values as follows (assuming each stimulus is presented Np times
% to Ns subjects):
%--------------------------------------------------------------------------
CV    = [1 1 0];
X     = kron(ones(Np,1),eye(Ns,Ns));
X0    = ones(size(X,1),1);
for i = 1:Nn
    
    % functionally specialised responses, randomly distributed over voxels
    %----------------------------------------------------------------------
    B    = c*diag(CV)*randn(Nc,Nv);
    
    % observation error
    %----------------------------------------------------------------------
    e    = randn(size(X,1),Nv)/32;
    
    % known confounds
    %----------------------------------------------------------------------
    B0   = randn(size(X0,2),Nv)/32;
    
    % known confounds
    %----------------------------------------------------------------------
    Y{i} = X*B + X0*B0 + e;
    
end

% canonical variates analysis (i.e. representational similarity analysis)
%==========================================================================
% We can now recover the canonical effects using CVA and accumulate the
% evidence for different contrasts or hypothesis matrices over subjects.
% Crucially, this analysis can either be specified directly in terms of the
% first-order contrasts – or the second order contrast matrices; i.e., the
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
    
    %  and accumulate the log evidence for each contrast
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
xlabel('stimulus'), ylabel('stimulus'), box off, axis square

% as inference about each contrast
%--------------------------------------------------------------------------
subplot(3,3,7), bar(RSA.F), xlabel('contrast'), ylabel('log evidence')
title('subject specific effects'), axis square

subplot(3,3,8), bar(sum(RSA.F,2)), xlabel('contrast'), ylabel('log evidence')
title('pooled'), axis square

% and in terms of a model comparison (i.e., the best hypothesis)
%--------------------------------------------------------------------------
subplot(3,3,9), bar(spm_softmax(sum(RSA.F,2))), xlabel('contrast'), 
ylabel('posterior probability'), title('model comparison'), axis square





