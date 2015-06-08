function DEMO_Lindley_paradox(pC,hE,hC,N,b)
% FORMAT DEMO_BAYES_FACTORS(pC,hE,hC,N)
% Demonstration Bayes factors and classical p-values
%--------------------------------------------------------------------------
% pC   - prior covariance             (e.g., 4)
% hE   - expectation of log precision (e.g., 1)
% hC   - covariance of log precision  (e.g., 1/8)
% N    - number of observations       (e.g., 16)
% b    - relative variance under alternate and null (e.g., 1/32)
%
% This demonstration routine uses a simple linear model to examine the
% relationship between free energy differences or log Bayes factors and
% classical F statistics. Using re-randomisation of a design matrix, it
% computes the null distribution over both statistics and plots them
% against each other.  There is a linear relationship, which allows one to
% evaluate the false-positive rate for any threshold on the Bayes factor.
% Ideally, one would like to see a positive log Bayes factor map to a
% classical threshold of p=0.05. The offset and slope of the linear
% relationship between the two statistics depends upon prior beliefs about
% the covariance of the parameters and the log precision. These can be
% changed by editing the code below (or supplying input arguments).
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_Lindley_paradox.m 6476 2015-06-08 09:37:01Z karl $


% set up
%--------------------------------------------------------------------------
rng('default')

try, pC; catch, pC = 1/8;    end
try, hE; catch, hE = 0;    end
try, hC; catch, hC = 1/8;    end

sigma_a = 1/8;               % weak effect size – alternative hypothesis
sigma_r = 1/16;              % reduced (null) effect size

% Model specification
%==========================================================================
M.nograph = 1;
M.noprint = 1;

M.IS = @(P,M,U) U*P;
M.pE = [0; 0];
M.pC = eye(2,2)*pC;
M.hE = hE;
M.hC = hC;

% re-randomisation
%--------------------------------------------------------------------------
Ns   = 1e3;
pE   = M.pE;                    % full prior expectations
pC   = M.pC;                    % full prior covariance
rC   = pC;                      % restricted or reduced priors
rC(1,1) = sigma_r^2;


k     = kron([1;0],ones(Ns/2,1)); % null and alterantive
N     = 8:8:128;                % umber of subjects
for n = 1:length(N)
    
    % design matrix and contrast
    %--------------------------------------------------------------------------
    X     = [randn(N(n),1) ones(N(n),1)];
    for i = 1:Ns
        
        % effect size - alternative or null
        %------------------------------------------------------------------
        if k(i); b = 0; else, b = sigma_a; end
        
        % generate data
        %------------------------------------------------------------------
        beta      = [b + randn(1)*sigma_r; 0];
        Y         = X*beta + randn(N(n),1);
        
        % Bayesian analysis (full comparison and model reduction)
        %------------------------------------------------------------------
        [qE,qC]   =  spm_nlsi_GN(M,X,Y);
        F(i,1)    = -spm_log_evidence(qE,qC,pE,pC,pE,rC);
        QE(i,1)   = qE(1);
        QC(i,1)   = qC(1);
        
        % classical analysis
        %------------------------------------------------------------------
        [t,df,qE] = spm_ancova(X,[],Y,[1; 0]);
        T(i,1)    = t;
        qe(i,1)   = qE(1);
        pe(i,1)   = beta(1);
        
    end
    
    
    % (linear) mapping between free energy difference and F ratio
    %--------------------------------------------------------------------------
    T   = T.^2;
    u   = u^2;
    j   = abs(F) < 32;
    b   = pinv([F(j) ones(size(F(j)))])*T(j);
    Fq  = (-32:32)';
    Tq  = [Fq, ones(size(Fq))]*b;
    
    
    
    subplot(2,2,2)
    plot(F,T,'.b','Markersize',8), hold on
    plot(Fq,Tq,'b'), hold on
    plot([3 3],[0 16],':r'), hold on
    plot([0 0],[0 16],'--r'), hold on
    plot([-32, 32],[u u],':k'), hold off
    xlabel('Free energy difference'), ylabel('Classical F-ratio')
    title('Classical and Bayesian statistics','FontSize',16)
    axis([-8 8 0 16])
    axis square


    % classical threshold and PPV
    %----------------------------------------------------------------------
    u      = spm_invTcdf(0.95,df(2));
    i      = find(~k);
    j      = find( k);
    
    FPR(n) = sum(T(j) > u)/length(j);
    PPV(n) = sum(T(i) > u)/sum(T > u);
    
    fpr(n) = sum(F(j) > 0)/length(j);
    ppv(n) = sum(F(i) > 0)/sum(F > 0);
    
    
    i      = find(T > u);
    ep(n)  = mean(pe(i));
    cp(n)  = var(pe(i));
    i      = find(F > 0);
    Ep(n)  = mean(pe(i));
    Cp(n)  = var(pe(i));

    
end




% show results
%==========================================================================
spm_figure('GetWin','Graphics_null');clf

subplot(2,2,1)
plot(N,(N*0 + 0.05),'--r'),  hold on
plot(N,(N*0 + 0.80),'--b'),  hold on
plot(N,FPR,'r',N,fpr,'--r'), hold on
plot(N,PPV,'b',N,ppv,'--b'), hold off
xlabel('Number of samples'), ylabel('Probability')
title('PPV and FPR','FontSize',16)
axis([0 N(end) 0 1]); axis square

subplot(2,2,3)
spm_plot_ci(ep,cp,N),          hold on
plot(N,N - N + sigma_r,'--'),  hold on
plot(N,N - N + sigma_a,'--r'), hold off
xlabel('Number of samples'), ylabel('effect science')
title('Detected effect size','FontSize',16)
axis square, axis([0 N(end) 0 1]); axis square

subplot(2,2,4)
spm_plot_ci(Ep,Cp,N),          hold on
plot(N,N - N + sigma_r,'--'),  hold on
plot(N,N - N + sigma_a,'--r'), hold off
xlabel('Number of samples'), ylabel('effect science')
title('Bayesian','FontSize',16)
axis square, axis([0 N(end) 0 1]); axis square




