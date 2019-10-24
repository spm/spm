function DEMO_Lindley_paradox(pC,hE,hC)
% FORMAT DEMO_BAYES_FACTORS(pC,hE,hC)
% Demonstration Bayes factors and classical p-values
%--------------------------------------------------------------------------
% pC   - prior covariance             (e.g., 4)
% hE   - expectation of log precision (e.g., 1)
% hC   - covariance of log precision  (e.g., 1/8)
%
% This demonstration routine uses a simple linear model to examine the
% relationship between free energy differences or log Bayes factors and
% classical F statistics.
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_Lindley_paradox.m 7679 2019-10-24 15:54:07Z spm $


% set up
%--------------------------------------------------------------------------
% rng('default')

try, pC; catch, pC = 1/4;  end
try, hE; catch, hE = 0;    end
try, hC; catch, hC = 1/8;  end

sigma_a = 1/8;               % weak effect size - alternative hypothesis
sigma_r = 1/64;              % reduced (null) effect size
epsilon = 1/16;

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
Ns   = 256;                       % number of samples
pE   = M.pE;                      % full prior expectations
pC   = M.pC;                      % full prior covariance
rC   = pC;                        % restricted or reduced priors
rC(1) = sigma_r^2;


k     = (1:Ns) < Ns/2;            % null and alternative (prevalence)
N     = 2.^(3:10);                % umber of subjects
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
        beta      = [rand(1)*b; 0];
        Y         = X*beta + randn(N(n),1);
        
        % Bayesian analysis (full comparison and model reduction)
        %------------------------------------------------------------------
        [qE,qC]   =  spm_nlsi_GN(M,X,Y);
        F(i,1)    = -spm_log_evidence(qE,qC,pE,pC,pE,rC);
        QE(i,1)   = qE(1);
        QC(i,1)   = qC(1);
        
        % classical analysis
        %------------------------------------------------------------------
        [t,df,qE] = spm_ancova(X,[],Y,[1 0;0 0]);
        T(i,1)    = t;
        qe(i,1)   = qE(1);
        P(i,1)    = beta(1);
        
    end
    

    % classical threshold and PPV
    %----------------------------------------------------------------------
    u      = spm_invFcdf(0.95,df);
    i      = find(P >  0);
    j      = find(P <= 0);
    
    FPR(n) = sum(T(j) > u)/length(j);
    PPV(n) = sum(T(i) > u)/sum(T > u);
    
    i      = find(P >  epsilon);
    j      = find(P <= epsilon);
    
    EPR(n) = sum(T(j) > u)/length(j);
    EPV(n) = sum(T(i) > u)/sum(T > u);
    
    fpr(n) = sum(F(j) > 0)/length(j);
    ppv(n) = sum(F(i) > 0)/sum(F > 0);
   
    i      = find(T > u);
    Ep(n)  = mean(P(i));
    Cp(n)  = var(P(i));
    i      = find(F > 0);
    ep(n)  = mean(P(i));
    cp(n)  = var(P(i));

    
end




% show results
%==========================================================================
spm_figure('GetWin','Graphics_null');clf

subplot(2,2,1)
plot(N,(N*0 + 0.05),'--r'),    hold on
plot(N,(N*0 + 0.80),'--b'),    hold on
plot(N,FPR, 'r',N,fpr,'--r'),  hold on
plot(N,EPR,'g',N,EPV,'g'), hold on
plot(N,PPV, 'b',N,ppv,'--b'), hold off
xlabel('Number of samples'), ylabel('Probability')
title('PPV and FPR','FontSize',16)
axis([0 N(end) 0 1]); axis square

subplot(2,2,3)
spm_plot_ci(Ep,Cp,N),          hold on
plot(N,N - N + sigma_r,'--'),  hold on
plot(N,N - N + sigma_a,'--r'), hold off
xlabel('Number of samples'), ylabel('effect science')
title('Detected effect size','FontSize',16)
axis square, axis([0 N(end) 0 1]); axis square

subplot(2,2,4)
spm_plot_ci(ep,cp,N),          hold on
plot(N,N - N + sigma_r,'--'),  hold on
plot(N,N - N + sigma_a,'--r'), hold off
xlabel('Number of samples'), ylabel('effect science')
title('Bayesian','FontSize',16)
axis square, axis([0 N(end) 0 1]); axis square


% (linear) mapping between free energy difference and F ratio
%--------------------------------------------------------------------------
j   = abs(F) < 32;
b   = pinv([F(j) ones(size(F(j)))])*T(j);
Fq  = (-32:32)';
Tq  = [Fq, ones(size(Fq))]*b;

subplot(2,2,2)
plot(F,T,'.b','Markersize',8), hold on
plot(Fq,Tq,'b'), hold on
plot([3 3],[0 16],':r'),  hold on
plot([0 0],[0 16],'--r'), hold on
plot([-32, 32],[u u],':k'), hold off
xlabel('Free energy difference'), ylabel('Classical F-ratio')
title('Classical and Bayesian statistics','FontSize',16)
axis([-8 8 0 16])
axis square

