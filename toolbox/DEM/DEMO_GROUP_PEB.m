function DEMO_GROUP_PEB
% Demonstration routine for empirical Bayes and Bayesian model reduction
%--------------------------------------------------------------------------
% This routine illustrates the use of Bayesian model reduction when
% inverting hierarchical (linear) models - it is essentially a software
% validation demo and proof of concept. It uses a parametric empirical
% Bayesian model (i.e., nested linear models) to eschew local minima issues
% and to assure the Laplace assumption is correct. In brief, the data are
% generated for multiple subjects, under a linear model with subject
% specific parameters at the first level and group specific parameters at
% the second. These model a group effect common to all subjects in a subset
% of parameters and differences in a further subset. In this demo, we
% consider the full hierarchical inversion of a multisubject  study by
% updating the priors at the first level, using the empirical priors from
% the second level. Crucially, this is done during the optimisation  at the
% first level (i.e., after  every iteration - or small number of iterations
% - at the first level.
%
% This provides a generic scheme for the hierarchical inversion of
% nonlinear and possibly dynamic models in which the first level
% optimisation is informed by the sufficient statistics of the second level
% (namely the empirical priors). This should be contrasted with the summary
% statistic approach, in which the second level optimisation, based upon
% the sufficient statistics of the first level (posteriors and priors) are
% computed after convergence of the first level. The results of inversion
% are compared in terms of the second level posteriors (and the second
% level free energy over iterations). Specifically, we compare a gold
% standard (PEB) inversion, with the summary statistic approach to
% empirical Bayes and the hierarchical inversion demonstrated in this
% routine.
% 
% The parameterisation of the models uses the format of DCM. This means
% parameters are specified as a structure with key parameters being in the
% fields A, B and C.
%
% See also: spm_dcm_bmr, spm_dcm_peb and spm_dcm_peb_bma
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_GROUP_PEB.m 6737 2016-03-03 12:05:51Z karl $


% set up
%==========================================================================
rng('default')

% model space - defined in terms of combinations of some parameters
%--------------------------------------------------------------------------
k     = spm_perm_mtx(3);
for i = 1:8;
    B{i} = k(i,:);
end

% model space
%--------------------------------------------------------------------------
mw  = 3;                              % true model (within)
mx  = 4;                              % true model (between)
Ns  = 16;                             % number of subjects
C   = 32;                             % within:between [co]variance ratio


% create subject-specifc GLM
%==========================================================================

% within subject effects:  condition specific effects 'B' (2 s.d.)
%--------------------------------------------------------------------------
pC          = 1/8;
sd          = sqrt(pC/C);
DCM.Ep.A    = randn(4,1)*sd;
DCM.Ep.B{1} = B{mw}*2*sd;
Np          = spm_length(DCM.Ep);
DCM.M.pE    = spm_zeros(DCM.Ep);
DCM.M.pC    = eye(Np,Np)*pC;

% between subject effects: constant and group difference
%--------------------------------------------------------------------------
X           = [ones(Ns,1) kron([-1;1],ones(Ns/2,1))];
DCM.Ex      = spm_zeros(DCM.Ep);
DCM.Ex.B{1} = B{mx}*2*sd;


% create subject-specifc DCM
%--------------------------------------------------------------------------
Ex    = spm_vec(DCM.Ex);
Ep    = spm_vec(DCM.Ep);
pC    = DCM.M.pC;
Cp    = sd*diag(~~spm_vec(Ep));
Ny    = 16;
for i = 1:Ns
    
    % generate data
    %----------------------------------------------------------------------
    Pp    = X(i,1)*Ep + X(i,2)*Ex + Cp*randn(Np,1);
    
    % generate data
    %----------------------------------------------------------------------
    Z{i,i} = randn(Ny,Np);
    y{i,1} = Z{i,i}*Pp + randn(Ny,1)/8;
      
    % design matrix and data
    %----------------------------------------------------------------------
    GCM{i,1}.xU    = Z{i,i};
    GCM{i,1}.xY.y  = y{i,1};
    GCM{i,1}.xY.X0 = [];
    
    % likelihood model and priors
    %----------------------------------------------------------------------
    GCM{i,1}.M.IS  = @(P,M,U) U*spm_vec(P);
    GCM{i,1}.M.pE  = DCM.M.pE;
    GCM{i,1}.M.pC  = pC;
    GCM{i,1}.Tp    = Pp;

end

% PEB (GLM) for inversion to provide a reference for BMR
%==========================================================================
Nx    = size(X,2);
Q     = spm_Ce(ones(1,Np));
for i = 1:Np
    Q{i} = kron(eye(Ns,Ns),Q{i})/128;
end
P{1}.X = spm_cat(Z);
P{1}.C = spm_Ce(ones(1,Ns)*Ny);
P{2}.X = kron(X,eye(Np,Np));
P{2}.C = Q;
P{3}.X = kron(zeros(Nx,1),zeros(Np,1));
P{3}.C = kron(eye(Nx,Nx),pC);

% Full hierarchical parametric empirical Bayes inversion
%--------------------------------------------------------------------------
[qP,~,F] = spm_PEB(spm_cat(y),P,1);

% record estimates as a reference
%--------------------------------------------------------------------------
PB.F  = F;
PB.Ep = qP{3}.E;
PB.Cp = qP{3}.C;

% repeat using non-linear empirical Bayes
%==========================================================================

% second level model
%--------------------------------------------------------------------------
M.X   = X;
M.pE  = DCM.M.pE;
M.pC  = DCM.M.pC;

% hierarchical inversion
%--------------------------------------------------------------------------
[gcm,peb,M] = spm_dcm_peb_fit(GCM,M);


% repeated using Bayesian model reduction summary statistic approach
%==========================================================================
GCM  = spm_dcm_fit(GCM);
PEB  = spm_dcm_peb(GCM,M);


% second level parameter estimates
%==========================================================================
spm_figure('GetWin','Figure 1'); clf

% estimated and true second level parameters
%--------------------------------------------------------------------------
subplot(2,2,1), spm_plot_ci(PB.Ep,PB.Cp), hold on, bar([Ep;Ex],1/2), hold off
xlabel('parameters'), ylabel('expectation'), 
title('Parametric Bayes','FontSize',16), axis square

subplot(2,2,2), spm_plot_ci(PEB.Ep(:),PEB.Cp), hold on, bar([Ep;Ex],1/2), hold off
xlabel('parameters'), ylabel('expectation'), 
title('Sufficient statistics','FontSize',16), axis square

subplot(2,2,3), spm_plot_ci(peb.Ep(:),peb.Cp), hold on, bar([Ep;Ex],1/2), hold off
xlabel('parameters'), ylabel('expectation'), 
title('Hierarchical inversion','FontSize',16), axis square

subplot(2,2,4), bar(gcm{1}.FEB)
xlabel('iterations'), ylabel('free energy'), 
title('second level free energy','FontSize',16), axis square

return


% Bayesian model reduction with and without hierarchical inversion
%==========================================================================

% define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),spm_length(DCM.Ep));
k     = spm_fieldindices(DCM.M.pE,'B');
for i = 1:length(B)
    K(i,k) = spm_vec(B{i})';
end

% defined model in terms of prior covariance
%--------------------------------------------------------------------------
Nm    = size(K,2);
for i = 1:Ns
    for j = 1:Nm
        gcm{i,j}      = gcm{i,1};
        GCM{i,j}      = GCM{i,1};
        gcm{i,j}.M.pC = diag(K(j,:))*M.pC*diag(K(j,:));
        GCM{i,j}.M.pC = diag(K(j,:))*M.pC*diag(K(j,:));
    end
end

rcm   = spm_dcm_bmr(gcm);
RCM   = spm_dcm_bmr(GCM);


% Free energies
%--------------------------------------------------------------------------
for i = 1:Ns
    for j = 1:Nm
        G(i,j,1) = rcm{i,j}.F - rcm{i,1}.F;
        G(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
    end
end

%  free energy model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');clf

p  = spm_softmax(sum(G(:,:,1))'); [m i] = max(p); 
subplot(2,2,1), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Iterative PEB','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p  = spm_softmax(sum(G(:,:,1))'); [m i] = max(p); 
subplot(2,2,2), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Bayesian model reduction','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square


