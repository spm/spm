function DEMO_BMR_PEB
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
% of parameters and differences in a further subset. The objective of
% empirical Bayesian inversion is to recover the group effects in terms of
% posteriors and perform Bayesian model comparison at the second (between
% subject) level.
%
% This provides empirical shrinkage priors at the first level, which can be
% used to compute the predictive posterior for any subject. In turn, the
% predictive posterior can be used for leave-one-out cross validation.
%
% The key aspect of this approach to empirical Bayesian modelling is that
% we use Bayesian model reduction throughout. In other words, after the
% subject-specific models have been inverted the data are discarded and we
% deal only with the free energies and posteriors for subsequent
% hierarchical analysis. This can be computationally very efficient when
% dealing with large first-level or complicated models: as in DCM.
% 
% The parameterisation of the models uses the format of DCM. This means
% parameters are specified as a structure with key parameters being in the
% fields A, B and C.
%
% See also: spm_dcm_bmr, spm_dcm_peb and spm_dcm_peb_bma
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_BMR_PEB.m 6737 2016-03-03 12:05:51Z karl $


% set up
%==========================================================================
corr = @(x,y) subsref(corrcoef(x,y),substruct('()',{1,2})); % Stats tbx
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
Nm  = length(B);                      % number of models
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


% (RFX) BMA - define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),spm_length(DCM.Ep));
k     = spm_fieldindices(DCM.M.pE,'B');
for i = 1:length(B)
    K(i,k) = spm_vec(B{i})';
end


% create subject-specifc DCM
%--------------------------------------------------------------------------
Ex    = spm_vec(DCM.Ex);
Ep    = spm_vec(DCM.Ep);
pC    = DCM.M.pC;
Cp    = sd*diag(~~spm_vec(Ep));
Ny    = 16;
for i = 1:Ns
    
    % report
    %----------------------------------------------------------------------
    fprintf('Creating subject %i\n',i)
    
    
    % generate data
    %----------------------------------------------------------------------
    Pp    = X(i,1)*Ep + X(i,2)*Ex + Cp*randn(Np,1);

    % generate data
    %----------------------------------------------------------------------
    Z{i,i} = randn(Ny,Np);
    y{i,1} = Z{i,i}*Pp + randn(Ny,1)/8;

    % invert models
    %----------------------------------------------------------------------
    P{1}.X = Z{i,i};
    P{1}.C = {eye(Ny,Ny)};
    P{2}.X = zeros(Np,1);
    P{2}.C = pC;
    for j = 1:Nm
        
        
        % defined model in terms of prior covariance
        %------------------------------------------------------------------
        P{2}.C        = diag(K(j,:))*pC*diag(K(j,:));
        [qP,~,F]      = spm_PEB(y{i,1},P,1);
        
        % store results in group array
        %------------------------------------------------------------------
        GCM{i,j}.M.pE = DCM.M.pE;
        GCM{i,j}.M.pC = P{2}.C;
        GCM{i,j}.Ep   = qP{2}.E;
        GCM{i,j}.Cp   = qP{2}.C;
        GCM{i,j}.B    = B(j);
        GCM{i,j}.F    = F;
        GCM{i,j}.Tp   = Pp;
        
    end
end

% PEB (GLM) for inversion to provide a reference for BMR
%==========================================================================
clear P

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
for i = 1:Nm
    for j = 1:Nm
        
        % model in terms of constant and differences
        %------------------------------------------------------------------
        pCi        = diag(K(i,:))*pC*diag(K(i,:));
        pCj        = diag(K(j,:))*pC*diag(K(j,:));
        P{3}.C     = blkdiag(pCi,pCj);
        [qP,~,F]   = spm_PEB(spm_cat(y),P,1);
        
        
        % record Estimates as a reference
        %------------------------------------------------------------------
        PB(i,j).F  = F;
        PB(i,j).Ep = qP{3}.E;
        PB(i,j).Cp = qP{3}.C;
        
        % free energy
        %------------------------------------------------------------------
        PF(i,j)    = F;
    end
end




% Bayesian model reduction - for each subject
%==========================================================================
RCM           = spm_dcm_bmr(GCM);
                    
% hierarchical (empirical Bayes) analysis using model reduction
%==========================================================================
[REB,PCM]     = spm_dcm_peb(RCM);

% BMA - (first level)
%--------------------------------------------------------------------------
bma   = spm_dcm_bma(GCM);
rma   = spm_dcm_bma(RCM);
pma   = spm_dcm_bma(PCM);

% second level model
%--------------------------------------------------------------------------
M     = struct('X',X);

% BMC - (first and second level) (with optimisation of hyperprior)
%--------------------------------------------------------------------------
[BMC,M] = spm_dcm_peb_test(RCM(:,1),M,{'B'});

% BMA - (second level)
%--------------------------------------------------------------------------
PEB     = spm_dcm_peb(RCM(:,1),M);
BMA     = spm_dcm_peb_bmc(PEB,RCM(1,:));


% posterior predictive density and cross validation
%==========================================================================
spm_dcm_loo(RCM(:,1),X,{'B'});


% show results
%==========================================================================
clear Q
for i = 1:Ns
    
    % Parameter averages
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Tp);
    Q(:,i,2) = spm_vec(bma.SUB(i).Ep);
    Q(:,i,3) = spm_vec(rma.SUB(i).Ep);
    Q(:,i,4) = spm_vec(pma.SUB(i).Ep);
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,1) = GCM{i,j}.F - GCM{i,1}.F;
        F(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
        F(i,j,3) = REB(j).F   - REB(1).F;
    end
    
end

% indices to select parameters
%--------------------------------------------------------------------------
iA    = spm_fieldindices(DCM.M.pE,'A');
iB    = spm_fieldindices(DCM.M.pE,'B');


% plot results: Bayesian model reduction vs. reduced models
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');clf

f  = F(:,:,1); f = f - max(f(:)) + 64; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (FFX)','FontSize',16)
axis square

f  = sum(f,1); f  = f - max(f) + 64; f(f < 0) = 0;
subplot(3,2,3), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

p  = exp(f - max(f)); p = p/sum(p); [m,i] = max(p); 
subplot(3,2,5), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f  = F(:,:,2); f = f - max(f(:)) + 64; f(f < 0) = 0;
subplot(3,2,2), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

f  = sum(f,1); f  = f - max(f) + 64; f(f < 0) = 0;
subplot(3,2,4), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square

p  = exp(f - max(f)); p = p/sum(p); [m,i] = max(p); 
subplot(3,2,6), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square



% first level parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');clf, ALim = 1/2;

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,2)));
str = sprintf('BMA: correlation = %-0.3f',r);
subplot(3,2,1), plot(Q(iA,:,1),Q(iA,:,2),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,2),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,3)));
str = sprintf('BMR: correlation = %-0.3f',r);
subplot(3,2,3), plot(Q(iA,:,1),Q(iA,:,3),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,3),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,4)));
str = sprintf('PEB: correlation = %-0.3f',r);
subplot(3,2,5), plot(Q(iA,:,1),Q(iA,:,4),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,4),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

f   = sum(F(:,:,1)); f = f - max(f(:)); f(f < -64) = -64;
p   = exp(f - max(f)); p = p/sum(p);
subplot(3,2,2), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f   = sum(F(:,:,2)); f = f - max(f(:)); f(f < -64) = -64;
p   = exp(f - max(f)); p = p/sum(p);
subplot(3,2,4), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f   = sum(F(:,:,3)); f = f - max(f(:)); f(f < -64) = -64;
p   = exp(f - max(f)); p = p/sum(p);
subplot(3,2,6), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (PEB)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square



% second level parameter estimates and Bayesian model comparison
%==========================================================================
spm_figure('GetWin','Figure 4'); clf

% and estimated second level parameters
%--------------------------------------------------------------------------
subplot(2,1,1), spm_plot_ci(BMA.Ep,BMA.Cp), hold on, bar([Ep;Ex],1/2), hold off
xlabel('parameters'), ylabel('expectation'), title('2nd level parameters','FontSize',16)
axis square

% random effects Bayesian model comparison
%--------------------------------------------------------------------------
[~,~,xp] = spm_dcm_bmc(RCM);

p   = full(spm_cat({REB.F})); p = exp(p - max(p)); p = p/sum(p);
subplot(2,2,3), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('posterior probability'), title('Random parameter effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p   = xp;
subplot(2,2,4), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('exceedance probability'), title('Random model effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square



% report log precision of random effects
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 5'); clf

subplot(2,2,1), spm_plot_ci(PEB.Eh,PEB.Ch),
title('Estimated log precision','FontSize',16)
axis square, a = axis;

subplot(2,2,2), bar(log((C - 0)/16))
title('True log precision','FontSize',16)
box off, axis square, axis(a)


% report (second level) model comparison using explicit PEB estimates (PF)
%--------------------------------------------------------------------------
spm_figure('Getwin','BMC - PEB'); clf

p = PF - max(PF(:));
p = exp(p);
p = p/sum(p(:));

subplot(3,2,1), imagesc(PF)
title('Free energy','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,3)
[m,i] = max(sum(p,2)); bar(sum(p,2)),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Commonalities','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,2), imagesc(p)
title('Posterior probabilities','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,4)
[m,i] = max(sum(p,1)); bar(sum(p,1))
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Differences','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square

subplot(3,2,6), imagesc(sparse(mw,mx,1,length(p),length(p)));
title('Correct solution','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

return


% Notes
%==========================================================================
hE    = linspace(-4,4,16);
hC    = 1/16;
clear Eh HF
for i = 1:length(hE)
    M.X     = X;
    M.hE    = hE(i);
    M.hC    = hC;
    PEB     = spm_dcm_peb(RCM(:,1),M);
    HF(i)   = PEB.F;
    Eh(:,i) = PEB.Eh;

end

subplot(2,2,1)
plot(hE,HF - max(HF))
subplot(2,2,2)
plot(hE,Eh)


% random field theory notes for number of the minimal
%--------------------------------------------------------------------------
W  = 1;                                                  % smoothness
u  = 1;                                                  % threshold
D  = 1:32;                                               % dimensionality

% expected number of maximum
%--------------------------------------------------------------------------
Em = exp(-u^2/2)*((2*pi).^(-(D + 1)/2)).*(W.^(-D)).*(u.^(D - 1));


