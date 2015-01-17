% function DEMO_PEB_PEB
% Test routine to check group DCM for electrophysiology
%--------------------------------------------------------------------------
% This routine illustrates the use of Bayesian model reduction when
% inverting hierarchical (dynamical) models; for example, multisubject DCM
% models. In this context, we have hierarchical models that are formally
% similar to parametric empirical Bayesian models – with the exception
% that the model of the first level can be nonlinear and dynamic. In brief,
% this routine shows how to finesse the brittleness of Bayesian model
% comparison at the single subject level by equipping the model with an
% extra (between subject) level. It illustrates the recovery of group
% effects on modulatory changes in effective connectivity (in the mismatch
% negativity paradigm) – based upon real data.
% 
% First, an EEG DCM (using empirical ggrand mean data) is inverted to
% find plausible group mean parameters. Single subject data are
% then generated using typical within and between subject variance (here, 
% group differences in the modulation of intrinsic connectivity. We then
% illustrate a variety of Bayesian model averaging and reduction procedures
% to recover the underlying group effects.
%
% See also: spm_dcm_bmr, spm_dcm_peb and spm_dcm_peb_bma
%__________________________________________________________________________
% Copyright (C) 2010-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_BMR_PEB.m 6305 2015-01-17 12:40:51Z karl $


% change to ddirectory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
close all, clear all
clc, rng('default'), drawnow


% set up
%==========================================================================

% model space - within subject effects
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


% (RFX) BMA – define the model space in terms of a matrix
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
        
        P{2}.C    = diag(K(j,:))*pC*diag(K(j,:));
        [qP,D,F]  = spm_PEB(y{i,1},P,1);
        
        GCM{i,j}.M.pE = DCM.M.pE;
        GCM{i,j}.M.pC = P{2}.C;
        GCM{i,j}.Ep   = qP{2}.E;
        GCM{i,j}.Cp   = qP{2}.C;
        GCM{i,j}.B    = B(j);
        GCM{i,j}.F    = F;
        GCM{i,j}.Pp   = Pp;
        
    end
end

% PEB (GLM)
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
        pCi        = diag(K(i,:))*pC*diag(K(i,:));
        pCj        = diag(K(j,:))*pC*diag(K(j,:));
        P{3}.C     = blkdiag(pCi,pCj);
        [qP,D,F]   = spm_PEB(spm_cat(y),P,1);
        
        PB(i,j).F  = F;
        PB(i,j).Ep = qP{3}.E;
        PB(i,j).Cp = qP{3}.C;
        
        PF(i,j)    = F;
        
    end
end


% show second level model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf

p = PF - max(PF(:));
p = exp(p);
p = p/sum(p(:));

spm_figure('Getwin','BMC - PEB'); clf

subplot(3,2,1), imagesc(PF)
title('Free energy','FontSize',16)
xlabel('Model (differences)','FontSize',12)
ylabel('Model (commonalities)','FontSize',12)
axis square

subplot(3,2,3)
[m i] = max(sum(p,1)); bar(sum(p,1)),
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
[m i] = max(sum(p,2)); bar(sum(p,2)),
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Differences','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (Nm + 1) 0 1]), axis square


% Bayesian model reduction – for each subject
%==========================================================================
[RCM,BMR] = spm_dcm_bmr(GCM);




% hierarchical (RFX) analysis
%==========================================================================
RCM{1,1}.M.eE = 4;
[REB,PCM] = spm_dcm_peb(RCM);

% BMA – (first level)
%--------------------------------------------------------------------------
bma   = spm_dcm_bma(GCM);
rma   = spm_dcm_bma(RCM);
pma   = spm_dcm_bma(PCM);

% BMA – (second level)
%--------------------------------------------------------------------------
PEB   = spm_dcm_peb(RCM(:,1),X);
BMA   = spm_dcm_peb_bmc(PEB,RCM(1,:));


% show results
%==========================================================================
clear Q
for i = 1:Ns
    
    % Parameter averages
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Pp);
    Q(:,i,2) = spm_vec(bma.SUB(i).Ep);
    Q(:,i,3) = spm_vec(rma.SUB(i).Ep);
    Q(:,i,4) = spm_vec(pma.SUB(i).Ep);
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,1) = GCM{i,j}.F;
        F(i,j,2) = RCM{i,j}.F;
        F(i,j,3) = REB(j).F;
    end
    
end

% select parameters
%--------------------------------------------------------------------------
iA    = spm_fieldindices(DCM.M.pE,'A');
iB    = spm_fieldindices(DCM.M.pE,'B');


% plot results: Bayesian model reduction vs. reduced models
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf

f  = F(:,:,1); f = f - max(f(:)); f(f < -64) = -64;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (FFX)','FontSize',16)
axis square

f  = sum(f,1); f  = f - min(f);
subplot(3,2,3), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

p  = exp(f - max(f)); p = p/sum(p); [m i] = max(p); 
subplot(3,2,5), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f  = F(:,:,2); f = f - max(f(:)); f(f < -64) = -64;
subplot(3,2,2), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

f  = sum(f,1); f  = f - min(f);
subplot(3,2,4), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square

p  = exp(f - max(f)); p = p/sum(p); [m i] = max(p); 
subplot(3,2,6), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square



% parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');clf

ALim = 1/2;

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,2)));
str = sprintf('BMA: correlation = %-0.2f',r);
subplot(3,2,1), plot(Q(iA,:,1),Q(iA,:,2),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,2),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,3)));
str = sprintf('BMR: correlation = %-0.2f',r);
subplot(3,2,3), plot(Q(iA,:,1),Q(iA,:,3),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,3),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

r   = corr(spm_vec(Q([iA; iB],:,1)),spm_vec(Q([iA; iB],:,4)));
str = sprintf('PEB: correlation = %-0.2f',r);
subplot(3,2,5), plot(Q(iA,:,1),Q(iA,:,4),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Q(iB,:,4),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title(str,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

f  = sum(F(:,:,1)); f = f - max(f(:)); f(f < -64) = -64;
p  = exp(f - max(f)); p = p/sum(p);
subplot(3,2,2), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f  = sum(F(:,:,2)); f = f - max(f(:)); f(f < -64) = -64;
p  = exp(f - max(f)); p = p/sum(p);
subplot(3,2,4), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

f  = sum(F(:,:,3)); f = f - max(f(:)); f(f < -64) = -64;
p  = exp(f - max(f)); p = p/sum(p);
subplot(3,2,6), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (PEB)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square




return



Ep  = spm_cat({BMA.SUB.Ep});
i   = ismember(BMA.Pind,iA);
j   = ismember(BMA.Pind,iB);

subplot(3,2,5), plot(Q(iA,:,1),Q(iA,:,1),'.c','MarkerSize',16), hold on
plot(Q(iB,:,1),Ep(j,:),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title('Parameters (RFX)','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

p   = sum(BMA.P,1);
subplot(3,2,6), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Posterior (RFX)','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('probability','FontSize',12)
axis([0 (length(p) + 1) 0 1]), axis square



% Notes
%==========================================================================
eE    = linspace(-4,4,16);
eC    = 1/2;
for i = 1:length(hE)
    RCM{1,1}.M.eE = eE(i);
    RCM{1,1}.M.eC = eC;
    PEB   = spm_dcm_peb(RCM(:,1));
    HF(i) = PEB.F;
    Eh(:,i) = PEB.Eh;

end

subplot(2,2,1)
plot(eE,HF)
subplot(2,2,2)
plot(eE,Eh)

W  = 1;
u  = 1;
D  = 1:32;

Em = exp(-u^2/2)*((2*pi).^(-(D + 1)/2)).*(W.^(-D)).*(u.^(D - 1));


