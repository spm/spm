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
% $Id: DEMO_DCM_PEB.m 6305 2015-01-17 12:40:51Z karl $


% change to ddirectory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
try
    cd('C:\home\spm\DCM\DCM tests')
catch
    cd('C:\Users\karl\Documents\SPM\DCM tests')
end
close all, clear all
clc
rng('default')
ALim = 3/4;

% set up
%==========================================================================
load DCM_MMN                               % base DCM

DCM.options.spatial  = 'IMG';
DCM.options.analysis = 'ERP';
DCM.options.model    = 'ERP';
DCM.options.Nmax     = 32;
DCM.options.DATA     = 1;
DCM.name             = 'DCM_GROUP';

% model space - within subject effects
%--------------------------------------------------------------------------
k     = spm_perm_mtx(3);
for i = 1:8;
    B{i}     = sparse(5,5);
    if k(i,1)
        B{i} = B{i} + sparse([1 2 3 4],[1 2 3 4],1,5,5);
    end
    if k(i,2)
        B{i} = B{i} + sparse([1 2],[3 4],1,5,5);
    end
    if k(i,3)
        B{i} = B{i} + sparse([3 4],[1 2],1,5,5);
    end
    B{i}     = full(B{i});
end


% model space
%--------------------------------------------------------------------------
mw  = 3;                              % true model (within)
mx  = 4;                              % true model (between)
Nm  = length(B);                      % number of models
Ns  = 16;                             % number of subjects
C   = 32;                             % within:between [co]variance ratio
occ = 64;                             % Ockham's window for display

% invert base model
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    DCM = rmfield(DCM,'M');
end
DCM.B = B(mw);
DCM   = spm_dcm_erp(DCM);

% create subject-specifc DCM
%==========================================================================

% within subject effects:  condition specific effects 'B' (2 s.d.)
%--------------------------------------------------------------------------
sd          = sqrt(DCM.M.pC.B{1}(1,1));
sd          = sd/sqrt(C);
DCM.Ep.B{1} = B{mw}*2*sd;

% between subject effects: constant, group difference and covariance
%--------------------------------------------------------------------------
X     = [ones(Ns,1) kron([0;1],ones(Ns/2,1)) randn(Ns,1)];
gE    = spm_zeros(DCM.Ep);
gE.B{1}(1,1) =  4*sd;
gE.B{1}(2,2) =  4*sd;
gE.B{1}(3,3) =  2*sd;
gE.B{1}(4,4) =  2*sd;

% create subject-specifc DCM
%--------------------------------------------------------------------------
DCM.options.DATA = 0;
DCM   = spm_dcm_erp_dipfit(DCM,1);

Np    = spm_length(DCM.M.pE);
Ng    = spm_length(DCM.M.gE);
Cp    = diag(spm_vec(DCM.M.pC))/C;
Cg    = diag(spm_vec(DCM.M.gC))/C;
for i = 1:Ns
    
    % report
    %----------------------------------------------------------------------
    fprintf('Creating subject %i\n',i)
    
    
    % generate data
    %----------------------------------------------------------------------
    ep  = spm_sqrtm(Cp)*randn(Np,1);
    Pp  = X(i,1)*spm_vec(DCM.Ep) + X(i,2)*spm_vec(gE) + ep;
    Pp  = spm_unvec(Pp,DCM.Ep);
    Pg  = spm_vec(DCM.Eg) + spm_sqrtm(Cg)*randn(Ng,1);
    Pg  = spm_unvec(Pg,DCM.Eg);
    
    % generate data
    %----------------------------------------------------------------------
    G   = feval(DCM.M.G, Pg,DCM.M);
    x   = feval(DCM.M.IS,Pp,DCM.M,DCM.xU);
    for c = 1:length(x)
        e    = spm_pinv(DCM.M.R)*DCM.R{c}*spm_pinv(DCM.M.U);
        e    = spm_phase_shuffle(full(e))*2;
        y{c} = x{c}*G' + e;
        y{c} = DCM.M.R*y{c};
    end
    
    % invert models
    %----------------------------------------------------------------------
    for j = 1:Nm
        GCM{i,j}          = rmfield(DCM,'M');
        GCM{i,j}.M.dipfit = DCM.M.dipfit;
        GCM{i,j}.B        = B(j);
        GCM{i,j}.xY.y     = y;
        GCM{i,j}.Tp       = Pp;
        GCM{i,j}.Tg       = Pg;
    end
end

% invert full models (first column)
%==========================================================================
GCM       = spm_dcm_fit(GCM);

% Bayesian model reduction – for each subject
%==========================================================================
[RCM,BMR] = spm_dcm_bmr(GCM);

% hierarchical (RFX) analysis
%==========================================================================
[PEB,PCM] = spm_dcm_peb(RCM);

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
for i = 1:Ns
    
    %  data – over subjects
    %----------------------------------------------------------------------
    Y(:,i,1) = GCM{i,1}.xY.y{1}*DCM.M.U(:,1);
    Y(:,i,2) = GCM{i,1}.xY.y{2}*DCM.M.U(:,1);
    
    % Parameter averages
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Tp);
    Q(:,i,2) = spm_vec(bma.SUB(i).Ep);
    Q(:,i,3) = spm_vec(rma.SUB(i).Ep);
    Q(:,i,4) = spm_vec(pma.SUB(i).Ep);
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,1) = GCM{i,j}.F;
        F(i,j,2) = RCM{i,j}.F;
        F(i,j,3) = PCM{i,j}.F;
    end
    
end

% select parameters
%--------------------------------------------------------------------------
pC    = GCM{1,1}.M.pC;
c     = spm_vec(pC);
iA    = spm_fieldindices(Pp,'A');
iB    = spm_fieldindices(Pp,'B');
iA    = iA(find(c(iA)));
iB    = iB(find(c(iB)));


% plot simulation data
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf

subplot(3,2,1)
plot(DCM.M.R*x{2}*G'), hold on
plot(x{2}*G',':'),     hold off
xlabel('pst'), ylabel('response'), title('Signal (single subject)','FontSize',16)
axis square, spm_axis tight,  a = axis;

subplot(3,2,2)
plot(DCM.M.R*e), hold on
plot(e,':'),     hold off
xlabel('pst'), ylabel('response'), title('Noise','FontSize',16)
axis square, spm_axis tight, axis(a)

subplot(3,2,3)
plot(Y(:,:,1)),     hold on
plot(Y(:,:,2),':'), hold off
xlabel('pst'), ylabel('response'), title('Group data','FontSize',16)
axis square, spm_axis tight

subplot(3,2,4)
plot(Y(:,:,1) - Y(:,:,2))
xlabel('pst'), ylabel('differential response'), title('Difference waveforms','FontSize',16)
axis square, spm_axis tight

i = spm_fieldindices(DCM.Ep,'B{1}(1,1)');
j = spm_fieldindices(DCM.Ep,'B{1}(2,2)');

subplot(3,2,5), 
plot(Q(i, ~X(:,2),1),Q(j, ~X(:,2),1),'.r','MarkerSize',32), hold on
plot(Q(i,~~X(:,2),1),Q(j,~~X(:,2),1),'.b','MarkerSize',32), hold off
xlabel('B{1}(1,1)'), ylabel('B{1}(2,2)'), title('Group effects','FontSize',16)
axis square

i = spm_fieldindices(DCM.Ep,'B{1}(3,3)');
j = spm_fieldindices(DCM.Ep,'B{1}(4,4)');

subplot(3,2,6), 
plot(Q(i, ~X(:,2),1),Q(j, ~X(:,2),1),'or','MarkerSize',8), hold on
plot(Q(i,~~X(:,2),1),Q(j,~~X(:,2),1),'ob','MarkerSize',8), hold off
xlabel('B{1}(3,3)'), ylabel('B{1}(4,4)'), title('Group effects','FontSize',16)
axis square


% plot results: Bayesian model reduction vs. reduced models
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');clf

subplot(3,2,1), image((F - max(F(:)))/4 + 64)
xlabel('model'), ylabel('subject'), title('Free energy','FontSize',16)
axis square

subplot(3,2,2), image((R - max(R(:)))*2 + 64)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

subplot(3,2,3)
GF  = sum(F); GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

subplot(3,2,4)
GF  = sum(R); GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square

subplot(3,2,5)
GF  = sum(F); GF  = exp(GF - max(GF)); GF = GF/sum(GF);
bar(GF), xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
spm_axis tight, axis square

subplot(3,2,6)
GF  = sum(R); GF  = exp(GF - max(GF)); GF = GF/sum(GF);
bar(GF), xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
spm_axis tight, axis square


% parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');clf

subplot(3,2,1), plot(Q(iA,:,3),Q(iA,:,1),'.c','MarkerSize',16), hold on
plot(Q(iB,:,3),Q(iB,:,1),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title('Parameters (BMA)','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

subplot(3,2,3), plot(Q(iA,:,3),Q(iA,:,2),'.c','MarkerSize',16), hold on
plot(Q(iB,:,3),Q(iB,:,2),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title('Parameters (BMR)','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

GF  = sum(F); GF  = exp(GF - max(GF)); p = GF/sum(GF);
subplot(3,2,2), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

GF  = sum(R); GF  = exp(GF - max(GF)); p = GF/sum(GF);
subplot(3,2,4), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

% a more detailed analysis of Bayesian model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf

for i = 1:Ns
    p       = exp(R(i,:));
    p       = p/sum(p);
    RR(i,:) = p;
end

subplot(3,2,1), image((R - max(R(:)))*2 + 64)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

subplot(3,2,2), imagesc(RR)
xlabel('model'), ylabel('subject'), title('Model posterior (BMR)','FontSize',16)
axis square

[p i] = max(RR(:,1));
[p j] = min(RR(:,1));
stri  = sprintf('Subject %i',i);
strj  = sprintf('Subject %i',j);

subplot(3,2,3), bar(RR(i,:))
xlabel('model'), ylabel('probability'), title(stri,'FontSize',16)
axis square, spm_axis tight

subplot(3,2,4), bar(RR(j,:))
xlabel('model'), ylabel('probability'), title(strj,'FontSize',16)
axis square, spm_axis tight

k   = spm_fieldindices(DCM.Ep,'B');
pE  = RCM{i,1}.Tp.B{1}; pE = spm_vec(pE);
qE  = RCM{i,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = RCM{i,1}.Cp(k,k); qC = diag(qC);
pE  = pE(find(qC));
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,5), spm_plot_ci(qE,qC), hold on, bar(pE,1/2), hold off
xlabel('parameter (B)'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, a = axis;

pE  = RCM{j,1}.Tp.B{1}; pE = spm_vec(pE);
qE  = RCM{j,1}.Ep.B{1}; qE = spm_vec(qE);
qC  = RCM{j,1}.Cp(k,k); qC = diag(qC);
pE  = pE(find(qC));
qE  = qE(find(qC));
qC  = qC(find(qC));

subplot(3,2,6), spm_plot_ci(qE,qC), hold on, bar(pE,1/2), hold off
xlabel('parameter (B)'), ylabel('expectation'), title('Parameters','FontSize',16)
axis square, axis(a);


% hierarchical (RFX) analysis
%==========================================================================
PEB = spm_dcm_peb(GCM(:,1),X(:,1),{'A','B'});

% (RFX) BMA – define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),length(PEB.Pind));
k     = spm_fieldindices(Pp,'B');
j     = find(ismember(PEB.Pind,k));
q     = find(spm_vec(B{1})');
for i = 1:length(B)
    m      = spm_vec(B{i});
    K(i,j) = m(q)';
end
    
% Bayesian model averaging
%--------------------------------------------------------------------------
BMA = spm_dcm_group_BMA(PEB,K);


spm_figure('GetWin','Figure 3');

Ep  = spm_cat({BMA.SUB.Ep});
i   = ismember(BMA.Pind,iA);
j   = ismember(BMA.Pind,iB);

subplot(3,2,5), plot(Q(iA,:,3),Ep(i,:),'.c','MarkerSize',16), hold on
plot(Q(iB,:,3),Ep(j,:),'.b','MarkerSize',16), hold off
xlabel('true parameter'), ylabel('Model average'), title('Parameters (RFX)','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

p   = BMA.P;
subplot(3,2,6), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Posterior (RFX)','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('probability','FontSize',12)
axis([0 (length(p) + 1) 0 1]), axis square


% inference
%==========================================================================

% classical inference of second level
%--------------------------------------------------------------------------
CVA = spm_cva(Q(iB,:,3)',X,[],[0 1 0]'); CP(1) = log(CVA.p);
CVA = spm_cva(Q(iB,:,1)',X,[],[0 1 0]'); CP(2) = log(CVA.p);
CVA = spm_cva(Q(iB,:,2)',X,[],[0 1 0]'); CP(3) = log(CVA.p);
CVA = spm_cva(Ep(j,:)',  X,[],[0 1 0]'); CP(4) = log(CVA.p);


% Bayesian model comparison
%--------------------------------------------------------------------------
field  = {'A','B'};
PEB = spm_dcm_group(GCM(:,1),X(:,[1 1 1]),field);  HF(1) = PEB.F
PEB = spm_dcm_group(GCM(:,1),X(:,[1 1 2]),field);  HF(2) = PEB.F
PEB = spm_dcm_group(GCM(:,1),X(:,[1 1 3]),field);  HF(3) = PEB.F
PEB = spm_dcm_group(GCM(:,1),X(:,[1 2 3]),field);  HF(4) = PEB.F

spm_figure('GetWin','Figure 5');

p  = HF; p  = exp(p - max(p)); p = p/sum(p);
subplot(2,2,1), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Posterior probability','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (length(p) + 1) 0 1]), axis square
set(gca,'XTickLabel',{'1&1','1&2','1&3','1&2&3'})


% (RFX) BMA – define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),length(PEB.Pind));
k     = spm_fieldindices(Pp,'B');
j     = find(ismember(PEB.Pind,k));
q     = find(spm_vec(B{1})');
for i = 1:length(B)
    m      = spm_vec(B{i});
    K(i,j) = m(q)';
end
    
% Bayesian model comparison of the second level
%--------------------------------------------------------------------------
BMA = spm_dcm_group_BMA(PEB,K);


i = spm_fieldindices(DCM.Ep,'B{1}(1,1)');
j = spm_fieldindices(DCM.Ep,'B{1}(2,2)');

subplot(3,2,5), plot(Q(i,~X(:,2),3),Q(j,~X(:,2),3),'.r','MarkerSize',32), hold on
plot(Q(i,~~X(:,2),3),Q(j,~~X(:,2),3),'.b','MarkerSize',32), hold off
xlabel('B{1}(1,1)'), ylabel('B{1}(2,2)'), title('Group effects','FontSize',16)
axis square

i = spm_fieldindices(DCM.Ep,'B{1}(3,3)');
j = spm_fieldindices(DCM.Ep,'B{1}(4,4)');

subplot(3,2,6), plot(Q(i,~X(:,2),3),Q(j,~X(:,2),3),'or','MarkerSize',8), hold on
plot(Q(i,~~X(:,2),3),Q(j,~~X(:,2),3),'ob','MarkerSize',8), hold off
xlabel('B{1}(3,3)'), ylabel('B{1}(4,4)'), title('Group effects','FontSize',16)
axis square












% reinvert (full) model with initialization; recursively
%==========================================================================
for k = 1:8
    
    % (FFX) BPA – over subjects
    %----------------------------------------------------------------------
    BPA   = spm_dcm_average(RCM(:,1),{},1);
    
    for i = 1:Ns
        for j = 1:Nm
            try
                RCM{i,j}  = rmfield(RCM{i,j},'M');
            end
        end
        
        % invert the full (first) model
        %------------------------------------------------------------------
        RCM{i,1}.M.P      = BPA.Ep;
        RCM{i,1}.M.dipfit = DCM.M.dipfit;
        RCM{i,1}          = spm_dcm_erp(RCM{i,1});
        
        % correlations and free energy
        %------------------------------------------------------------------
        qp     = spm_vec(RCM{i,1}.Ep);
        pp     = spm_vec(RCM{i,1}.Tp);
        rho    = corr([qp(iq) pp(iq)]);
        cor(i) = rho(1,2);
        f(i)   = RCM{i,1}.F;
    end
    
    COR(k) = mean(cor)
    IF(k)  = mean(f)
end

% Bayesian model reduction
%--------------------------------------------------------------------------
for i = 1:Ns
    [q,P]    = spm_dcm_search_eeg(RCM(i,:));
    RCM(i,:) = P;
end


% (FFX) BMA – over subjects
%----------------------------------------------------------------------
rma   = spm_dcm_bma(RCM);
for i = 1:Ns
    
    % Parameter averages
    %------------------------------------------------------------------
    QP(:,i) = spm_vec(rma.mEps{i});
    
    % and free energy
    %----------------------------------------------------------------------
    for j = 1:Nm
        RF(i,j)   = RCM{i,j}.F;
    end
    
end


% select parameters
%--------------------------------------------------------------------------
qp     = QP(iq,:);


% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');clf

subplot(3,2,1), plot(Q(:,:,3),qp,'o'), hold on
xlabel('true parameter'), ylabel('Model average'), title('Parameters','FontSize',16)
axis square

subplot(3,2,2), plot(Q(:,:,3),Q(:,:,2),'o'), hold on
xlabel('true parameter'), ylabel('BMA (BMR)'), title('Parameters','FontSize',16)
axis square

subplot(3,2,3), imagesc(RF)
xlabel('model'), ylabel('subject'), title('Free energy','FontSize',16)
axis square

subplot(3,2,4), imagesc(R)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

subplot(3,2,5)
GF  = sum(RF);
GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

subplot(3,2,6)
GF  = sum(R);
GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('FFX Free energy'), title('BMR Free energy','FontSize',16)
spm_axis tight, axis square
