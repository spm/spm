function E = ROBOT_DCM_LAPLACE
% Test routine to check group DCM for electrophysiology
% This routine illustrates the use of Bayesian model reduction when
% inverting hierarchical (dynamical) models; for example, multisubject DCM
% models. In this context, we have hierarchical models that are formally
% similar to parametric empirical Bayesian models - with the exception
% that the model of the first level can be nonlinear and dynamic. In brief,
% this routine shows how to finesse the brittleness of Bayesian model
% comparison at the single subject level by equipping the model with an
% extra (between subject) level. It illustrates the recovery of group
% effects on modulatory changes in effective connectivity (in the mismatch
% negativity paradigm) - based upon real data.
% 
% First, an EEG DCM (using empirical ggrand mean data) is inverted to
% find plausible group mean parameters. Single subject data are
% then generated using typical within and between subject variance (here, 
% group differences in the modulation of intrinsic connectivity. We then
% illustrate a variety of Bayesian model averaging and reduction procedures
% to recover the underlying group effects.

% $Id: ROBOT_DCM_LAPLACE.m 6299 2015-01-08 12:56:00Z guillaume $


% change to directory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
try
    cd('/home/spm/tests/MEEG')
catch
    cd('C:\Users\karl\Documents\SPM\DCM tests')
end
close all
delete(get(0,'Children'))
rng('default')
ALim = 3/4;

E = {};

% set up
%==========================================================================
load DCM_MMN                               % base DCM

DCM.options.spatial  = 'ECD';
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
Pm  = 2;                              % true model
Nm  = length(B);                      % number of models
Ns  = 16;                             % number of subjects
C   = 16;                             % within:between [co]variance ratio
occ = 64;                             % Ockham's window for display

% invert base model
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    DCM = rmfield(DCM,'M');
end
DCM.B = B(Pm);
DCM   = spm_dcm_erp(DCM);

% create subject-specifc DCM
%==========================================================================

% between subject effects: constant, group difference and covariance
%--------------------------------------------------------------------------
X     = [ones(Ns,1) kron([0;1],ones(Ns/2,1)) randn(Ns,1)];
gE    = spm_zeros(DCM.Ep);
gE.B{1}(1,1) = 1/4;
gE.B{1}(2,2) = 1/4;
gE.B{1}(3,3) = -1/4;
gE.B{1}(4,4) = -1/4;

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
    fprintf('\nCreating subject %i\n',i)
    
    
    % generate data
    %----------------------------------------------------------------------
    Ep  = X(i,2)*spm_vec(gE);
    Pp  = spm_vec(DCM.Ep) + Ep + spm_sqrtm(Cp)*randn(Np,1);
    Pp  = spm_unvec(Pp,DCM.Ep);
    Pg  = spm_vec(DCM.Eg) + spm_sqrtm(Cg)*randn(Ng,1);
    Pg  = spm_unvec(Pg,DCM.Eg);
    
    % generate data
    %----------------------------------------------------------------------
    G   = feval(DCM.M.G, Pg,DCM.M);
    x   = feval(DCM.M.IS,Pp,DCM.M,DCM.xU);
    for c = 1:length(x)
        e    = spm_pinv(DCM.M.R)*DCM.R{c}*spm_pinv(DCM.M.U);
        e    = spm_phase_shuffle(full(e));
        y{c} = x{c}*G' + 2*e;
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
        GCM{i,j}          = spm_dcm_erp(GCM{i,j});
    end
end

% Bayesian model reduction - over subjects
%--------------------------------------------------------------------------
for i = 1:Ns
    [q,P]    = spm_dcm_search_eeg(GCM(i,:));
    RCM(i,:) = P;
end

% get results
%==========================================================================

% BMA - over subjects
%--------------------------------------------------------------------------
bma   = spm_dcm_bma(GCM);
rma   = spm_dcm_bma(RCM);
for i = 1:Ns
    
    %  data - over subjects
    %----------------------------------------------------------------------
    Y(:,i,1) = GCM{i,1}.xY.y{1}*DCM.M.U(:,1);
    Y(:,i,2) = GCM{i,1}.xY.y{2}*DCM.M.U(:,1);
    
    % Parameter averages
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(bma.mEps{i});
    Q(:,i,2) = spm_vec(rma.mEps{i});
    Q(:,i,3) = spm_vec(GCM{i,1}.Tp);
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j)   = GCM{i,j}.F;
        R(i,j)   = RCM{i,j}.F;
    end
    
end

% select parameters
%--------------------------------------------------------------------------
Pp    = GCM{1,1}.M.pC;
c     = spm_vec(Pp);
iA    = spm_fieldindices(Pp,'A');
iB    = spm_fieldindices(Pp,'B');
iA    = iA(find(c(iA)));
iB    = iB(find(c(iB)));

% relative free energy
%--------------------------------------------------------------------------
F = F - F(:,1)*ones(1,Nm);


% plot simulation data
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf

subplot(3,2,1), plot(DCM.M.R*x{2}*G'), hold on
plot(x{2}*G',':'), hold off
xlabel('pst'), ylabel('response'), title('Signal','FontSize',16)
axis square, spm_axis tight,  a = axis;

subplot(3,2,2), plot(DCM.M.R*e), hold on
plot(e,':'), hold off
xlabel('pst'), ylabel('response'), title('Noise','FontSize',16)
axis square, spm_axis tight, axis(a)

subplot(3,2,3), plot(Y(:,:,1)), hold on
plot(Y(:,:,2),':'), hold off
xlabel('pst'), ylabel('response'), title('Data','FontSize',16)
axis square, spm_axis tight

subplot(3,2,4), plot(Y(:,:,1) - Y(:,:,2))
xlabel('pst'), ylabel('Condition effects'), title('Data','FontSize',16)
axis square, spm_axis tight

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


% plot results
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
HDM = spm_dcm_group(GCM(:,1),X(:,1),{'A','B'});

% (RFX) BMA - define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),length(HDM.Pind));
k     = spm_fieldindices(Pp,'B');
j     = find(ismember(HDM.Pind,k));
q     = find(spm_vec(B{1})');
for i = 1:length(B)
    m      = spm_vec(B{i});
    K(i,j) = m(q)';
end
    
% Bayesian model averaging
%--------------------------------------------------------------------------
BMA = spm_dcm_group_BMA(HDM,K);


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
HDM = spm_dcm_group(GCM(:,1),X(:,[1 1 1]),field);  HF(1) = HDM.F
HDM = spm_dcm_group(GCM(:,1),X(:,[1 1 2]),field);  HF(2) = HDM.F
HDM = spm_dcm_group(GCM(:,1),X(:,[1 1 3]),field);  HF(3) = HDM.F
HDM = spm_dcm_group(GCM(:,1),X(:,[1 2 3]),field);  HF(4) = HDM.F

spm_figure('GetWin','Figure 5');

p  = HF; p  = exp(p - max(p)); p = p/sum(p);
subplot(2,2,1), bar(p),[m i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
title('Posterior probability','FontSize',16)
xlabel('Model','FontSize',12)
ylabel('Probability','FontSize',12)
axis([0 (length(p) + 1) 0 1]), axis square
set(gca,'XTickLabel',{'1&1','1&2','1&3','1&2&3'})


% (RFX) BMA - define the model space in terms of a matrix
%--------------------------------------------------------------------------
K     = ones(length(B),length(HDM.Pind));
k     = spm_fieldindices(Pp,'B');
j     = find(ismember(HDM.Pind,k));
q     = find(spm_vec(B{1})');
for i = 1:length(B)
    m      = spm_vec(B{i});
    K(i,j) = m(q)';
end
    
% Bayesian model comparison of the second level
%--------------------------------------------------------------------------
BMA = spm_dcm_group_BMA(HDM,K);


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
    
    % (FFX) BPA - over subjects
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


% (FFX) BMA - over subjects
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
