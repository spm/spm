% function DEMO_DCM_PEB
% Test routine to check group DCM for electrophysiology
%--------------------------------------------------------------------------
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
% First, an EEG DCM (using empirical grand mean data) is inverted to
% find plausible group mean parameters. Single subject data are
% then generated using typical within and between subject variance (here, 
% group differences in the modulation of intrinsic connectivity. We then
% illustrate a variety of Bayesian model averaging and reduction procedures
% to recover the underlying group effects.
%
% See also: spm_dcm_bmr, spm_dcm_peb and spm_dcm_peb_bma
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_DCM_PEB.m 6329 2015-02-05 19:25:52Z karl $


% change to directory with empirical data
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
mw  = 2;                              % true model (within)
mx  = 4;                              % true model (between)
Nm  = length(B);                      % number of models
Ns  = 16;                             % number of subjects
C   = 16;                             % within:between [co]variance ratio

% invert base model
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    DCM = rmfield(DCM,'M');
end
DCM.B   = B(mw);
DCM     = spm_dcm_erp(DCM);

% create subject-specifc DCM
%==========================================================================

% within subject effects:  condition specific effects 'B' (2 s.d.)
%--------------------------------------------------------------------------
sd          = sqrt(DCM.M.pC.B{1}(1,1));
sd          = sd/sqrt(C);
DCM.Ep.B{1} = B{mw}*2*sd;

% between subject effects: constant, group difference and covariance
%--------------------------------------------------------------------------
X           = [ones(Ns,1) kron([-1;1],ones(Ns/2,1)) randn(Ns,1)];
DCM.Ex      = spm_zeros(DCM.Ep);
DCM.Ex.B{1} = B{mx}*2*sd;

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
    Pp  = X(i,1)*spm_vec(DCM.Ep) + X(i,2)*spm_vec(DCM.Ex) + ep;
    Pp  = spm_unvec(Pp,DCM.Ep);
    Pg  = spm_vec(DCM.Eg) + spm_sqrtm(Cg)*randn(Ng,1);
    Pg  = spm_unvec(Pg,DCM.Eg);
    
    % generate data
    %----------------------------------------------------------------------
    G   = feval(DCM.M.G, Pg,DCM.M);
    x   = feval(DCM.M.IS,Pp,DCM.M,DCM.xU);
    for c = 1:length(x)
        y{c} = x{c}*G';
        e    = spm_conv(randn(size(y{c})),8,0);
        e    = e*mean(std(y{c}))/mean(std(e))/8;
        y{c} = y{c} + e;
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

% The following section contains the key analyses
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% second level model
%--------------------------------------------------------------------------
M     = struct('X',X);
M.pE  = GCM{1}.M.pE;
M.pC  = GCM{1}.M.pC;

% invert rreduced models (standard inversion)
%==========================================================================
% GCM   = spm_dcm_fit(GCM);

% Bayesian model reduction (avoiding local minima over models)
%==========================================================================
RCM   = spm_dcm_bmr(GCM);

% hierarchical (empirical Bayes) analysis using model reduction
%==========================================================================

% BMC – first level
%--------------------------------------------------------------------------
[BMC,PEB,PCM] = spm_dcm_bmc_peb(RCM,M,{'A','B'});

% BMA – first level
%--------------------------------------------------------------------------
bma   = spm_dcm_bma(GCM);
rma   = spm_dcm_bma(RCM);
pma   = spm_dcm_bma(PCM);

% BMC/BMA – second level
%--------------------------------------------------------------------------
BMA   = spm_dcm_peb_bmc(PEB,GCM(1,:));
PMA   = spm_dcm_bmr_all(PEB,{'A','B'});

% posterior predictive density and LOO cross validation
%==========================================================================
spm_dcm_loo(PCM,X,{'A','B'});

% full hierarchical inversion
%--------------------------------------------------------------------------
% [dcm,peb,Fpeb] = spm_dcm_peb_fit(GCM,M,{'A','B'});


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



% extract and plot results
%==========================================================================
clear Q
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
        F(i,j,1) = GCM{i,j}.F - GCM{i,1}.F;
        F(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
    end
    
end

% indices to plot parameters
%--------------------------------------------------------------------------
pC    = GCM{1,1}.M.pC;
c     = spm_vec(pC);
iA    = spm_fieldindices(Pp,'A');
iB    = spm_fieldindices(Pp,'B');
iA    = iA(find(c(iA)));
iB    = iB(find(c(iB)));


% plot simulation data
%==========================================================================
spm_figure('GetWin','Figure 1');clf

subplot(3,2,1)
plot(DCM.M.R*x{2}*G','k'), hold on
plot(x{2}*G',':k'),     hold off
xlabel('pst'), ylabel('response'), title('Signal (single subject)','FontSize',16)
axis square, spm_axis tight,  a = axis;

subplot(3,2,2)
plot(DCM.M.R*e,'k'), hold on
plot(e,':k'),     hold off
xlabel('pst'), ylabel('response'), title('Noise','FontSize',16)
axis square, spm_axis tight, axis(a)

p =  (X(:,2) > 0);
q = ~(X(:,2) > 0);

subplot(3,2,3)
plot(Y(:,q,1),'r'),  hold on
plot(Y(:,q,2),':r'), hold on
plot(Y(:,p,1),'b'),  hold on
plot(Y(:,p,2),':b'), hold off
xlabel('pst'), ylabel('response'), title('Group data','FontSize',16)
axis square, spm_axis tight

subplot(3,2,4)
plot(Y(:,q,1) - Y(:,q,2),'r'), hold on
plot(Y(:,p,1) - Y(:,p,2),'b'), hold off
xlabel('pst'), ylabel('differential response'), title('Difference waveforms','FontSize',16)
axis square, spm_axis tight

i = spm_fieldindices(DCM.Ep,'B{1}(1,1)');
j = spm_fieldindices(DCM.Ep,'B{1}(2,2)');

subplot(3,2,5)
plot(Q(i,q,1),Q(j,q,1),'.r','MarkerSize',32), hold on
plot(Q(i,p,1),Q(j,p,1),'.b','MarkerSize',32), hold off
xlabel('B{1}(1,1)'), ylabel('B{1}(2,2)'), title('Group effects','FontSize',16)
axis square

i = spm_fieldindices(DCM.Ep,'B{1}(3,3)');
j = spm_fieldindices(DCM.Ep,'B{1}(4,4)');

subplot(3,2,6)
plot(Q(i,q,1),Q(j,q,1),'or','MarkerSize',8), hold on
plot(Q(i,p,1),Q(j,p,1),'ob','MarkerSize',8), hold off
xlabel('B{1}(3,3)'), ylabel('B{1}(4,4)'), title('Group effects','FontSize',16)
axis square


% plot results: Bayesian model reduction vs. reduced models
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf

occ = 512;
f   = F(:,:,1); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (FFX)','FontSize',16)
axis square

f   = sum(f,1); f  = f - max(f) + occ; f(f < 0) = 0;
subplot(3,2,3), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

p   = softmax(f'); [m,i] = max(p); 
subplot(3,2,5), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

occ = 128;
f   = F(:,:,2); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,2), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

f   = sum(f,1); f  = f - max(f) + occ; f(f < 0) = 0;
subplot(3,2,4), bar(f), xlabel('model'), ylabel('Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square

p   = softmax(f'); [m,i] = max(p); 
subplot(3,2,6), bar(p)
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square


% a more detailed analysis of Bayesian model comparison
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf

pp = softmax(F(:,:,2)')';

f  = F(:,:,2); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

subplot(3,2,2), imagesc(pp)
xlabel('model'), ylabel('subject'), title('Model posterior (BMR)','FontSize',16)
axis square

[p i] = max(pp(:,1));
[p j] = min(pp(:,1));
stri  = sprintf('Subject %i',i);
strj  = sprintf('Subject %i',j);

subplot(3,2,3), bar(pp(i,:))
xlabel('model'), ylabel('probability'), title(stri,'FontSize',16)
axis square, spm_axis tight

subplot(3,2,4), bar(pp(j,:))
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


% first level parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4');clf, ALim = 1/2;

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


p   = spm_softmax(sum(F(:,:,1))');
subplot(3,2,2), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (FFX)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p   = spm_softmax(sum(F(:,:,2))');
subplot(3,2,4), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (BMR)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p   = BMC.Pw;
subplot(3,2,6), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (PEB)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square



% second level parameter estimates and Bayesian model comparison
%==========================================================================
spm_figure('GetWin','Figure 5');clf

% and estimated second level parameters
%--------------------------------------------------------------------------
Tp  = spm_vec(DCM.Ep);
Tx  = spm_vec(DCM.Ex);
Tp  = Tp(PEB.Pind);
Tx  = Tx(PEB.Pind);

i   = spm_fieldindices(DCM.Ep,'B');
i   = ismember(PEB.Pind,i);
j   = find([i; i]);
i   = find(i);

subplot(2,2,1), spm_plot_ci(PEB.Ep(j),PEB.Cp(j,j)), hold on, bar([Tp(i);Tx(i)],1/2), hold off
xlabel('parameters'), ylabel('expectation'), title('2nd level parameters','FontSize',16)
axis square


% random effects Bayesian model comparison
%==========================================================================
spm_figure('GetWin','Figure 6');clf
[~,~,xp] = spm_dcm_bmc(RCM);

p   = BMC.Pw;
subplot(2,2,1), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('posterior probability'), title('Random parameter effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square

p   = xp;
subplot(2,2,2), bar(p),[m,i] = max(p); 
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('exceedance probability'), title('Random model effects','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square


return


% classical inference
%==========================================================================

% classical inference of second level
%--------------------------------------------------------------------------
CVA = spm_cva(Q(iB,:,1)',X,[],[0 1 0]'); CP(1) = log(CVA.p);
CVA = spm_cva(Q(iB,:,2)',X,[],[0 1 0]'); CP(2) = log(CVA.p);
CVA = spm_cva(Q(iB,:,3)',X,[],[0 1 0]'); CP(3) = log(CVA.p);


% Bayesian model comparison
%--------------------------------------------------------------------------
field  = {'B'};
i  = 1;
PB = spm_dcm_peb(GCM(:,i),X(:,[1    ]),field);  BF(1) = PB.F;
PB = spm_dcm_peb(GCM(:,i),X(:,[1 2  ]),field);  BF(2) = PB.F;
PB = spm_dcm_peb(GCM(:,i),X(:,[1 3  ]),field);  BF(3) = PB.F;
PB = spm_dcm_peb(GCM(:,i),X(:,[1 2 3]),field);  BF(4) = PB.F;


% Notes
%==========================================================================
hE    = 0;
hC    = linspace(-4,4,16);
clear Eh HF
for i = 1:length(hC)
    M.X     = X(:,1:2);
    M.hE    = hE;
    M.hC    = exp(hC(i));
    PEB     = spm_dcm_peb(GCM(:,1),M,{'A','B'});
    HF(i)   = PEB.F;
    Eh(:,i) = PEB.Eh;

end

subplot(2,2,1)
plot(hC,HF - max(HF))
subplot(2,2,2)
plot(hC,Eh)

