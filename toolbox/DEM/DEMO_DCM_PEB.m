function DEMO_DCM_PEB
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
% $Id: DEMO_DCM_PEB.m 6759 2016-03-27 19:45:17Z karl $


% change to directory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
rng('default')
corr = @(x,y) subsref(corrcoef(x,y),substruct('()',{1,2})); % Stats tbx


% Set up
%==========================================================================
try
    load DEMO_DCM_PEB_MAT
catch
    
    % get base DCM and generate data
    %----------------------------------------------------------------------
    load DCM_MMN
    try
        D  = spm_eeg_load(DCM.xY.Dfile);
        if ~spm_existfile(D.inv{end}.forward.vol)
            D.inv{end}.forward.vol = spm_file(D.inv{end}.forward.vol,...
                'path',fullfile(spm('Dir'),'canonical'));
            save(D);
        end
    end
    
    DCM.options.spatial  = 'ECD';
    DCM.options.analysis = 'ERP';
    DCM.options.model    = 'ERP';
    DCM.options.Nmax     = 32;
    DCM.options.DATA     = 1;
    DCM.name             = 'DCM_GROUP';
    
    % model space - within subject effects
    %----------------------------------------------------------------------
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
    DCM.Ep.B{1} = [
        0.1  0    0   0   0;
        0    0.1  0   0   0;
        0.3  0    0.3 0   0;
        0    0.3  0   0.3 0;
        0    0    0   0   0];
    
    % between subject effects: constant, group difference and covariance
    %--------------------------------------------------------------------------
    X           = [ones(Ns,1) kron([-1;1],ones(Ns/2,1)) randn(Ns,1)];
    DCM.Ex      = spm_zeros(DCM.Ep);
    DCM.Ex.B{1} = -B{mx}*2*sd;
    
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
        xY{i}    = y;
        
        % specify models
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
    
    % invert rreduced models (standard inversion)
    %==========================================================================
    GCM   = spm_dcm_fit(GCM);
    
    
    % remove redundant fields
    %----------------------------------------------------------------------
    KEEP  = {'M','Ep','Cp','F','pE','pC','Tp','Ex','U','R'};
    str   = fieldnames(DCM);
    for i = 1:numel(str)
        if ~ismember(str{i},KEEP)
            DCM = rmfield(DCM,str{i});
        end
    end
    str   = fieldnames(DCM.M);
    for i = 1:numel(str)
        if ~ismember(str{i},KEEP)
            DCM.M = rmfield(DCM.M,str{i});
        end
    end
    
    for j = 1:numel(GCM)
        str   = fieldnames(GCM{j});
        for i = 1:numel(str)
            if ~ismember(str{i},KEEP)
                GCM{j} = rmfield(GCM{j},str{i});
            end
        end
        str   = fieldnames(GCM{j}.M);
        for i = 1:numel(str)
            if ~ismember(str{i},KEEP)
                GCM{j}.M = rmfield(GCM{j}.M,str{i});
            end
        end
        
    end
    
    % save variables for subsequent analysis
    %----------------------------------------------------------------------
    save DEMO_DCM_PEB_MAT DCM GCM X Ns G mw xY x e
    
end


% The following section contains the key analyses
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% Bayesian model reduction (avoiding local minima over models)
%==========================================================================
RCM   = spm_dcm_bmr(GCM);

% hierarchical (empirical Bayes) model reduction
%==========================================================================
[peb,PCM] = spm_dcm_peb(RCM,[],{'A','B'});


% BMA - first level
%--------------------------------------------------------------------------
bma  = spm_dcm_bma(GCM);
rma  = spm_dcm_bma(RCM);
pma  = spm_dcm_bma(PCM);

% BMC/BMA - second level
%==========================================================================

% second level model
%--------------------------------------------------------------------------
M    = struct('X',X);

% randomisation analysis
%--------------------------------------------------------------------------
spm_dcm_peb_rnd(RCM(:,mw),M,{'A','B'});

% BMC - search over first and second level effects
%--------------------------------------------------------------------------
[BMC,PEB] = spm_dcm_bmc_peb(RCM,M,{'A','B'});

% BMA - exhaustive search over second level parameters
%--------------------------------------------------------------------------
PEB.gamma = 1/128;
BMA       = spm_dcm_peb_bmc(PEB);

% overlay true values
%--------------------------------------------------------------------------
Tp        = spm_vec(DCM.Ep);           % true second level paramters
Tx        = spm_vec(DCM.Ex);           % true second level paramters

subplot(3,3,1),hold on, bar(Tp(BMA.Pind),1/2), hold off
subplot(3,3,4),hold on, bar(Tp(BMA.Pind),1/2), hold off
subplot(3,3,2),hold on, bar(Tx(BMA.Pind),1/2), hold off
subplot(3,3,5),hold on, bar(Tx(BMA.Pind),1/2), hold off
drawnow

if 0
    % alternative search over restricted model space at second level
    %----------------------------------------------------------------------
    spm_dcm_peb_bmc(PEB,RCM(1,:));
    
end

% posterior predictive density and LOO cross validation
%==========================================================================
spm_dcm_loo(RCM(:,1),X,{'B'});

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



% extract and plot results
%==========================================================================
[Ns Nm] = size(GCM);
clear P Q R
for i = 1:Ns
    
    % data - over subjects
    %----------------------------------------------------------------------
    Y(:,i,1) = xY{i}{1}*DCM.M.U(:,1);
    Y(:,i,2) = xY{i}{2}*DCM.M.U(:,1);
    
    % Parameter estimates under true model
    %----------------------------------------------------------------------
    P(:,i,1) = spm_vec(GCM{i,mw}.Tp);
    P(:,i,2) = spm_vec(GCM{i,mw}.Ep);
    P(:,i,3) = spm_vec(RCM{i,mw}.Ep);
    P(:,i,4) = spm_vec(PCM{i,mw}.Ep);
    
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
        F(i,j,3) = PCM{i,j}.F - PCM{i,1}.F;
    end
    
end


% classical inference
%==========================================================================

% indices to plot parameters
%--------------------------------------------------------------------------
pC    = GCM{1,1}.M.pC;
iA    = spm_find_pC(pC,pC,'A');
iB    = spm_find_pC(pC,pC,'B');


% classical inference of second level
%--------------------------------------------------------------------------
i   = [iA;iB];
CVA = spm_cva(Q(i,:,1)',X,[],[0 1 0]'); CP(1) = log(CVA.p); CR(1) = corr(Tx(i),CVA.V);
CVA = spm_cva(Q(i,:,2)',X,[],[0 1 0]'); CP(2) = log(CVA.p); CR(2) = corr(Tx(i),CVA.V);
CVA = spm_cva(Q(i,:,3)',X,[],[0 1 0]'); CP(3) = log(CVA.p); CR(3) = corr(Tx(i),CVA.V);
CVA = spm_cva(Q(i,:,4)',X,[],[0 1 0]'); CP(4) = log(CVA.p); CR(4) = corr(Tx(i),CVA.V);

spm_figure('GetWin','CVA');clf
subplot(2,2,1), bar(spm_en([Tx(i) CVA.V]))
title('Canonical vector','FontSize',16)
xlabel('parameter'), ylabel('weight'), axis square

subplot(2,2,2), bar(spm_en([X(:,2) CVA.v]))
title('Canonical variate','FontSize',16)
xlabel('parameter'), ylabel('weight'), axis square
legend({'RFX','True'})

subplot(2,1,2), bar(abs(CR(2:end)))
title('Correlations with true values','FontSize',16)
xlabel('inversion scheme'), ylabel('correlation'), axis square
set(gca,'XTickLabel',{'FFX','BMR','PEB'})

% correlations with true values
%==========================================================================
for i = 1:3
    R(i)    = corr(spm_vec(P([iA; iB],:,1)),spm_vec(P([iA; iB],:,i + 1)));
    Rstr{i} = num2str(R(i));
end

% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Correlations');clf
subplot(2,2,1)
plot(spm_vec(P([iA; iB],:,2)),spm_vec(P([iA; iB],:,3)),'.','MarkerSize',8)
title('Parameter estimates','FontSize',16)
xlabel('mean (standard inverion)'), ylabel('BMR'), axis square

subplot(2,2,2), bar(R)
title('Correlations with true values','FontSize',16)
xlabel('inversion scheme'), ylabel('correlation'), axis square
text((1:3) - 1/8,R/2,Rstr(:),'Color','w','Fontsize',10)
set(gca,'XTickLabel',{'FFX','BMR','PEB'})



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
plot(Q(i,q,1),Q(j,q,1),'.r','MarkerSize',24), hold on
plot(Q(i,p,1),Q(j,p,1),'.b','MarkerSize',24), hold off
xlabel('B{1}(1,1)'), ylabel('B{1}(2,2)'), title('Group effects','FontSize',16)
axis square

i = spm_fieldindices(DCM.Ep,'B{1}(3,3)');
j = spm_fieldindices(DCM.Ep,'B{1}(4,4)');

subplot(3,2,6)
plot(Q(i,q,1),Q(j,q,1),'.r','MarkerSize',24), hold on
plot(Q(i,p,1),Q(j,p,1),'.b','MarkerSize',24), hold off
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

f   = F(:,:,2); f = f - max(f(:)) + occ; f(f < 0) = 0;
subplot(3,2,1), imagesc(f)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

pp  = softmax(f')';
subplot(3,2,2), imagesc(pp)
xlabel('model'), ylabel('subject'), title('Model posterior (BMR)','FontSize',16)
axis square

[p,i] = max(pp(:,1));
[p,j] = min(pp(:,1));
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

p   = spm_softmax(sum(F(:,:,3))');
subplot(3,2,6), bar(p),[m,i] = max(p);
text(i - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
xlabel('model'), ylabel('probability'), title('Posterior (PEB)','FontSize',16)
axis([0 (length(p) + 1) 0 1]), axis square


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



% Notes
%==========================================================================
hE    = linspace(-4,4,16);
hC    = 1;
clear Eh HF
for i = 1:length(hE)
    M.X     = X(:,1:2);
    M.hE    = hE(i);
    M.hC    = 1;
    PEB     = spm_dcm_peb(GCM(:,1),M,{'A','B'});
    HF(i)   = PEB.F;
    Eh(:,i) = PEB.Eh;
    
end

subplot(2,2,1)
plot(hE,HF - max(HF))
title('Free-energy','FontSize',16)
xlabel('Prior')
subplot(2,2,2)
plot(hE,Eh)
xlabel('Prior')
title('log-preciion','FontSize',16)

