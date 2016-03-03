function DEMO_DCM_PEB_FIT
% Test routine to check group DCM for electrophysiology
%--------------------------------------------------------------------------
% This routine illustrates the use of Bayesian model reduction when
% inverting hierarchical (dynamical) models; for example, multisubject DCM
% models. In this demonstration empirical Bayesian model reduction is
% applied recursively in an attempt to finesse the local  minimum problem
% with a nonlinear DCMs. The estimates are compared against standard
% Bayesian model reduction, in terms of the subject specific estimates and
% Bayesian model comparison  (and averages) at the between subject level.
%
% This demo considers a single goup (e.g., of normal subjects) and the
% differences between the group average using emprical Bayesian reduction
% and the Bayesian reduction of the (grand) average response.
%
% See also: DEMO_DCM_PEB_REC.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_DCM_PEB_FIT.m 6737 2016-03-03 12:05:51Z karl $


% change to directory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
% try
%     cd(fullfile(spm('Dir'),'tests','data','MEEG'))
% catch
%     cd('C:\Users\karl\Documents\SPM\DCM tests')
% end
close all, clear all
clc
rng('default')
MODEL = 'ERP';


% set up
%==========================================================================
load DCM_MMN                               % base DCM

DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'ERP';
DCM.options.model    = MODEL;
DCM.options.Nmax     = 64;
DCM.options.DATA     = 1;
DCM.options.Nmodes   = 8;
DCM.name             = 'DCM_GROUP';

% model space - within subject effects
%--------------------------------------------------------------------------
b{1}  = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0]; % intrinsic lower
b{2}  = [0 0 0 0 0;
    0 0 0 0 0;
    0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 0]; % intrinsic higher
b{3}  = [0 0 0 0 0;
    0 0 0 0 0;
    1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 0 0]; % forward
b{4}  = [0 0 1 0 0;
    0 0 0 1 0;
    0 0 0 0 0;
    0 0 0 0 0;
    0 0 0 0 0]; % and backward

k     = spm_perm_mtx(length(b));
for i = 1:size(k,1)
    B{i}  = zeros(size(b{1}));
    for j = 1:size(k,2)
        B{i} = B{i} + k(i,j)*b{j};
    end
end

% model space
%--------------------------------------------------------------------------
switch MODEL
    case{'ERP'}
        mw   = 2;                      % true model (within)
        Nm   = length(B);              % number of models
        Ns   = 8;                      % number of subjects
        C    = 16;                     % within:between [co]variance ratio
        Se   = 8;                      % s.d. of noise
        
        DCM.options.onset = 80;
        DCM.Ep.B{1} = [ ...
            0.0900    0         0         0         0
            0    0.1316         0         0         0
            0.1847    0    0.7691         0         0
            0    0.2937         0    0.2906         0
            0         0         0         0         0];
        
    case{'CMC'}
        mw   = 4;                      % true model (within)
        Nm   = length(B);              % number of models
        Ns   = 8;                      % number of subjects
        C    = 16;                     % within:between [co]variance ratio
        Se   = 8;                      % s.d. of noise
        
        DCM.A{3}   = DCM.A{3}*0;       % preclude modulatory effects
        DCM.options.onset = 100;
        
    otherwise
        
end


% invert empirical data
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    DCM = rmfield(DCM,'M');
end
DCM.B   = B(mw);
DCM     = spm_dcm_erp(DCM);

% create subject-specifc DCM
%==========================================================================
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
    Pp  = spm_vec(DCM.Ep) + spm_sqrtm(Cp)*randn(Np,1);
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
        e    = e*mean(std(y{c}))/mean(std(e))/Se;
        y{c} = y{c} + e;
        y{c} = DCM.M.R*y{c};
    end
    
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

% one model
%--------------------------------------------------------------------------
switch MODEL
    case{'CMC'}
        GCM = GCM(:,mw);
    otherwise
        GCM = GCM(:,1);
end

% creative grand average
%--------------------------------------------------------------------------
clear y Tp Tg
gcm   = GCM{1};
for i = 1:Ns
    y(:,i)  = spm_vec(GCM{i,1}.xY.y);
    Tp(:,i) = spm_vec(GCM{i,1}.Tp);
    Tg(:,i) = spm_vec(GCM{i,1}.Tg);
end

% principal eigenvariates
%--------------------------------------------------------------------------
gcm.xY.y = spm_unvec(mean(y,2),gcm.xY.y);
gcm.Tp   = spm_unvec(mean(Tp,2),gcm.Tp);
gcm.Tg   = spm_unvec(mean(Tg,2),gcm.Tg);

% indices for plotting
%--------------------------------------------------------------------------
iAB   = spm_find_pC(GCM{1},[],{'A','B'});
iA    = spm_find_pC(GCM{1},[],{'A'});
iB    = spm_find_pC(GCM{1},[],{'B'});


% The following section contains the key analyses
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% inversion of grand average
%==========================================================================
gcm = spm_dcm_fit(gcm);

% Bayesian model comparison and average
%--------------------------------------------------------------------------
[rcm,bmr,bma] = spm_dcm_bmr(gcm,{'A','B'});

% overlay true values and change graphic type
%--------------------------------------------------------------------------
subplot(3,2,3),hold on, bar(Tp(iAB),1/2), hold off
subplot(3,2,4),hold on, bar(Tp(iAB),1/2), hold off
set(gcf,'Tag','BMR - grand average','name','BMR - grand average');


% empirical Bayesian group inversion
%==========================================================================
for i = 1:Ns
    GCM{i,1}.M.dipfit = DCM.M.dipfit;
end

[RCM,REB,M,HCM] = spm_dcm_peb_fit(GCM);
HCM    = spm_dcm_reduce(HCM,M.bE,M.bC);


% BMC/BMA - second level
%==========================================================================

% Empirical Bayes and model comparison (specified model space)
%--------------------------------------------------------------------------
[PEB,PCM] = spm_dcm_peb(HCM(:,2),[],{'A','B'});
BMA       = spm_dcm_peb_bmc(PEB);

% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','BMR - all');
Tp  = spm_vec(DCM.Ep);
subplot(3,2,3),hold on, bar(Tp(iAB),1/2), hold off
subplot(3,2,4),hold on, bar(Tp(iAB),1/2), hold off



%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% extract and plot results over iterations
%==========================================================================
clear Q T
for i = 1:Ns
    
    % Parameter estimates
    %----------------------------------------------------------------------
    T(:,i) = spm_vec(GCM{i,1}.Tp);
    for j = 1:size(HCM,2)
        Q(:,i,j) = spm_vec(HCM{i,j}.Ep);
    end
    
end
i     = iB;
for j = 1:size(HCM,2)
    R(j) = corr(spm_vec(T(i,:)),spm_vec(Q(i,:,j)));
end

% free energy over iterations of empirical Bayes
%--------------------------------------------------------------------------
spm_figure('GetWin','Inversion'); j = 1:4;

subplot(3,2,1), bar(RCM{1}.FEB(j))
xlabel('iteration'), ylabel('(second level) free energy')
title('Free energy','FontSize',16)
axis square

subplot(3,2,2), bar(R(j))
xlabel('iteration'), ylabel('correlation')
title('Correlation with truth','FontSize',16)
axis square

subplot(3,2,3), bar(RCM{1}.EEB(j),'b')
xlabel('iteration'), ylabel('correlation')
title('log precision','FontSize',16)
axis square

subplot(3,2,4), bar(RCM{1}.HEB(j),'c')
xlabel('iteration'), ylabel('correlation')
title('Posterior uncertainty','FontSize',16)
axis square


% extract and plot results
%==========================================================================
clear Q
for i = 1:Ns
    
    % data - over subjects
    %----------------------------------------------------------------------
    Y(:,i,1) = GCM{i,1}.xY.y{1}*DCM.M.U(:,1);
    Y(:,i,2) = GCM{i,1}.xY.y{2}*DCM.M.U(:,1);
    
    % Parameter estimates
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Tp);
    Q(:,i,2) = spm_vec(PCM{i,1}.Ep);
    
end


% first level parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','parameters');clf, ALim = 2;

% indices to plot parameters
%--------------------------------------------------------------------------
Qp   = spm_vec(bma.Ep);
r    = corr(spm_vec(Q(iB,:,1)),spm_vec(Q(iB,:,2)));
tstr = sprintf('Empirical Bayes: cor = %-0.2f',r);

subplot(2,1,1)
plot(Q(iA,:,1),Q(iA,:,2),'.c','MarkerSize',12), hold on
plot(Q(iB,:,1),Q(iB,:,2),'.b','MarkerSize',12), hold on
plot([-1 1],[-1 1],'-.'), hold off
xlabel('true parameter'), ylabel('estimate')
title(tstr,'FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

% group means
%--------------------------------------------------------------------------
jA  = 1:length(iA);
jB  = (1:length(iB)) + jA(end);

subplot(2,2,3)
plot(Tp(iAB),BMA.Ep,  '.k','MarkerSize',16), hold on
plot(Tp(iAB),Qp(iAB), '.b','MarkerSize',16), hold on
plot([-1 1],[-1 1],'-.'), hold off
xlabel('true parameter'), ylabel('estimate')
title(' Bayesian estimators','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square

subplot(2,2,4)
plot(Qp(iA),BMA.Ep(jA),'.c','MarkerSize',16), hold on
plot(Qp(iB),BMA.Ep(jB),'.b','MarkerSize',16), hold on
plot([-1 1],[-1 1],'-.'), hold off
xlabel('Grand average estimate'), ylabel('empirical Bayesian estimate')
title('Empirical and grand average','FontSize',16)
axis([-1 1 -1 1]*ALim), axis square


% plot simulation data
%==========================================================================
spm_figure('GetWin','Data');clf

pst = DCM.xY.pst;
g(:,1) = gcm{1}.xY.y{1}*DCM.M.U(:,1);
g(:,2) = gcm{1}.xY.y{2}*DCM.M.U(:,1);

subplot(2,2,1)
plot(pst,DCM.M.R*x{2}*G','k'), hold on
plot(pst,x{2}*G',':k'),        hold off
xlabel('pst (ms)'), ylabel('response'), title('Signal (single subject)','FontSize',16)
axis square, spm_axis tight,  a = axis;

subplot(2,2,2)
plot(pst,DCM.M.R*e,'k'), hold on
plot(pst,e,':k'),        hold off
xlabel('pst (ms)'), ylabel('response'), title('Noise','FontSize',16)
axis square, spm_axis tight, axis(a)

subplot(2,2,3)
plot(pst,Y(:,:,1), 'k'), hold on
plot(pst,Y(:,:,2),':k'), hold on
plot(pst,g(:,1),'r','LineWidth',2), hold on
plot(pst,g(:,2),'r','LineWidth',2), hold on
xlabel('pst (ms)'), ylabel('response'), title('Group data','FontSize',16)
axis square, spm_axis tight, a = axis;

subplot(2,2,4)
plot(pst,Y(:,:,2) - Y(:,:,1),'k'), hold on
plot(pst,g(:,2) - g(:,1),'r','LineWidth',2), hold on
xlabel('pst (ms)'), ylabel('differential response'), title('Difference waveforms','FontSize',16)
axis square, spm_axis tight, axis(a)

