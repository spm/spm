% function DEMO_DCM_PEB_REC
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
% First, an EEG DCM (using empirical grand mean data) is inverted to
% find plausible group mean parameters. Single subject data are
% then generated using typical within and between subject variance (here, 
% group differences in the modulation of intrinsic connectivity. We then
% illustrate empirical Bayesian model reduction in a recursive mode.
%
% See also: DEMO_DCM_PEB.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Peter Zeidman
% $Id: DEMO_DCM_PEB_REC.m 6737 2016-03-03 12:05:51Z karl $


% change to directory with empirical data
%--------------------------------------------------------------------------
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'
%--------------------------------------------------------------------------
try
    cd(fullfile(spm('Dir'),'tests','data','MEEG'))
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
Ep      = DCM.Ep;

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
Tp          = spm_vec(DCM.Ep);           % true second level paramters
Tx          = spm_vec(DCM.Ex);           % true second level paramters

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

% The following section contains the key analyses
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


% second level model
%--------------------------------------------------------------------------
M     = struct('X',X);

% invert rreduced models (standard inversion)
%==========================================================================
GCM   = spm_dcm_fit(GCM);

% Bayesian model reduction (avoiding local minima over models)
%==========================================================================
RCM   = spm_dcm_bmr(GCM);

% hierarchical (empirical Bayes) model reduction
%==========================================================================
[peb,PCM] = spm_dcm_peb(RCM,[],'all');

% alternative (more robust but expensive) iterative inversion
%--------------------------------------------------------------------------
for i = 1:Ns
    GCM{i,1}.M.dipfit = DCM.M.dipfit;
end
ECM   = spm_dcm_peb_fit(GCM);


spm_figure('GetWin','SI');
subplot(3,2,6), bar(ECM{1}.FEB)
xlabel('recursion'), ylabel(' (second level) free energy')
title('Free energy','FontSize',16)
axis square


% BMC/BMA - second level
%==========================================================================

% BMC - search over first and second level effects
%--------------------------------------------------------------------------
[BMC,PEB] = spm_dcm_bmc_peb(PCM,M,{'A','B'});
set(gcf,'Name',[get(gcf,'Name') ' - standard' ]);
set(gcf,'Tag' ,[get(gcf,'Tag' ) ' - standard' ]);

[EMC,EEB] = spm_dcm_bmc_peb(ECM,M,{'A','B'});
set(gcf,'Name',[get(gcf,'Name') ' - recursive']);
set(gcf,'Tag' ,[get(gcf,'Tag' ) ' - recursive' ]);



%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

% extract and plot results
%==========================================================================
clear Q
for i = 1:Ns
        
    % Parameter estimates
    %----------------------------------------------------------------------
    Q(:,i,1) = spm_vec(GCM{i,1}.Tp);
    Q(:,i,2) = spm_vec(GCM{i,1}.Ep);
    Q(:,i,3) = spm_vec(RCM{i,1}.Ep);
    Q(:,i,4) = spm_vec(PCM{i,1}.Ep);
    Q(:,i,5) = spm_vec(ECM{i,1}.Ep);
    
    % Free energies
    %----------------------------------------------------------------------
    for j = 1:Nm
        F(i,j,1) = GCM{i,j}.F - GCM{i,1}.F;
        F(i,j,2) = RCM{i,j}.F - RCM{i,1}.F;
        F(i,j,3) = PCM{i,j}.F - PCM{i,1}.F;
        F(i,j,4) = ECM{i,j}.F - ECM{i,1}.F;
    end
    
end



% first level parameter estimates and Bayesian model averages
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf, ALim = 2/3;

% indices to plot parameters
%--------------------------------------------------------------------------
pC    = GCM{1,1}.M.pC;
iA    = spm_find_pC(pC,pC,'A');
iB    = spm_find_pC(pC,pC,'B');
str   = {'FFX','BMR','PEB','RFX'};
for i = 1:4
    
    r    = corr(spm_vec(Q([iB],:,1)),spm_vec(Q([iB],:,i + 1)));
    tstr = sprintf('%s: cor = %-0.2f',str{i},r);
    subplot(4,2,(i - 1)*2 + 1)
    plot(Q(iA,:,1),Q(iA,:,i + 1),'.c','MarkerSize',12), hold on
    plot(Q(iB,:,1),Q(iB,:,i + 1),'.b','MarkerSize',12), hold off
    xlabel('true parameter'), ylabel('Model average')
    title(tstr,'FontSize',16)
    axis([-1 1 -1 1]*ALim), axis square
    
    p   = spm_softmax(sum(F(:,:,i))');
    subplot(4,2,(i - 1)*2 + 2)
    bar(p),[m,j] = max(p); 
    text(j - 1/4,m/2,sprintf('%-2.0f%%',m*100),'Color','w','FontSize',8)
    xlabel('model'), ylabel('probability')
    title(['Posterior ('  str{i} ')'],'FontSize',16)
    axis([0 (length(p) + 1) 0 1]), axis square

end
