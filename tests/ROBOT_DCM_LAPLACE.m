function E = ROBOT_DCM_LAPLACE
% test routine to check group DCM for electrophysiology
%==========================================================================
%   options.analysis     - 'ERP','CSD', 'IND' or 'TFM
%   options.model        - 'ERP','SEP','CMC','LFP','NNM' or 'MFM'
%   options.spatial      - 'ECD','LFP' or 'IMG'

% $Id: ROBOT_DCM_LAPLACE.m 6296 2015-01-02 16:20:19Z guillaume $

% tests of spatial models: 'ECD', 'LFP' or 'IMG'
%==========================================================================
try
    cd('/home/spm/tests/DCM/MEEG')
catch
    cd('C:\Users\karl\Documents\SPM\DCM tests')
end
close all
delete(get(0,'Children'))
rng('default')

E = {};

% set up
%==========================================================================
load DCM_MMN                             % base DCM

DCM.options.spatial  = 'ECD';
DCM.options.analysis = 'ERP';
DCM.options.model    = 'ERP';
DCM.options.Nmax     = 32;
DCM.options.DATA     = 1;
DCM.name             = 'DCM_GROUP';

% model space
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
Pm  = 3;                                  % true model
Nm  = length(B);                          % number of models
Ns  = 8;                                  % number of subjects
C   = 16;                                 % within:between [co]variance

% invert base model
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    DCM = rmfield(DCM,'M');
end
DCM.B = B(Pm);
DCM   = spm_dcm_erp(DCM);



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
    Pp  = spm_vec(DCM.Ep) + spm_sqrtm(Cp)*randn(Np,1);
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
        L(i,j,:) = GCM{i,j}.L;
        R(i,j)   = RCM{i,j}.F;
    end
    
end

% select parameters
%--------------------------------------------------------------------------
i     = spm_fieldindices(Pp,'A','B');
c     = spm_vec(DCM.M.pC);
iq    = i(find(c(i)));
Q     = Q(iq,:,:);


% relative free energy
%--------------------------------------------------------------------------
F = F - F(:,1)*ones(1,Ns);


% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Data');clf

subplot(2,2,1), plot(Y(:,:,1)), hold on
plot(Y(:,:,2),'-.'), hold off
xlabel('pst'), ylabel('response'), title('Data','FontSize',16)
axis square

subplot(2,2,2), plot(Y(:,:,1) - Y(:,:,2))
xlabel('pst'), ylabel('condition effects'), title('Data','FontSize',16)
axis square


% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');clf

subplot(3,2,1), plot(Q(:,:,3),Q(:,:,1),'o'), hold on
xlabel('true parameter'), ylabel('Model average'), title('Parameters','FontSize',16)
axis square

subplot(3,2,2), plot(Q(:,:,3),Q(:,:,2),'o'), hold on
xlabel('true parameter'), ylabel('BMA (BMR)'), title('Parameters','FontSize',16)
axis square

subplot(3,2,3), imagesc(F)
xlabel('model'), ylabel('subject'), title('Free energy','FontSize',16)
axis square

subplot(3,2,4), imagesc(R)
xlabel('model'), ylabel('subject'), title('Free energy (BMR)','FontSize',16)
axis square

subplot(3,2,5)
GF  = sum(F);
GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('Free energy'), title('Free energy (FFX)','FontSize',16)
spm_axis tight, axis square

subplot(3,2,6)
GF  = sum(R);
GF  = GF - min(GF);
bar(GF), xlabel('model'), ylabel('FFX Free energy'), title('Free energy (BMR)','FontSize',16)
spm_axis tight, axis square


% reinvert (full) model with initialization; recursively
%==========================================================================
clear COR IF cor f
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
