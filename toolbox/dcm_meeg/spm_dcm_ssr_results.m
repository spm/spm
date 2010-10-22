function [DCM] = spm_dcm_ssr_results(DCM,Action)
% Results for ERP Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_erp_results(DCM,'spectral data');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (A)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (B)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (C)');
% FORMAT spm_dcm_erp_results(DCM,'trial-specific effects');
% FORMAT spm_dcm_erp_results(DCM,'Input');
% FORMAT spm_dcm_erp_results(DCM,'Cross-spectral density');
% FORMAT spm_dcm_erp_results(DCM,'Dipoles');
%                
%___________________________________________________________________________
%
% DCM is a causal modelling procedure for dynamical systems in which
% causality is inherent in the differential equations that specify the model.
% The basic idea is to treat the system of interest, in this case the brain,
% as an input-state-output system.  By perturbing the system with known
% inputs, measured responses are used to estimate various parameters that
% govern the evolution of brain states.  Although there are no restrictions
% on the parameterisation of the model, a bilinear approximation affords a
% simple re-parameterisation in terms of effective connectivity.  This
% effective connectivity can be latent or intrinsic or, through bilinear
% terms, model input-dependent changes in effective connectivity.  Parameter
% estimation proceeds using fairly standard approaches to system
% identification that rest upon Bayesian inference.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ssr_results.m 4096 2010-10-22 19:40:34Z karl $
 
 
% get figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf

% placespectral features in xY.y
%--------------------------------------------------------------------------
DCM.xY.y  = spm_cond_units(DCM.xY.csd,'csd');

% trial data
%--------------------------------------------------------------------------
xY  = DCM.xY;                   % data
nt  = length(xY.y);             % Nr trial types
nf  = size(xY.y{1},1);          % Nr frequency bins
nm  = size(xY.y{1},2);          % Nr spatial modes
Hz  = xY.Hz;                    % PST

% switch
%--------------------------------------------------------------------------
switch(lower(Action))    
    
case{lower('spectral data')}
    
    % spm_dcm_ssr_results(DCM,'Data');
    %----------------------------------------------------------------------
    co = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz = xY.Hz;
    q  = max(spm_vec(xY.y));
    nm = min(nm,4);
    
    for k = 1:nt
        str{k} = sprintf('trial %i',k);
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j),cla
            for k = 1:nt
                plot(Hz,xY.y{k}(:,i,j),'color',co{k}), hold on
                set(gca,'YLim',[0 q])
            end
        end
 
        % spectral density
        %------------------------------------------------------------------
        subplot(2,2,3)
        for k = 1:nt
            plot(Hz,xY.y{k}(:,i,i),'color',co{i}), hold on
            set(gca,'YLim',[0 q])
        end
    end
    
    title('spectral density over modes')
    xlabel('Frequency (Hz)')
    ylabel('CSD')
    axis square
    return
    
end
 
% post inversion parameters
%--------------------------------------------------------------------------
nu  = length(DCM.B);          % Nr experimental inputs
ns  = size(DCM.A{1},2);       % Nr of sources

 
% switch
%--------------------------------------------------------------------------
switch(lower(Action))    
    
case{lower('Coupling (A)')}
    
    % spm_dcm_ssr_results(DCM,'coupling (A)');
    %----------------------------------------------------------------------
    str = {'Forward','Backward','Lateral'};
    for  i = 1:3
        
        % images
        %------------------------------------------------------------------
        subplot(4,3,i)
        imagesc(exp(DCM.Ep.A{i}))
        title(str{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square
    
        % table
        %------------------------------------------------------------------
        subplot(4,3,i + 3)
        text(0,1/2,num2str(full(exp(DCM.Ep.A{i})),' %.2f'),'FontSize',8)
        axis off,axis square
 
    
        % PPM
        %------------------------------------------------------------------
        subplot(4,3,i + 6)
        image(64*DCM.Pp.A{i})
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        title('PPM')
        axis square
    
        % table
        %------------------------------------------------------------------
        subplot(4,3,i + 9)
        text(0,1/2,num2str(DCM.Pp.A{i},' %.2f'),'FontSize',8)
        axis off, axis square
        
    end
    
case{lower('Coupling (C)')}
    
    % spm_dcm_ssr_results(DCM,'coupling (C)');
    %----------------------------------------------------------------------
    
    % images
    %----------------------------------------------------------------------
    subplot(2,4,1)
    imagesc(exp(DCM.Ep.C))
    title('Factors','FontSize',10)
    set(gca,'XTick',[1:nu],'XTickLabel','Input','FontSize',8)
    set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname, 'FontSize',8)
    axis square
    
    % PPM
    %----------------------------------------------------------------------
    subplot(2,4,3)
    image(64*DCM.Pp.C)
    title('Factors','FontSize',10)
    set(gca,'XTick',[1:nu],'XTickLabel','Input','FontSize',8)
    set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname, 'FontSize',8)
    axis square
    title('PPM')
    
    % table
    %----------------------------------------------------------------------
    subplot(2,4,2)
    text(0,1/2,num2str(full(exp(DCM.Ep.C)),' %.2f'),'FontSize',8)
    axis off
 
    % table
    %----------------------------------------------------------------------
    subplot(2,4,4)
    text(0,1/2,num2str(DCM.Pp.C,' %.2f'),'FontSize',8)
    axis off
 
 
case{lower('Coupling (B)')}
    
    % spm_dcm_ssr_results(DCM,'coupling (B)');
    %----------------------------------------------------------------------
    for i = 1:nu
        
        % images
        %------------------------------------------------------------------
        subplot(4,nu,i)
        imagesc(exp(DCM.Ep.B{i}))
        title(DCM.xU.name{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square
 
        % tables
        %------------------------------------------------------------------
        subplot(4,nu,i + nu)
        text(0,1/2,num2str(full(exp(DCM.Ep.B{i})),' %.2f'),'FontSize',8)
        axis off
        axis square
        
        % PPM
        %------------------------------------------------------------------
        subplot(4,nu,i + 2*nu)
        image(64*DCM.Pp.B{i})
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        title('PPM')
        axis square
 
        % tables
        %------------------------------------------------------------------
        subplot(4,nu,i + 3*nu)
        text(0,1/2,num2str(DCM.Pp.B{i},' %.2f'),'FontSize',8)
        axis off
        axis square
        
    end
    
case{lower('trial-specific effects')}
    
    % spm_dcm_ssr_results(DCM,'trial-specific effects');
    %----------------------------------------------------------------------
    for i = 1:ns
        for j = 1:ns
 
            % ensure connection is enabled
            %--------------------------------------------------------------
            q     = 0;
            for k = 1:nu
                q = q | DCM.B{k}(i,j);
            end
 
            % plot trial-specific effects
            %--------------------------------------------------------------
            if q
                B     = zeros(nt,1);
                for k = 1:nu
                    B = B + DCM.xU.X(:,k)*DCM.Ep.B{k}(i,j);
                end
                
                subplot(ns,ns,(i - 1)*ns + j)
                bar(exp(B)*100,'c')
                title([DCM.Sname{j}, ' to ' DCM.Sname{i}],'FontSize',10)
                xlabel('trial',  'FontSize',8)
                ylabel('strength (%)','FontSize',8)
                set(gca,'XLim',[0 nt + 1])
                axis square
 
            end
        end
    end
    
case{lower('Input')}
    
    % spectrum of innovations or noise (Gu)
    %----------------------------------------------------------------------
    try
        Gu   = exp(DCM.Ep.a)*xY.Hz.^(-1)*2;    % spectral density of (AR) input
        Gu   = Gu + exp(DCM.Ep.b);             % spectral density of IID input
    catch
        Gu   = exp(DCM.Ep.a(1))*xY.Hz.^(-1);    % spectral density of (AR) input
        Gu   = Gu + exp(DCM.Ep.a(2));           % spectral density of IID input
    end
    
    % plot spectral density of innovations
    % ---------------------------------------------------------------------
    subplot(2,1,1)
    plot(xY.Hz,Gu)
    xlabel('frquency (Hz)')
    title('spectrum of innovations or noise')
    axis square, grid on
    
case{lower('Cross-spectral density')}
    
    % spm_dcm_ssr_results(DCM,'Cross-spectral density');
    %----------------------------------------------------------------------
    co = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz = xY.Hz;
    q  = max(spm_vec(DCM.Hc));
    nm = min(nm,4);
    
    tstr = {};
    mstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('predicted: trial %i',k);
        tstr{end + 1} = sprintf('observed: trial %i',k);
    end
    for k = 1:nm
        mstr{end + 1} = sprintf('predicted: mode %i',k);
        mstr{end + 1} = sprintf('observed: mode %i',k);
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j),cla
            for k = 1:nt
                plot(Hz,DCM.Hc{k}(:,i,j),'color',co{k}), hold on
                plot(Hz,DCM.Hc{k}(:,i,j) + DCM.Rc{k}(:,i,j),':','color',co{k})
                set(gca,'YLim',[0 q])
            end
        end
 
        % legend
        %------------------------------------------------------------------      
        if i == nm && j == nm
            legend(tstr)
        end
        
        % spectral density
        %------------------------------------------------------------------
        subplot(2,2,3)
        for k = 1:nt
            plot(Hz,DCM.Hc{k}(:,i,i),'color',co{i}), hold on
            plot(Hz,DCM.Hc{k}(:,i,i) + DCM.Rc{k}(:,i,i),':','color',co{i})
            set(gca,'YLim',[0 q])
        end
    end
   
    title({'Spectral density over modes';'(in channel-space)'},'FontSize',16)
    xlabel('Frequency (Hz)')
    ylabel('root CSD')
    axis square
    legend(mstr)
    
    
    
case{lower('Dipoles')}
    
    % return if LFP
    % ---------------------------------------------------------------------
    if strcmp(lower(DCM.xY.modality),'lfp')
        warndlg('There are no ECDs for these LFP data')
        return
    end
    
    % plot dipoles
    % ---------------------------------------------------------------------
    try
        P            = DCM.Ep;   
        np           = size(P.L,2)/size(P.Lpos,2);
        sdip.n_seeds = 1;
        sdip.n_dip   = np*ns;
        sdip.Mtb     = 1;
        sdip.j{1}    = full(P.L);
        sdip.j{1}    = sdip.j{1}(:);
        sdip.loc{1}  = kron(ones(1,np),full(P.Lpos));
        spm_eeg_inv_ecd_DrawDip('Init', sdip)
    end
 
end
