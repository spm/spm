function [DCM] = spm_dcm_tfm_results(DCM,Action,fig)
% Results for CSD (SSR) Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_tfm_results(DCM,'induced responses');
% FORMAT spm_dcm_tfm_results(DCM,'induced and evoked responses');
% FORMAT spm_dcm_tfm_results(DCM,'Coupling (A)');
% FORMAT spm_dcm_tfm_results(DCM,'Coupling (B)');
% FORMAT spm_dcm_tfm_results(DCM,'Coupling (C)');
% FORMAT spm_dcm_tfm_results(DCM,'trial-specific effects');
% FORMAT spm_dcm_tfm_results(DCM,'Endogenous input');
% FORMAT spm_dcm_tfm_results(DCM,'Exogenous input');
% FORMAT spm_dcm_tfm_results(DCM,'Transfer functions');
% FORMAT spm_dcm_tfm_results(DCM,'induced predictions')
% FORMAT spm_dcm_tfm_results(DCM,'induced and evoked predictions')
% FORMAT spm_dcm_tfm_results(DCM,'induced predictions - sources')
% FORMAT spm_dcm_tfm_results(DCM,'induced and evoked predictions - sources')
% FORMAT spm_dcm_tfm_results(DCM,'Dipoles');
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
% $Id: spm_dcm_tfm_results.m 6234 2014-10-12 09:59:10Z karl $
 
 
% get figure
%--------------------------------------------------------------------------
if nargin < 3, spm_figure('GetWin','Graphics'); end, colormap(gray), clf
 
% get action if neccessary
%--------------------------------------------------------------------------
if nargin == 1
    
    str    = {'induced responses',...
              'induced and evoked responses',...
              'Coupling (A)',...
              'Coupling (B)',...
              'Coupling (C)',...
              'trial-specific effects',...
              'Endogenous input',...
              'Exogenous input',...
              'Transfer functions',...
              'induced predictions',...
              'induced and evoked predictions',...
              'induced predictions - sources',...
              'induced and evoked predictions - sources',...
              'Dipoles'};
    
    s      = listdlg('PromptString','Select an option:',...
        'SelectionMode','single',...
        'ListString',str);
    Action = str{s};
end
 
% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
    case{lower('induced responses')}
        
        % induced responses (channels)
        % -----------------------------------------------------------------
        for c = 1:length(DCM.xY.csd)
            spm_figure('GetWin',sprintf('induced responses: condition %i',c));
            spm_dcm_tfm_image(DCM.xY.csd{c},DCM.xY.pst,DCM.xY.Hz)
        end
        
        
    case{lower('induced and evoked responses')}
        
        % induced and evoked responses (channels)
        % -----------------------------------------------------------------
        spm_figure('GetWin','induced and evoked responses');
        spm_dcm_tfm_response(DCM.xY,DCM.xY.pst,DCM.xY.Hz)
        
end
 
% dimensions
%--------------------------------------------------------------------------
nu  = length(DCM.B);          % Nr experimental inputs
ns  = size(DCM.A{1},2);       % Nr of sources
 
% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
    case{lower('Coupling (A)')}
        
        % spm_dcm_tfm_results(DCM,'coupling (A)');
        %------------------------------------------------------------------
        spm_figure('GetWin','Graphics');
        
        str = {'Forward (i)','Forward (ii)','Backward (i)','Backward (ii)'};
        m   = length(DCM.Ep.A);
        for  i = 1:m
            
            % images
            %--------------------------------------------------------------
            subplot(4,m,i)
            imagesc(exp(DCM.Ep.A{i}))
            title(str{i},'FontSize',10)
            set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
            set(gca,'XTick',[])
            xlabel('from','FontSize',8)
            ylabel('to','FontSize',8)
            axis square
            
            % table
            %--------------------------------------------------------------
            subplot(4,m,i + m)
            text(0,1/2,num2str(full(exp(DCM.Ep.A{i})),' %.2f'),'FontSize',8)
            axis off,axis square
            
            
            % PPM
            %--------------------------------------------------------------
            subplot(4,m,i + m + m)
            image(64*DCM.Pp.A{i})
            set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
            set(gca,'XTick',[])
            title('PPM')
            axis square
            
            % table
            %--------------------------------------------------------------
            subplot(4,m,i + m + m + m)
            text(0,1/2,num2str(DCM.Pp.A{i},' %.2f'),'FontSize',8)
            axis off, axis square
            
        end
        
    case{lower('Coupling (C)')}
        
        % spm_dcm_tfm_results(DCM,'coupling (C)');
        %------------------------------------------------------------------
        
        % images
        %------------------------------------------------------------------
        subplot(2,4,1)
        imagesc(exp(DCM.Ep.C))
        title('Factors','FontSize',10)
        set(gca,'XTick',1:nu,'XTickLabel','Input','FontSize',8)
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname, 'FontSize',8)
        axis square
        
        % PPM
        %-----------------------------------------------------------------
        subplot(2,4,3)
        image(64*DCM.Pp.C)
        title('Factors','FontSize',10)
        set(gca,'XTick',1:nu,'XTickLabel','Input','FontSize',8)
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname, 'FontSize',8)
        axis square
        title('PPM')
        
        % table
        %------------------------------------------------------------------
        subplot(2,4,2)
        text(0,1/2,num2str(full(exp(DCM.Ep.C)),' %.2f'),'FontSize',8)
        axis off
        
        % table
        %------------------------------------------------------------------
        subplot(2,4,4)
        text(0,1/2,num2str(DCM.Pp.C,' %.2f'),'FontSize',8)
        axis off
        
        
    case{lower('Coupling (B)')}
        
        % spm_dcm_tfm_results(DCM,'coupling (B)');
        %------------------------------------------------------------------
        for i = 1:nu
            
            % images
            %--------------------------------------------------------------
            subplot(4,nu,i)
            imagesc(exp(DCM.Ep.B{i}))
            title(DCM.xU.name{i},'FontSize',10)
            set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
            set(gca,'XTick',[])
            xlabel('from','FontSize',8)
            ylabel('to','FontSize',8)
            axis square
            
            % tables
            %--------------------------------------------------------------
            subplot(4,nu,i + nu)
            text(0,1/2,num2str(full(exp(DCM.Ep.B{i})),' %.2f'),'FontSize',8)
            axis off
            axis square
            
            % PPM
            %--------------------------------------------------------------
            subplot(4,nu,i + 2*nu)
            image(64*DCM.Pp.B{i})
            set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
            set(gca,'XTick',[])
            title('PPM')
            axis square
            
            % tables
            %--------------------------------------------------------------
            subplot(4,nu,i + 3*nu)
            text(0,1/2,num2str(DCM.Pp.B{i},' %.2f'),'FontSize',8)
            axis off
            axis square
            
        end
        
    case{lower('trial-specific effects')}
        
        % spm_dcm_tfm_results(DCM,'trial-specific effects');
        %------------------------------------------------------------------
        nt    = size(DCM.xU.X,1);
        for i = 1:ns
            for j = 1:ns
                
                % ensure connection is enabled
                %----------------------------------------------------------
                q     = 0;
                for k = 1:nu
                    q = q | DCM.B{k}(i,j);
                end
                
                % plot trial-specific effects
                %----------------------------------------------------------
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
        
    case{lower('Endogenous input')}
        
        % spectrum of innovations or noise (Gu)
        %------------------------------------------------------------------
        [Gu,Gs,Gn,Hz] = spm_csd_mtf_gu(DCM.Ep,DCM.M);
        
        
        % plot spectral density of innovations
        % -----------------------------------------------------------------
        subplot(2,1,1)
        plot(Hz,Gu)
        xlabel('frequency (Hz)')
        title('Spectrum of innovations','FontSize',16)
        axis square, grid on
        legend(DCM.Sname)
        
        % plot spectral density of noise
        % -----------------------------------------------------------------
        subplot(2,2,3)
        plot(Hz,Gs)
        xlabel('frequency (Hz)')
        title('Channel-specific noise')
        axis square, grid on
        
        % plot spectral density of noise
        % -----------------------------------------------------------------
        subplot(2,2,4)
        plot(Hz,Gn)
        xlabel('frequency (Hz)')
        title('Non-specific noise')
        axis square, grid on
        
    case{lower('Exogenous input')}
        
        % plot data
        % -----------------------------------------------------------------
        pst   = DCM.pst - DCM.pst(1);
        U     = spm_erp_u(pst,DCM.Ep,DCM.M);
        
        subplot(2,1,1)
        plot(DCM.xY.pst,U)
        xlabel('time (ms)')
        title('input')
        axis square, grid on
        for i = 1:length(DCM.M.ons)
            str{i} = sprintf('input (%i)',i);
        end
        legend(str)
        
        
    case{lower('Transfer functions')}
        
        for c = 1:length(DCM.DTF)
            str = sprintf('Directed transfer functions (among sources): condition %i',c);
            spm_figure('GetWin',str);
            spm_dcm_tfm_transfer(DCM.DTF{c},DCM.xY.pst,DCM.xY.Hz)
        end
        
        
    case{lower('induced predictions')}
        
        % induced predictions (channels)
        % -----------------------------------------------------------------
        for c = 1:length(DCM.xY.csd)
            str = sprintf('induced predictions (among channels): condition %i',c);
            spm_figure('GetWin',str);
            spm_dcm_tfm_image(DCM.csd{c},DCM.xY.pst,DCM.xY.Hz)
        end
        
        
    case{lower('induced and evoked predictions')}
        
        % induced and evoked predictions (channels)
        % -----------------------------------------------------------------
        spm_figure('GetWin','induced and evoked predictions');
        xY.csd = DCM.csd;
        xY.erp = DCM.erp;
        spm_dcm_tfm_response(xY,DCM.xY.pst,DCM.xY.Hz,DCM.xY)
        
    case{lower('induced predictions - sources')}
        
        % induced predictions (sources)
        % -----------------------------------------------------------------
        for c = 1:length(DCM.xY.csd)
            str = sprintf('induced predictions (among sources): condition %i',c);
            spm_figure('GetWin',str);
            spm_dcm_tfm_image(DCM.CSD{c},DCM.xY.pst,DCM.xY.Hz)
        end
        
        
    case{lower('induced and evoked predictions - sources')}
        
        % induced and evoked predictions (sources)
        % -----------------------------------------------------------------
        str    = 'induced and evoked predictions (among sources)';
        spm_figure('GetWin',str);
        xY.csd = DCM.CSD;
        xY.erp = DCM.ERP;
        spm_dcm_tfm_response(xY,DCM.xY.pst,DCM.xY.Hz)
        
        
    case{lower('Dipoles')}
        
        % return if LFP
        % -----------------------------------------------------------------
        if strcmpi(DCM.xY.modality,'lfp')
            warndlg('There are no ECDs for these LFP data')
            return
        end
        
        % plot dipoles
        % -----------------------------------------------------------------
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
