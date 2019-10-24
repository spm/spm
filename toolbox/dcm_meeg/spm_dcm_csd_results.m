function spm_dcm_csd_results(DCM,Action,fig)
% Results for CSD (SSR) Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_csd_results(DCM,'spectral data');
% FORMAT spm_dcm_csd_results(DCM,'Coupling (A)');
% FORMAT spm_dcm_csd_results(DCM,'Coupling (B)');
% FORMAT spm_dcm_csd_results(DCM,'Coupling (C)');
% FORMAT spm_dcm_csd_results(DCM,'trial-specific effects');
% FORMAT spm_dcm_csd_results(DCM,'Input');
% FORMAT spm_dcm_csd_results(DCM,'Transfer functions');
% FORMAT spm_dcm_csd_results(DCM,'Cross-spectra (sources)')
% FORMAT spm_dcm_csd_results(DCM,'Cross-spectra (channels)')
% FORMAT spm_dcm_csd_results(DCM,'Coherence (sources)')
% FORMAT spm_dcm_csd_results(DCM,'Coherence (channels)')
% FORMAT spm_dcm_csd_results(DCM,'Covariance (sources)')
% FORMAT spm_dcm_csd_results(DCM,'Covariance (channels)')
% FORMAT spm_dcm_csd_results(DCM,'Dipoles');
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
% $Id: spm_dcm_csd_results.m 7679 2019-10-24 15:54:07Z spm $
 
 
% get figure
%--------------------------------------------------------------------------
if nargin < 3
    spm_figure('GetWin','Graphics');
    colormap(gray), clf
else
    figure(fig)
end

 
% get action if neccessary
%--------------------------------------------------------------------------
if nargin == 1
    str    = {'spectral data',...
              'Coupling (A)',...
              'Coupling (B)',...
              'Coupling (C)',...
              'trial-specific effects',...
              'Input',...
              'Transfer functions',...
              'Cross-spectra (sources)',...
              'Cross-spectra (channels)',...
              'Coherence (sources)',...
              'Coherence (channels)',...
              'Covariance (sources)',...
              'Covariance (channels)',...
              'Dipoles'};
          
    s      = listdlg('PromptString','Select an option:',...
                     'SelectionMode','single',...
                     'ListString',str);
    Action = str{s};
end
 
 
% place spectral features in xY.y
%--------------------------------------------------------------------------
if ischar(DCM), load(DCM); end
DCM.xY.y  = DCM.xY.csd;
 
% trial data
%--------------------------------------------------------------------------
xY  = DCM.xY;                   % data
nt  = length(xY.y);             % Nr trial types
nm  = size(xY.y{1},2);          % Nr spatial modes
Hz  = xY.Hz;                    % PST
 
% switch
%--------------------------------------------------------------------------
switch(lower(Action))    
    
case{lower('spectral data')}
    
    % spm_dcm_csd_results(DCM,'Data');
    %----------------------------------------------------------------------
    co    = {'b', 'r', 'g', 'm', 'y', 'k', 'c', 'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    lt    = {'-', ':', '-.', '-', ':', '-.','-', ':', '-.','-', ':', '-.'};
    Hz    = xY.Hz;
    nm    = min(nm,4);
    try
        name = DCM.xY.name;
    catch
        name = DCM.Sname;
    end
    
    tstr = {};
    mstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('trial %i: real', k);
        tstr{end + 1} = sprintf('trial %i: image',k);
        tstr{end + 1} = sprintf('trial %i: abs',  k);
    end
    for i = 1:nm
        for k = 1:nt
            mstr{end + 1} = sprintf('mode %i trial %i',i,k);
        end
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j),cla
            for k = 1:nt
                plot(Hz,real(xY.y{k}(:,i,j)),'-.','color',co{k}), hold on
                plot(Hz,imag(xY.y{k}(:,i,j)),':', 'color',co{k}), hold on
                plot(Hz, abs(xY.y{k}(:,i,j)),     'color',co{k}), hold on
                title(sprintf('%s to %s',name{j},name{i}))
                axis square, spm_axis tight
            end
        end
 
        % spectral density
        %------------------------------------------------------------------
        subplot(2,2,3)
        for k = 1:nt
            plot(Hz,abs(xY.y{k}(:,i,i)),'color',co{i},'LineStyle',lt{k}), hold on
            axis square, spm_axis tight
        end
    end
    
    title('Spectral density over modes')
    xlabel('Frequency (Hz)')
    ylabel('CSD')
    axis square, spm_axis tight
    legend(mstr)
    
    subplot(nm,nm,nm)
    legend(tstr)
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
    
    % spm_dcm_csd_results(DCM,'coupling (A)');
    %----------------------------------------------------------------------
    str = {'Forward (i)','Forward (ii)','Backward (i)','Backward (ii)'};
    m   = length(DCM.Ep.A);
    for  i = 1:m
        
        % images
        %------------------------------------------------------------------
        subplot(4,m,i)
        imagesc(exp(DCM.Ep.A{i}))
        title(str{i},'FontSize',10)
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square
    
        % table
        %------------------------------------------------------------------
        subplot(4,m,i + m)
        text(0,1/2,num2str(full(exp(DCM.Ep.A{i})),' %.2f'),'FontSize',8)
        axis off,axis square
 
    
        % PPM
        %------------------------------------------------------------------
        subplot(4,m,i + m + m)
        image(64*DCM.Pp.A{i})
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        title('PPM')
        axis square
    
        % table
        %------------------------------------------------------------------
        subplot(4,m,i + m + m + m)
        text(0,1/2,num2str(DCM.Pp.A{i},' %.2f'),'FontSize',8)
        axis off, axis square
        
    end
    
case{lower('Coupling (C)')}
    
    % spm_dcm_csd_results(DCM,'coupling (C)');
    %----------------------------------------------------------------------
    
    % images
    %----------------------------------------------------------------------
    subplot(2,4,1)
    imagesc(exp(DCM.Ep.C))
    title('Factors','FontSize',10)
    set(gca,'XTick',1:nu,'XTickLabel','Input','FontSize',8)
    set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname, 'FontSize',8)
    axis square
    
    % PPM
    %----------------------------------------------------------------------
    subplot(2,4,3)
    image(64*DCM.Pp.C)
    title('Factors','FontSize',10)
    set(gca,'XTick',1:nu,'XTickLabel','Input','FontSize',8)
    set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname, 'FontSize',8)
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
    
    % spm_dcm_csd_results(DCM,'coupling (B)');
    %----------------------------------------------------------------------
    for i = 1:nu
        
        % images
        %------------------------------------------------------------------
        subplot(4,nu,i)
        imagesc(exp(DCM.Ep.B{i}))
        title(DCM.xU.name{i},'FontSize',10)
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
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
        set(gca,'YTick',1:ns,'YTickLabel',DCM.Sname,'FontSize',8)
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
    
    % spm_dcm_csd_results(DCM,'trial-specific effects');
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
    [Gu,Gs,Gn,Hz] = spm_csd_mtf_gu(DCM.Ep,DCM.M);
 
    
    % plot spectral density of innovations
    % ---------------------------------------------------------------------
    subplot(2,1,1)
    plot(Hz,Gu)
    xlabel('frequency (Hz)')
    title('Spectrum of innovations','FontSize',16)
    axis square, grid on, spm_axis scale
    legend(DCM.Sname)
 
    % plot spectral density of noise
    % ---------------------------------------------------------------------
    subplot(2,2,3)
    plot(Hz,Gs)
    xlabel('frequency (Hz)')
    title('Channel-specific noise')
    axis square, grid on, spm_axis scale
    legend(xY.name)
 
    % plot spectral density of noise
    % ---------------------------------------------------------------------
    subplot(2,2,4)
    plot(Hz,Gn)
    xlabel('frequency (Hz)')
    title('Non-specific noise')
    axis square, grid on, spm_axis scale
    
    
    
case{lower('Transfer functions')}
    
    % spm_dcm_csd_results(DCM,'Cross-spectral density');
    %----------------------------------------------------------------------
    co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz   = DCM.Hz(:);
    name = DCM.Sname;
    nm   = length(name);
    
    
    % spectrum of innovations or noise (Gu)
    %----------------------------------------------------------------------
    Gu   = spm_csd_mtf_gu(DCM.Ep,Hz);
    
    
    tstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('trial %i: abs',k);
        tstr{end + 1} = sprintf('trial %i: phase (sec) ',k);
    end
    
    for i = 1:nm
        for j = 1:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j)
            for k = 1:nt
                
                dtf = DCM.dtf{k}(:,i,j);
                fsd = unwrap(angle(dtf))./(2*pi*Hz);
                dtf = abs(dtf.*sqrt(Gu(:,j)));
                
                plot(Hz,dtf,    'color',co{k}), hold on
                plot(Hz,fsd,':','color',co{k}), hold on
                title(sprintf('Spectral transfer: %s to %s',name{j},name{i}))
                xlabel('frequency Hz')
                axis square, spm_axis tight
            end
        end
 
        % legend
        %------------------------------------------------------------------      
        if i == nm && j == nm
            legend(tstr)
        end
    end
   
    
case{lower('Cross-spectra (sources)')}
    
    % spm_dcm_csd_results(DCM,'Cross-spectral density');
    %----------------------------------------------------------------------
    co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz   = DCM.Hz;
    name = DCM.Sname;
    nm   = length(name);
    
    tstr = {};
    mstr = {};
    for k = 1:nt
        mstr{end + 1} = sprintf('trial: %i',k);
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j)
            for k = 1:nt
                plot(Hz,abs(DCM.Hs{k}(:,i,j)),'color',co{k}), hold on
                title(sprintf('CSD: %s to %s',name{j},name{i}))
                xlabel('frequency Hz')
                axis square, spm_axis tight
            end
        end
 
        % legend
        %------------------------------------------------------------------      
        if i == nm && j == nm
            legend(mstr)
        end
        
        % spectral density
        %------------------------------------------------------------------
        subplot(2,2,3)
        for k = 1:nt
            plot(Hz,abs(DCM.Hs{k}(:,i,i)),'color',co{i}), hold on
            axis square, spm_axis tight
            tstr{end + 1} = sprintf('%s: trial %i',name{i},k);
        end
    end
   
    title({'Spectral density (source-space)'},'FontSize',16)
    xlabel('frequency (Hz)')
    ylabel('abs(CSD)')
    axis square, spm_axis tight
    legend(tstr)
    
    
case{lower('Cross-spectra (channels)')}
    
    % spm_dcm_csd_results(DCM,'Cross-spectral density');
    %----------------------------------------------------------------------
    co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz   = DCM.Hz;
    nm   = min(nm,4);
    
    tstr = {};
    mstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('predicted: trial %i',k);
        tstr{end + 1} = sprintf('observed: trial %i',k);
    end
    for i = 1:nm
        for k = 1:nt
            mstr{end + 1} = sprintf('predicted: mode %i trial %i',i,k);
            mstr{end + 1} = sprintf( 'observed: mode %i trial %i',i,k);
        end
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j)
            for k = 1:nt
                plot(Hz,abs(DCM.Hc{k}(:,i,j)),'color',co{k}), hold on
                plot(Hz,abs(DCM.Hc{k}(:,i,j) + DCM.Rc{k}(:,i,j)),':','color',co{k})
                title(sprintf('mode %i to %i',j,i))
                xlabel('frequency Hz')
                axis square, spm_axis tight
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
            plot(Hz,abs(DCM.Hc{k}(:,i,i)),'color',co{i}), hold on
            plot(Hz,abs(DCM.Hc{k}(:,i,i) + DCM.Rc{k}(:,i,i)),':','color',co{i})
            axis square, spm_axis tight
        end
    end
   
    title({'Spectral density over modes';'(in channel-space)'},'FontSize',16)
    xlabel('frequency (Hz)')
    ylabel('abs(CSD)')
    legend(mstr)
    
    
case{lower('Coherence (sources)')}
    
    % get coherence
    %----------------------------------------------------------------------
    [coh fsd] = spm_csd2coh(DCM.Hs,DCM.Hz);
    
    % spm_dcm_csd_results(DCM, 'Coherence (sources)');
    %----------------------------------------------------------------------
    co    = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    Hz    = DCM.Hz;
    name  = DCM.Sname;
    nm    = length(name);
    tstr  = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('trial %i',k);
    end
    
    for i = 1:nm
        for j = 1:nm
 
            if j > i
                
                % for each trial type - coherence
                %----------------------------------------------------------
                subplot(nm,nm,(i - 1)*nm + j),cla
                for k = 1:nt
                    plot(Hz,coh{k}(:,i,j),'color',co{k}), hold on
                    title(sprintf('Coh: %s to %s',name{j},name{i}))
                    xlabel('frequency Hz')
                    axis square, spm_axis tight
                end
            
            end
            
            if j < i
                % for each trial type - delay
                %--------------------------------------------------------------
                subplot(nm,nm,(i - 1)*nm + j),cla
                for k = 1:nt
                    plot(Hz,1000*fsd{k}(:,i,j),'color',co{k}), hold on
                    title(sprintf('Delay (ms) %s to %s',name{j},name{i}))
                    xlabel('Frequency Hz')
                    axis square, spm_axis tight
                end  
            end
           
        end
    end
    legend(tstr)
    
    
case{lower('Coherence (channels)')}
     
    % get coherence
    %----------------------------------------------------------------------
    [coh fsd] = spm_csd2coh(DCM.Hc,DCM.Hz);
    [COH FSD] = spm_csd2coh(DCM.xY.y,DCM.Hz);
 
    
    % spm_dcm_csd_results(DCM,'Coherence (channels)');
    %----------------------------------------------------------------------
    co    = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    nm    = min(nm,4);
    tstr  = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('predicted: trial %i',k);
        tstr{end + 1} = sprintf('observed: trial %i',k);
    end
 
    
    for i = 1:nm
        for j = 1:nm
 
            if j > i
                
                % for each trial type - coherence
                %----------------------------------------------------------
                subplot(nm,nm,(i - 1)*nm + j)
                for k = 1:nt
                    plot(Hz,coh{k}(:,i,j),'color',co{k}), hold on
                    plot(Hz,COH{k}(:,i,j),':','color',co{k}), hold on
                    title(sprintf('Coh: %i to %i',j,i))
                    xlabel('frequency Hz')
                    axis square, spm_axis tight
                end
            
            end
            
            if j < i
                % for each trial type - delay
                %----------------------------------------------------------
                subplot(nm,nm,(i - 1)*nm + j)
                for k = 1:nt
                    plot(Hz,1000*fsd{k}(:,i,j),'color',co{k}), hold on
                    plot(Hz,1000*FSD{k}(:,i,j),':','color',co{k}), hold on
                    title(sprintf('Delay (ms) %i to %i',j,i))
                    xlabel('frequency Hz')
                    axis square, spm_axis tight
                end  
            end
           
        end
    end
    legend(tstr)
    
    
case{lower('Covariance (sources)')}
    
    % get covariances
    %----------------------------------------------------------------------
    [ccf pst] = spm_csd2ccf(DCM.Hs,DCM.Hz);
    
    
    % spm_dcm_csd_results(DCM,'Cross-covariance (sources)');
    %----------------------------------------------------------------------
    co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    pst  = 1000*pst;
    name = DCM.Sname;
    nm   = length(name);
    q    = max(spm_vec(ccf));
    p    = min(spm_vec(ccf));
    tstr = {};
    mstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('trial %i',k);
    end
    for k = 1:nm
        mstr{end + 1} = sprintf('source %i',k);
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j),cla
            plot([0 0],[p q],':'), hold on
            for k = 1:nt
                plot(pst,ccf{k}(:,i,j),'color',co{k}), hold on
                title(sprintf('%s to %s',name{j},name{i}))
                xlabel('lag (ms)')
                axis tight square, set(gca,'XLim',[-128 128],'YLim',[p q])
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
            plot(pst,ccf{k}(:,i,i),'color',co{i}), hold on
            axis tight square, set(gca,'XLim',[-128 128],'YLim',[p q])
        end
    end
   
    title({'Auto-covariance';'(in source-space)'},'FontSize',16)
    xlabel('lag (ms)')
    ylabel('auto-covariance')
    axis square
    legend(mstr)
    
    
case{lower('Covariance (channels)')}
    
    % get covariances
    %----------------------------------------------------------------------
    [ccf pst] = spm_csd2ccf(DCM.Hc,DCM.Hz);
    
    
    % spm_dcm_csd_results(DCM,'Cross-covariance (sources)');
    %----------------------------------------------------------------------
    co   = {'b', 'r', 'g', 'm', 'y', 'k', 'c'};
    pst  = 1000*pst;
    nm   = min(nm,4);
    q    = max(spm_vec(ccf));
    p    = min(spm_vec(ccf));
    tstr = {};
    mstr = {};
    for k = 1:nt
        tstr{end + 1} = sprintf('trial %i',k);
    end
    for k = 1:nm
        mstr{end + 1} = sprintf('channel %i',k);
    end
    
    for i = 1:nm
        for j = i:nm
 
            % for each trial type
            %--------------------------------------------------------------
            subplot(nm,nm,(i - 1)*nm + j),cla
            plot([0 0],[p q],':'), hold on
            for k = 1:nt
                plot(pst,ccf{k}(:,i,j),'color',co{k}), hold on
                title(sprintf('mode %i to %i',j,i))
                xlabel('lag (ms)')
                axis tight square, set(gca,'XLim',[-128 128],'YLim',[p q])
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
            plot(pst,ccf{k}(:,i,i),'color',co{i}), hold on
            axis tight square, set(gca,'XLim',[-128 128],'YLim',[p q])
        end
    end
   
    title({'Auto-covariance';'(in channel-space)'},'FontSize',16)
    xlabel('Lag (ms)')
    ylabel('auto-covariance')
    axis square
    legend(mstr)
    
    
case{lower('Dipoles')}
    
    % return if LFP
    % ---------------------------------------------------------------------
    if strcmpi(DCM.xY.modality,'lfp')
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
