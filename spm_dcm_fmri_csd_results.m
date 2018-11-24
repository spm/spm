function spm_dcm_fmri_csd_results(DCM,action,fig)
% Review an estimated DCM for BOLD CSD
% FORMAT spm_dcm_fmri_csd_results(DCM,action,fig)
%
% Action:
%     'Spectral data'
%     'Coupling (A)'
%     'Coupling (C)'
%     'Inputs'
%     'Outputs'
%     'Transfer functions'
%     'Cross-spectra (BOLD)'
%     'Cross-spectra (neural)'
%     'Coherence (neural)'
%     'Covariance (neural)'
%     'Kernels'
%     'Functional connectivity'
%     'Location of regions'
%     'Quit'
%__________________________________________________________________________
% Copyright (C) 2013-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_csd_results.m 7497 2018-11-24 17:00:25Z karl $


%-Input arguments
%==========================================================================
 
%-Get DCM structure
%--------------------------------------------------------------------------
if nargin < 1
    [DCM, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
end
if ~isstruct(DCM)
    load(DCM);
end
 
%-Get action if necessary
%--------------------------------------------------------------------------
if nargin < 2 || isempty(action)
    str = {'Timeseries data',...
        'Spectral data',...
        'Coupling (A)',...
        'Coupling (C)',...
        'Inputs',...
        'Transfer functions',...
        'Cross-spectra (BOLD)',...
        'Cross-spectra (neural)',...
        'Coherence (neural)',...
        'Covariance (neural)',...
        'Kernels',...
        'Functional connectivity',...
        'Location of regions',...
        'Quit'};
    
    action = str{spm_input('review',1,'m',str)};
end
 
%-Get figure
%--------------------------------------------------------------------------
if nargin < 3 || isempty(fig)
    Fgraph = spm_figure('GetWin','Graphics');
else
    Fgraph = fig;
    spm_figure('Focus',Fgraph);
end
colormap(gray)
if ~strcmpi(action,'Quit')
    spm_figure('Clear',Fgraph);
end

 
%-CSD data
%--------------------------------------------------------------------------
xY   = DCM.Y;                  % data
Hz   = DCM.Hz;                 % frequencies
name = {DCM.xY.name};          % names
ns   = size(DCM.a,1);          % number of regions
nu   = size(DCM.c,2);          % number of exogenous inputs
U    = 0.9;                    % p-value threshold for display
ns   = min(ns,8);              % bounded number of regions

 
%-Exponentiate coupling parameters if a two-state model
%--------------------------------------------------------------------------
try
    if DCM.options.two_state
        Ep.A = full(exp(DCM.Ep.A));
        Ep.B = full(exp(DCM.Ep.B));
        Ep.D = full(exp(DCM.Ep.D));
        Ep.C = full(DCM.Ep.C);
        disp('NB: The (A,B,D) parameters are scale parameters')
    else
        Ep.A = full(DCM.Ep.A);
        Ep.B = full(DCM.Ep.B);
        Ep.D = full(DCM.Ep.D);
        Ep.C = full(DCM.Ep.C);
    end
end
 
 
% switch
%--------------------------------------------------------------------------
switch lower(action)
            
    %======================================================================
    % Timeseries data
    %======================================================================
    case lower('Timeseries data')
        
        % graph
        %------------------------------------------------------------------
        x     = (1:DCM.v)*DCM.Y.dt;
        for i = 1:ns
            subplot(ns,1,i);
            plot(x,DCM.Y.y(:,i));
            title([DCM.Y.name{i} ': responses'],'FontSize',16);
            xlabel('time {seconds}');
        end
        
    %======================================================================
    % spectral data
    %======================================================================
    case lower('Spectral data')
        
        for i = 1:ns
            for j = i:ns
                subplot(ns,ns,(i - 1)*ns + j),cla
                plot(Hz,real(xY.csd(:,i,j))), hold on
                plot(Hz,imag(xY.csd(:,i,j)),':'), hold on
                title(sprintf('%s to %s',name{j},name{i}))
                axis square, spm_axis tight
            end
            
            % spectral density
            %--------------------------------------------------------------
            subplot(2,2,3)
            csd(:,i) = abs(xY.csd(:,i,i));
        end
        
        % spectral density
        %------------------------------------------------------------------
        subplot(2,2,3)
        plot(Hz,csd)
        title('Spectral density over modes')
        xlabel('Frequency (Hz)')
        ylabel('CSD')
        axis square, spm_axis tight
        legend(name),legend('boxoff')
        
        subplot(ns,ns,ns)
        legend('real','imag'),legend('boxoff')
    
    %======================================================================
    % Location of regions
    %======================================================================
    case lower('Location of regions')
        
        % graph
        %------------------------------------------------------------------
        subplot(2,1,1)
        spm_dcm_graph(DCM.xY);
        title('Regional locations','FontSize',16)
        
        % table
        %------------------------------------------------------------------
        subplot(2,1,2)
        y = 0;
        line([0 4],[y y])
        y = y - 1;
        text(0.0,y,'Name',         'FontSize',14)
        text(1.2,y,'Voxels',       'FontSize',14)
        text(2.0,y,'Location (mm)','FontSize',14)
        y = y - 1;
        line([0 4],[y y],'LineWidth',4)
        y = y - 1;
        for i = 1:length(DCM.xY)
            N    = size(DCM.xY(i).XYZmm, 2);
            L    = DCM.xY(i).xyz;
            r    = DCM.xY(i).spec;
            text(0.0,y,name{i}, 'FontWeight','bold',        'FontSize',12)
            text(1.5,y,sprintf('%0.0f',N),               'FontSize',12)
            text(2.0,y,sprintf('%-4.0f %-4.0f %-4.0f',L),'FontSize',12)
            y = y - 1;
        end
        line([0 4],[y y])
        axis off square
        
        
    %======================================================================
    % Inputs
    %======================================================================
    case lower('Inputs')
        
        % inputs
        %------------------------------------------------------------------
        try
            t     = (1:length(DCM.U.u))*DCM.U.dt;
            for i = 1:nu
                
                subplot(nu,1,i)
                plot(t,full(DCM.U.u(:,i)))
                title(['Input - ' DCM.U.name{i}],'FontSize',16)
                ylabel('event density {Hz}')
                spm_axis tight
                
            end
            xlabel('time {seconds}')
        end
        
    %======================================================================
    % Parameter estimates (A)
    %======================================================================
    case lower('Coupling (A)')
        
        % intrinsic effects
        %------------------------------------------------------------------
        subplot(2,1,1)
        a(:,:,1) = DCM.Pp.A;
        a(:,:,2) = Ep.A;
        spm_dcm_display(DCM.xY,a)
        title(sprintf('%s P(coupling > 0)','fixed'),'FontSize',16)
        
        
        % intrinsic interactions
        %------------------------------------------------------------------
        subplot(2,2,3)
        bar(a(:,:,2))
        title('A - fixed effects','FontSize',16)
        set(gca,'XTick',(1:ns),'XTickLabel',name)
        xlabel('target region')
        ylabel('strength (Hz)')
        legend(name),legend('boxoff')
        
        
        % intrinsic interactions - probabilities
        %------------------------------------------------------------------
        subplot(2,2,4)
        P = a(:,:,1);
        bar(P), hold on, plot([0 (length(P) + 1)],[U U],'-.r'), hold off
        title('A - probability','FontSize',16)
        set(gca,'XTick',1:ns,'XTickLabel',name)
        xlabel('target region')
        ylabel(sprintf('P(A > %0.2f)',0))
        
        
    %======================================================================
    % Parameter estimates (C)
    %======================================================================
    case lower('Coupling (C)')
        
        if spm_length(Ep.C)
            
            % intrinsic effects
            %--------------------------------------------------------------
            subplot(2,1,1)
            c(:,:,1) = DCM.Pp.C;
            c(:,:,2) = Ep.C;
            spm_dcm_display(DCM.xY,[],c)
            title(sprintf('%s P(coupling > 0)','fixed'),'FontSize',16)
            
            
            % intrinsic interactions
            %--------------------------------------------------------------
            subplot(2,2,3)
            bar(c(:,:,2))
            title('C - exogenous effects','FontSize',16)
            set(gca,'XTick',(1:ns),'XTickLabel',name)
            xlabel('target region')
            ylabel('strength (Hz)')
            legend(DCM.U.name),legend('boxoff')
            
            
            % intrinsic interactions - probabilities
            %--------------------------------------------------------------
            subplot(2,2,4)
            bar(c(:,:,1))
            title('C - probability','FontSize',16)
            set(gca,'XTick',(1:ns),'XTickLabel',name)
            xlabel('target region')
            ylabel(sprintf('P(C > %0.2f)',0))
            
        end
        
        
    %======================================================================
    % Transfer functions and cross spectra
    %======================================================================
    case lower('Transfer functions')
        
        for i = 1:ns
            for j = 1:ns
                
                subplot(ns,ns,(i - 1)*ns + j)
                dtf = abs(DCM.dtf(:,i,j));
                plot(Hz,dtf), hold on
                title(sprintf('Spectral transfer: %s to %s',name{j},name{i}))
                xlabel('frequency Hz')
                axis square, spm_axis tight
                
            end
        end
        
        
    %======================================================================
    % Transfer functions and cross spectra
    %======================================================================
    case lower('Cross-spectra (BOLD)')
        
        Gc    = DCM.Hc;                     % predicted
        Rc    = Gc + DCM.Rc;                % observed
        for i = 1:ns
            for j = i:ns
                subplot(ns,ns,(i - 1)*ns + j)
                plot(Hz,real(Gc(:,i,j)),'r'),  hold on
                plot(Hz,imag(Gc(:,i,j)),'c')
                plot(Hz,real(Rc(:,i,j)),'r:') 
                plot(Hz,imag(Rc(:,i,j)),'c:'), hold off
                title(sprintf('CSD: %s to %s',name{j},name{i}))
                xlabel('frequency Hz')
                axis square, spm_axis tight
            end
            
            % spectral density
            %--------------------------------------------------------------
            Hc(:,i) = abs(Gc(:,i,i));
            Yc(:,i) = abs(Rc(:,i,i));
            
        end
        
        subplot(ns,ns,ns)
        legend({'pred - real','pred - imag','obs - real','obs - imag'})
        legend('boxoff')
        
        subplot(2,2,3)
        plot(Hz,Hc), hold on
        plot(Hz,Yc,':'), hold off
        title({'Spectral density (BOLD)'},'FontSize',16)
        xlabel('frequency (Hz)')
        ylabel('abs(CSD)')
        axis square, spm_axis tight
        legend(name),legend('boxoff')
        
        
        
    %======================================================================
    % Cross-spectra (neural)
    %====================================================================== 
    case lower('Cross-spectra (neural)')
        
        for i = 1:ns
            for j = i:ns
                subplot(ns,ns,(i - 1)*ns + j) 
                plot(Hz,real(DCM.Hs(:,i,j)),'r'), hold on
                plot(Hz,imag(DCM.Hs(:,i,j)),'c'), hold off
                title(sprintf('CSD: %s to %s',name{j},name{i}))
                xlabel('frequency Hz')
                axis square, spm_axis tight
                
            end
            
            % spectral density
            %--------------------------------------------------------------
            Hs(:,i) = abs(DCM.Hs(:,i,i));
            
        end
        subplot(ns,ns,ns)
        legend({'real','imag'})
        legend('boxoff')
        
        subplot(2,2,3)
        plot(Hz,Hs), hold on
        title({'Spectral density (neural)'},'FontSize',16)
        xlabel('frequency (Hz)')
        ylabel('abs(CSD)')
        axis square, spm_axis tight
        legend(name),legend('boxoff')
        
        
    %======================================================================
    % Coherence (neural)
    %======================================================================  
    case lower('Coherence (neural)')
        
        
        for i = 1:ns
            for j = 1:ns
                
                % coherence
                %---------------------------------------------------------
                if j > i
                    subplot(ns,ns,(i - 1)*ns + j),cla
                    plot(Hz,DCM.coh(:,i,j))
                    title(sprintf('Coh: %s to %s',name{j},name{i}))
                    xlabel('frequency Hz')
                    axis square, spm_axis tight
                end
                
                % delays
                %----------------------------------------------------------
                if j < i
                    subplot(ns,ns,(i - 1)*ns + j),cla  
                    plot(Hz,DCM.fsd(:,i,j))
                    title(sprintf('Delay (s) %s to %s',name{j},name{i}))
                    xlabel('Frequency Hz')
                    axis square, spm_axis tight         
                end
                
            end
        end
        
        
    %======================================================================
    % Cross covariance functions
    %======================================================================
    case lower('Covariance (neural)')
        
        for i = 1:ns
            for j = i:ns
            
                subplot(ns,ns,(i - 1)*ns + j),cla 
                plot(DCM.pst,DCM.ccf(:,i,j))
                title(sprintf('%s to %s',name{j},name{i}))
                xlabel('lag (s)')
                axis tight square
                
            end
            
            % cross covariance functions
            %--------------------------------------------------------------
            ccf(:,i) = DCM.ccf(:,i,i);
            
        end
        
        % cross covariance functions
        %--------------------------------------------------------------
        subplot(2,2,3)
        plot(DCM.pst,ccf)
        title({'Auto-covariance';'(neural)'},'FontSize',16)
        xlabel('lag (s)')
        ylabel('auto-covariance')
        axis square
        
    %======================================================================
    % Functional connectivity
    %======================================================================
    case lower('Functional connectivity')
        
        % predicted responses
        %------------------------------------------------------------------
        [Gy,Hz,S,Gu,Gn] = spm_csd_fmri_mtf(DCM.Ep,DCM.M,DCM.U);
        
        % spectral density
        %------------------------------------------------------------------
        for i = 1:ns   
            Yc(:,i) = abs(DCM.Y.csd(:,i,i));          % observed
            Hc(:,i) = abs(Gy(:,i,i));                 % predicted
            Nn(:,i) = abs(Gn(:,i,i));                 % Non-neuronal
            Gc(:,i) = abs(Gy(:,i,i) - Gn(:,i,i));     % Neuronal
        end
        
        
        % zero-lag correlations (functional connectivity)
        %------------------------------------------------------------------
        subplot(2,4,1)
        ccf = spm_csd2ccf(DCM.Y.csd,Hz,DCM.Y.dt);
        ccf = ccf((size(ccf,1) + 1)/2,:,:);
        cc  = spm_cov2corr(squeeze(ccf));
        image(cc*64), title('Correlations','FontSize',16), axis square
        
        subplot(2,4,2)
        ccf = spm_csd2ccf(Gy,Hz,DCM.Y.dt);
        ccf = ccf((size(ccf,1) + 1)/2,:,:);
        cc  = spm_cov2corr(squeeze(ccf));
        image(cc*64), title('Predicted','FontSize',16), axis square
        
        subplot(2,4,3)
        ccf = spm_csd2ccf(Gn,Hz,DCM.Y.dt);
        ccf = ccf((size(ccf,1) + 1)/2,:,:);
        cc  = spm_cov2corr(squeeze(ccf));
        image(cc*64), title('Non-neuronal','FontSize',16), axis square
        
        subplot(2,4,4)
        ccf = spm_csd2ccf(Gy - Gn,Hz,DCM.Y.dt);
        ccf = ccf((size(ccf,1) + 1)/2,:,:);
        cc  = spm_cov2corr(squeeze(ccf));
        image(cc*64), title('Neuronal','FontSize',16), axis square
        
        % spectral decomposition
        %------------------------------------------------------------------
        subplot(2,2,3)
        plot(Hz,Hc),    hold on
        plot(Hz,Yc,':'),hold off
        title('Predicted and observed spectra','FontSize',16)
        axis square, spm_axis tight
        legend(name),legend('boxoff')
        
        subplot(2,2,4)
        plot(Hz,Hc,':'), hold on
        plot(Hz,Gc),     hold on
        plot(Hz,Nn,'-.'), hold off
        title('Neuronal and non-neuronal part','FontSize',16)
        axis square, spm_axis tight
        legend(name),legend('boxoff')

        
    %======================================================================
    % Kernels
    %======================================================================
    case lower('Kernels')
        
        % input effects
        %------------------------------------------------------------------
        x     = (1:DCM.M.N)*DCM.M.dt;
        for i = 1:ns
            
            % input effects - neural
            %--------------------------------------------------------------
            y = DCM.K1(:,:,i);
            subplot(ns,2,2*(i - 1) + 1)
            plot(x,y)
            set(gca,'XLim',[0 16])
            axis square
            title(['neural responses to ' name{i}])
            xlabel('time {seconds}')
            
            
            % input effects - BOLD
            %--------------------------------------------------------------
            y = DCM.K1(:,:,i);
            k = DCM.H1(:,:,i);
            subplot(ns,2,2*(i - 1) + 2)
            plot(x,y,':'), hold on
            plot(x,k)
            set(gca,'XLim',[0 16])
            axis square
            title('BOLD responses','FontSize',12)
            
        end
        legend(name),legend('boxoff')
        
        
        
    %======================================================================
    % Quit
    %======================================================================
    case lower('Quit')
        return;
        
    %======================================================================
    % Unknown action
    %======================================================================
    otherwise
        disp('Unknown option.')
end
 
% return to menu
%--------------------------------------------------------------------------
if nargin < 2
    spm_dcm_fmri_csd_results(DCM);
end
