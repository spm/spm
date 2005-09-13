function [DCM] = spm_dcm_erp_results(DCM,Action)
% Results for ERP Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_erp_results(DCM,'ERPs (channel)');
% FORMAT spm_dcm_erp_results(DCM,'ERPs (sources)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (A)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (B)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (C)');
% FORMAT spm_dcm_erp_results(DCM,'Input');
% FORMAT spm_dcm_erp_results(DCM,'Response');
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
% 
%_______________________________________________________________________
% %W% Karl Friston %E%


% get figure handle
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
figure(Fgraph)
clf
try
    nt     = length(DCM.H);
    nu     = size(DCM.C,2);
    nc     = size(DCM.H{1},2);
    ns     = size(DCM.K{1},2);
end

% switch
%------------------------------------------------------------------------
switch(lower(Action))
    
case{lower('ERPs (channel)')}
    
    % spm_dcm_erp_results(DCM,'ERPs (channel)');
    %--------------------------------------------------------------------
    co = {'k', 'k'};
    lo = {'-', '--'};
    T     = [1:size(DCM.H{1},1)]*DCM.xY.dt*1000;
    for i = 1:nc
        subplot(ceil(nc/2),2,i), hold on
        str   = {};
        for j = 1:nt
            plot(T,DCM.H{j}(:,i), lo{1},...
                'Color', co{j},...
                'LineWidth',2);
            str{end + 1} = sprintf('trial %i (predicted)',j);
            plot(T,DCM.H{j}(:,i) + DCM.R{j}(:,i), lo{2},...
                'Color',co{j});
            str{end + 1} = sprintf('trial %i (observed)',j);
                        set(gca, 'XLim', [0 T(end)]);

        end
        hold off
        title(sprintf('channel %i',i))
        grid on
        axis square
    end
    xlabel('time (ms)')
    legend(str)
    
    
case{lower('ERPs (sources)')}
    
    % spm_dcm_erp_results(DCM,'ERPs (sources)');
    %--------------------------------------------------------------------
    T     = [1:size(DCM.H{1},1)]*DCM.xY.dt*1000;
    
    mx = max(max(cat(2, DCM.K{:})));
    mi = min(min(cat(2, DCM.K{:})));
    
    for i = 1:ns
        subplot(ceil(ns/2),2,i), hold on
        str   = {};
        for j = 1:nt
            plot(T,DCM.K{j}(:,i), ...
                'Color',[1 1 1] - j/nt, ...
                'LineWidth',1);
            str{end + 1} = sprintf('trial %i',j);
        end
        set(gca, 'YLim', [mi mx]);
        hold off
        title(DCM.Sname{i})
        grid on
        axis square
    end
    xlabel('time (ms)')
    legend(str)
    
    
case{lower('Coupling (A)')}
    
    % spm_dcm_erp_results(DCM,'coupling (A)');
    %--------------------------------------------------------------------
    str = {'Forward','Backward','Lateral'};
    for  i =1:3
        
        % images
        %-----------------------------------------------------------
        subplot(4,3,i)
        imagesc(exp(DCM.Qp.A{i}))
        title(str{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square
    
        % tables
        %--------------------------------------------------------------------
        subplot(4,3,i + 3)
        text(0,1/2,num2str(full(exp(DCM.Qp.A{i})),' %.2f'),'FontSize',8)
        axis off,axis square

    
        % images
        %-----------------------------------------------------------
        subplot(4,3,i + 6)
        imagesc(exp(DCM.Pp.A{i}))
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        axis square
    
        % tables
        %--------------------------------------------------------------------
        subplot(4,3,i + 9)
        text(0,1/2,num2str(DCM.Pp.A{i},' %.2f'),'FontSize',8)
        axis off, axis square
        
    end
    
case{lower('Coupling (C)')}
    
    % spm_dcm_erp_results(DCM,'coupling (A)');
    %--------------------------------------------------------------------
    
    % images
    %-----------------------------------------------------------
    subplot(2,4,1)
    imagesc(exp(DCM.Qp.C))
    title('Factors','FontSize',10)
    set(gca,'XTick',[1:nu],'XTickLabel',DCM.U.name,'FontSize',8)
    set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname, 'FontSize',8)
    axis square
    
    subplot(2,4,3)
    imagesc(DCM.Pp.C)
    title('Factors','FontSize',10)
    set(gca,'XTick',[1:nu],'XTickLabel',DCM.U.name,'FontSize',8)
    set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname, 'FontSize',8)
    axis square
    
    % table
    %--------------------------------------------------------------------
    subplot(2,4,2)
    text(0,1/2,num2str(full(exp(DCM.Qp.C)),' %.2f'),'FontSize',8)
    axis off

    % table
    %--------------------------------------------------------------------
    subplot(2,4,4)
    text(0,1/2,num2str(DCM.Pp.C,     ' %.2f'),'FontSize',8)
    axis off

 
case{lower('Coupling (B)')}
    
    % spm_dcm_erp_results(DCM,'coupling (B)');
    %--------------------------------------------------------------------
    for i = 1:nu
        
        % images
        %-----------------------------------------------------------
        subplot(4,nu,i)
        imagesc(exp(DCM.Qp.B{i}))
        title(DCM.U.name{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square

        % tables
        %--------------------------------------------------------------------
        subplot(4,nu,i + nu)
        text(0,1/2,num2str(full(exp(DCM.Qp.B{i})),' %.2f'),'FontSize',8)
        axis off
        axis square
        
        % images
        %-----------------------------------------------------------
        subplot(4,nu,i + 2*nu)
        imagesc(DCM.Pp.B{i})
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        axis square

        % tables
        %--------------------------------------------------------------------
        subplot(4,nu,i + 3*nu)
        text(0,1/2,num2str(DCM.Pp.B{i},' %.2f'),'FontSize',8)
        axis off
        axis square
        
    end
    
case{lower('Input')}
    
    % plot data
    % --------------------------------------------------------------------
    xU   = DCM.xU;
    t    = xU.dt:xU.dt:xU.dur;
    U    = spm_erp_u(DCM.Qp,t);
    
    subplot(2,1,1)
    plot(t*1000,U)
    xlabel('time (ms)')
    title('input')
    axis square, grid on

    
case{lower('Response')}
    
    % plot data
    % --------------------------------------------------------------------
    xY    = DCM.Y;
    n     = length(xY.xy);
    try
        dt = xY.dt;
    catch
        dt = 1;
    end
    for i = 1:n
        subplot(n,1,i)
        plot([1:size(xY.xy{i},1)]*dt,xY.xy{i})
        xlabel('time (ms)')
        title(sprintf('ERP %i',i))
        axis square, grid on
    end

end
