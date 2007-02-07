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
%__________________________________________________________________________
% %W% Karl Friston %E


% get figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

if ~strcmp(lower(Action), 'dipoles')
    figure(Fgraph)
    clf
end

try
    nt     = length(DCM.H);          % Nr of trials
    nu     = length(DCM.B);          % Nr inputs
    nc     = size(DCM.H{1},2);       % Nr modes
    ns     = size(DCM.K{1},2);       % Nr of sources
end

% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
case{lower('ERPs (channel)')}
    
    xY    = DCM.xY;
    n     = length(xY.xy);
    
    try
        t = xY.Time;
    catch
        t = xY.dt*[1:n];
    end


    % spm_dcm_erp_results(DCM,'ERPs (channel)');
    %----------------------------------------------------------------------
    co = {'b', 'r', 'g', 'm', 'y', 'k'};
    lo = {'-', '--'};
    
    for i = 1:nc
        subplot(ceil(nc/2),2,i), hold on
        str   = {};
        for j = 1:nt
            plot(t,DCM.H{j}(:,i), lo{1},...
                'Color', co{j},...
                'LineWidth',2);
            str{end + 1} = sprintf('trial %i (predicted)',j);
            plot(t,DCM.H{j}(:,i) + DCM.R{j}(:,i), lo{2},...
                'Color',co{j});
            str{end + 1} = sprintf('trial %i (observed)',j);
                        set(gca, 'XLim', [t(1) t(end)]);

        end
        hold off
        title(sprintf('channel %i',i))
        grid on
        axis square
        try
            axis(A);
        catch
            A = axis;
        end
    end
    xlabel('time (ms)')
    legend(str)
    
    
case{lower('ERPs (sources)')}
    
    % spm_dcm_erp_results(DCM,'ERPs (sources)');
    %----------------------------------------------------------------------
    xY    = DCM.xY;
    n     = length(xY.xy);
    
    try
        t = xY.Time;
    catch
        t = xY.dt*[1:n];
    end

    mx = max(max(cat(2, DCM.K{:})));
    mi = min(min(cat(2, DCM.K{:})));
    
    for i = 1:ns
        subplot(ceil(ns/2),2,i), hold on
        str   = {};
        for j = 1:nt
            plot(t, DCM.K{j}(:,i), ...
                'Color',[1 1 1] - j/nt, ...
                'LineWidth',1);
            str{end + 1} = sprintf('trial %i',j);
        end
        set(gca, 'YLim', [mi mx], 'XLim', [t(1) t(end)]);
        hold off
        title(DCM.Sname{i})
        grid on
        axis square
    end
    xlabel('time (ms)')
    legend(str)
    
    
case{lower('Coupling (A)')}
    
    % spm_dcm_erp_results(DCM,'coupling (A)');
    %----------------------------------------------------------------------
    str = {'Forward','Backward','Lateral'};
    for  i =1:3
        
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
        axis square
    
        % table
        %------------------------------------------------------------------
        subplot(4,3,i + 9)
        text(0,1/2,num2str(DCM.Pp.A{i},' %.2f'),'FontSize',8)
        axis off, axis square
        
    end
    
case{lower('Coupling (C)')}
    
    % spm_dcm_erp_results(DCM,'coupling (C)');
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
    
    % table
    %--------------------------------------------------------------------
    subplot(2,4,2)
    text(0,1/2,num2str(full(exp(DCM.Ep.C)),' %.2f'),'FontSize',8)
    axis off

    % table
    %--------------------------------------------------------------------
    subplot(2,4,4)
    text(0,1/2,num2str(DCM.Pp.C,' %.2f'),'FontSize',8)
    axis off

 
case{lower('Coupling (B)')}
    
    nu = nt-1; % nr of modulatory inputs
    % spm_dcm_erp_results(DCM,'coupling (B)');
    %--------------------------------------------------------------------
    for i = 1:nu
        
        % images
        %-----------------------------------------------------------
        subplot(4,nu,i)
        imagesc(exp(DCM.Ep.B{i}))
        title(DCM.xU.name{i},'FontSize',10)
        set(gca,'YTick',[1:ns],'YTickLabel',DCM.Sname,'FontSize',8)
        set(gca,'XTick',[])
        xlabel('from','FontSize',8)
        ylabel('to','FontSize',8)
        axis square

        % tables
        %--------------------------------------------------------------------
        subplot(4,nu,i + nu)
        text(0,1/2,num2str(full(exp(DCM.Ep.B{i})),' %.2f'),'FontSize',8)
        axis off
        axis square
        
        % PPM
        %-----------------------------------------------------------
        subplot(4,nu,i + 2*nu)
        image(64*DCM.Pp.B{i})
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
    tU    = 0:xU.dt:xU.dur;
    [U N] = spm_erp_u(DCM.Ep,tU);
    
    xY    = DCM.xY;
    n     = length(xY.xy);
    try
        t = xY.Time;
    catch
        t = xY.dt*[1:n];
    end

    subplot(2,1,1)
    plot(t,U,t,N,':')
    xlabel('time (ms)')
    title('input')
    axis square, grid on
    legend({'input','nonspecific'})
    
case{lower('Response')}
    
    % plot data
    % --------------------------------------------------------------------
    xY    = DCM.xY;
    n     = length(xY.xy);
    try
        t = xY.Time;
        E = DCM.M.E*DCM.M.E';
    catch
        t = xY.dt*[1:n];
        E = 1;
    end
    
    for i = 1:n
        subplot(n,2,2*i - 1)
        plot(xY.Time,xY.xy{i}*E)
        xlabel('time (ms)')
        title(sprintf('Observed ERP %i',i))
        axis square, grid on, A = axis;
        try
            if isfield(DCM,'Hc')
                subplot(n,2,2*i - 0)
                plot(xY.Time,DCM.Hc{i})
                xlabel('time (ms)')
                title(sprintf('Predicted ERP %i',i))
                axis(A); axis square, grid on
            end
        end
    end

case{lower('Dipoles')}
    
    try
        P = DCM.Ep;        
        sdip.n_seeds = 1;
        sdip.n_dip  = ns;
        sdip.Mtb    = 1;
        sdip.j{1}   = full(P.Lmom);
        sdip.j{1}   = sdip.j{1}(:);
        sdip.loc{1} = full(P.Lpos);
        spm_eeg_inv_ecd_DrawDip('Init', sdip)
    catch
        warndlg('spm_eeg_inv_visu3D_api to view results')
        return
    end
    case{lower('Spatial overview')}
        spm_dcm_erp_viewspatial(DCM)
end
