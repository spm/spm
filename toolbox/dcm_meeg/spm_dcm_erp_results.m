function [DCM] = spm_dcm_erp_results(DCM,Action)
% Results for ERP Dynamic Causal Modeling (DCM)
% FORMAT spm_dcm_erp_results(DCM,'ERPs (mode)');
% FORMAT spm_dcm_erp_results(DCM,'ERPs (sources)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (A)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (B)');
% FORMAT spm_dcm_erp_results(DCM,'Coupling (C)');
% FORMAT spm_dcm_erp_results(DCM,'trial-specific effects');
% FORMAT spm_dcm_erp_results(DCM,'Input');
% FORMAT spm_dcm_erp_results(DCM,'Response');
% FORMAT spm_dcm_erp_results(DCM,'Response (image)');
% FORMAT spm_dcm_erp_results(DCM,'Data');
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
% $Id: spm_dcm_erp_results.m 4367 2011-06-15 17:04:45Z vladimir $


% get Action if necessary
%--------------------------------------------------------------------------
if nargin < 2
    
    str{1} = 'ERPs (mode)';
    str{2} = 'ERPs (sources)';
    str{3} = 'Coupling (A)';
    str{4} = 'Coupling (B)';
    str{5} = 'Coupling (C)';
    str{6} = 'trial-specific effects';
    str{7} = 'Input';
    str{8} = 'Response';
    str{9} = 'Response (image)';
    str{10} = 'Scalp maps';
    str{11} = 'Data';
    
    s      = listdlg('PromptString','Select an option:',...
        'SelectionMode','single',...
        'ListString',str);
    
    Action = str{s};
    
end

% get figure handle
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
colormap(gray)
figure(Fgraph)
clf

% trial data
%--------------------------------------------------------------------------
xY  = DCM.xY;                   % data
nt  = length(xY.xy);            % Nr trial types
ne  = size(xY.xy{1},2);         % Nr electrodes
nb  = size(xY.xy{1},1);         % Nr time bins
t   = xY.pst;                   % PST

% plot data
%--------------------------------------------------------------------------
switch(lower(Action))
    case{lower('Data')}
        try
            A = [0 0];
            for i = 1:nt
                
                % confounds if specified
                %--------------------------------------------------------------
                try
                    X0 = spm_orth(xY.X0(1:nb,:),'norm');
                    R  = speye(nb,nb) - X0*X0';
                catch
                    R  = speye(nb,nb);
                end
                
                % plot data
                %--------------------------------------------------------------
                subplot(nt,2,(i - 1)*2 + 1)
                plot(t,R*xY.xy{i})
                xlabel('time (ms)')
                try
                    title(sprintf('Observed response (code:%i)',xY.code(i)))
                catch
                    title(sprintf('Observed response %i',i))
                end
                axis square
                a    = axis;
                A(1) = min(A(1),a(3));
                A(2) = max(A(2),a(4));
                
                % image data
                %--------------------------------------------------------------
                subplot(nt,2,(i - 1)*2 + 2)
                imagesc([1:ne],t,R*xY.xy{i})
                xlabel('channels');ylabel('peri-stimulus time (ms)')
                axis square
                try
                    title(sprintf('Observed response (code:%i)',xY.code(i)))
                catch
                    title(sprintf('Observed response %i',i))
                end
            end
            
            % set axis
            %------------------------------------------------------------------
            for i = 1:nt
                subplot(nt,2,(i - 1)*2 + 1)
                set(gca,'YLim',A)
            end
        end
        return
end

% post inversion parameters
%--------------------------------------------------------------------------
nu  = length(DCM.B);          % Nr inputs
nc  = size(DCM.H{1},2);       % Nr modes
ns  = size(DCM.A{1},1);       % Nr of sources
np  = size(DCM.K{1},2)/ns;    % Nr of population per source

% switch
%--------------------------------------------------------------------------
switch(lower(Action))
    
    case{lower('ERPs (mode)')}
        
        % spm_dcm_erp_results(DCM,'ERPs (mode)');
        %----------------------------------------------------------------------
        co = {'b', 'r', 'g', 'm', 'y', 'k'};
        lo = {'-', '--'}; A = [0 0];
        
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
            title(sprintf('mode %i',i))
            grid on
            axis square
            a    = axis;
            A(1) = min(A(1),a(3));
            A(2) = max(A(2),a(4));
        end
        xlabel('time (ms)')
        legend(str)
        
        % set axis
        %----------------------------------------------------------------------
        for i = 1:nc
            subplot(ceil(nc/2),2,i)
            set(gca,'YLim',A)
        end
        
    case{lower('ERPs (sources)')}
        
        % spm_dcm_erp_results(DCM,'ERPs (sources)');
        %----------------------------------------------------------------------
        col   = {'b','r','g','m','y','c'}; A = [0 0];
        for i = 1:ns
            str   = {};
            subplot(ceil(ns/2),2,i), hold on
            
            % if J maps from states to sources, use J (a matrix)
            %------------------------------------------------------------------
            if strcmpi(DCM.options.model,'DEM')
                
                for j = find(DCM.Eg.J(i,:))
                    for k = 1:nt
                        plot(t, DCM.K{k}(:,j) , ...
                            'Color',col{k}, ...
                            'LineWidth',1);
                        str{end + 1} = sprintf('trial %i (pop. %i)',k,j);
                    end
                end
                
            else
                
                % otherwise assume normal form for states (source x states)
                %--------------------------------------------------------------
                for j = 1:np
                    for k = 1:nt
                        if j == np
                            plot(t, DCM.K{k}(:,i + ns*(j - 1)), ...
                                'Color',col{k}, ...
                                'LineWidth',2);
                        else
                            plot(t, DCM.K{k}(:,i + ns*(j - 1)), ':', ...
                                'Color',col{k}, ...
                                'LineWidth',2);
                        end
                        str{end + 1} = sprintf('trial %i (pop. %i)',k,j);
                    end
                end
            end
            
            hold off
            title(DCM.Sname{i},'FontSize',16)
            grid on
            axis square
            a    = axis;
            A(1) = min(A(1),a(3));
            A(2) = max(A(2),a(4));
        end
        xlabel('time (ms)','FontSize',14)
        legend(str)
        
        % set axis
        %----------------------------------------------------------------------
        for i = 1:ns
            subplot(ceil(ns/2),2,i)
            axis([t(1) t(end) A(1) A(2)]);
        end
        
    case{lower('Coupling (A)')}
        
        % spm_dcm_erp_results(DCM,'coupling (A)');
        %----------------------------------------------------------------------
        if ~isfield(DCM.Ep,'A'), return, end
        str   = {'Forward','Backward','Lateral'};
        
        for i = 1:length(DCM.Ep.A)
            
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
        
        % spm_dcm_erp_results(DCM,'coupling (C)');
        %----------------------------------------------------------------------
        if ~isfield(DCM.Ep,'C'), return, end
        
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
        
        % spm_dcm_erp_results(DCM,'coupling (B)');
        %----------------------------------------------------------------------
        if ~isfield(DCM.Ep,'B'), return, end
        
        for i = 1:length(DCM.Ep.B)
            
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
        
        % spm_dcm_erp_results(DCM,'trial-specific effects');
        %----------------------------------------------------------------------
        if ~isfield(DCM.Ep,'B'), return, end
        
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
        
        % plot data
        % ---------------------------------------------------------------------
        xU    = DCM.xU;
        tU    = t - t(1);
        U     = spm_erp_u(tU/1000,DCM.Ep,DCM.M);
        
        subplot(2,1,1)
        plot(t,U)
        xlabel('time (ms)')
        title('input')
        axis square, grid on
        for i = 1:length(DCM.M.ons)
            str{i} = sprintf('input (%i)',i);
        end
        legend(str)
        
    case{lower('Response')}
        
        % get spatial projector
        % ---------------------------------------------------------------------
        try
            U = DCM.M.E';
        catch
            U = 1;
        end
        
        % plot data
        % ---------------------------------------------------------------------
        try
            A     = [];
            for i = 1:nt
                subplot(nt,2,2*i - 1)
                plot(t,(DCM.H{i} + DCM.R{i})*U)
                xlabel('time (ms)')
                try
                    title(sprintf('Observed (adjusted-code:%i)',xY.code(i)))
                catch
                    title(sprintf('Observed (adjusted) %i',i))
                end
                A(end + 1,:) = axis;
                
                subplot(nt,2,2*i - 0)
                plot(t,DCM.H{i}*U)
                xlabel('channels')
                ylabel('time (ms)')
                title('Predicted')
                A(end + 1,:) = axis;
            end
            a(1)  = min(A(:,1));
            a(2)  = max(A(:,2));
            a(3)  = min(A(:,3));
            a(4)  = max(A(:,4));
            for i = 1:nt
                subplot(nt,2,2*i - 1)
                axis(a); axis square, grid on
                subplot(nt,2,2*i - 0)
                axis(a); axis square, grid on
            end
            
        end
        
        
    case{lower('Response (image)')}
        
        % get spatial projector
        % ---------------------------------------------------------------------
        try
            U = DCM.M.E';
        catch
            U = 1;
        end
        
        
        % plot data
        % ---------------------------------------------------------------------
        try
            for i = 1:nt
                subplot(nt,2,2*i - 1)
                imagesc([1:ne],t,(DCM.H{i} + DCM.R{i})*U)
                xlabel('time (ms)')
                try
                    title(sprintf('Observed (adjusted-code:%i)',xY.code(i)))
                catch
                    title(sprintf('Observed (adjusted) %i',i))
                end
                axis square, grid on, A = axis;
                
                subplot(nt,2,2*i - 0)
                imagesc([1:ne],t,DCM.H{i}*U)
                xlabel('channels')
                ylabel('time (ms)')
                title('Predicted')
                axis(A); axis square, grid on
            end
        end
        
        
    case{lower('Scalp maps')}
        
        % get spatial projector
        % ---------------------------------------------------------------------
        try
            U = DCM.M.E';
        catch
            U = 1;
        end
        
        try
            pos = DCM.xY.coor2D;
        catch
            [xy, label] = spm_eeg_project3D(DCM.M.dipfit.sens, DCM.xY.modality);
            [sel1, sel2] = spm_match_str(DCM.xY.name, label);
            pos = xy(:, sel2);
        end
        
        ns = 5; %number of time frames
        
        in = [];
        in.type = DCM.xY.modality;
        in.f = Fgraph;
        in.noButtons = 1;
        in.cbar = 0;
        in.plotpos = 0;
        
        % plot data
        % ---------------------------------------------------------------------
        
        for i = 1:nt
            Yo  = (DCM.H{i} + DCM.R{i})*U;
            Yp  = DCM.H{i}*U;
            
            for j = 1:ns
                ind = ((j-1)*floor(nb/ns)+1):j*floor(nb/ns);
                
                in.max = max(abs(mean(Yo(ind, :))));
                in.min = -in.max;
                
                in.ParentAxes = subplot(nt*2,ns,(i - 1)*2*ns + j);
                spm_eeg_plotScalpData(mean(Yo(ind, :))', pos , DCM.xY.name, in);
                title(sprintf('observed\n%s\n%.0f ms', DCM.xY.code{i}, mean(t(ind))));
                
                in.ParentAxes = subplot(nt*2,ns,(i - 1)*2*ns + ns + j);
                spm_eeg_plotScalpData(mean(Yp(ind, :))', pos, DCM.xY.name, in);
                title(sprintf('predicted\n%s\n%.0f ms', DCM.xY.code{i}, mean(t(ind))));
            end
        end
        
        
        
    case{lower('Dipoles')}
        
        
        % return if LFP
        % ---------------------------------------------------------------------
        if strcmpi(DCM.xY.modality,'lfp')
            warndlg('There are no ECDs for these LFP data')
            return
        end
        
        % plot dipoles
        % ---------------------------------------------------------------------
        switch DCM.options.spatial
            
            case{'ECD'}
                P            = DCM.Eg;
                P.L          = spm_cat(P.L);
                np           = size(P.L,2)/size(P.Lpos,2);
                sdip.n_seeds = 1;
                sdip.n_dip   = np*ns;
                sdip.Mtb     = 1;
                sdip.j{1}    = full(P.L);
                sdip.j{1}    = sdip.j{1}./repmat(sqrt(sum(sdip.j{1}.^2)), 3, 1);
                sdip.j{1}    = sdip.j{1}(:);
                sdip.loc{1}  = kron(ones(1,np),full(P.Lpos));
                spm_eeg_inv_ecd_DrawDip('Init', sdip)
                
            case{'LFP'}
                warndlg('This is a LFP spatial model')
                
            case{'IMG'}
                warndlg('use the render API button to see reconstruction')
                
            otherwise
                return
        end
        
end
