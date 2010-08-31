function spm_dcm_review(DCM,action)
% Review an estimated DCM
% FORMAT spm_dcm_review(DCM,action)
%
% DCM    - DCM filename
%
% action -
% 'fixed connections'
% ['    effects of ' DCM.U.name{i}];
% 'contrast of connections'
% 'location of regions'
% 'inputs'
% 'outputs'
% 'kernels'
% 'estimates of states'
% 'estimates of parameters'
% 'estimates of precisions'
% ['   hidden states: ' DCM.Y.name{i}]
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_review.m 4057 2010-08-31 15:23:39Z guillaume $


%-Get DCM structure
%==========================================================================
if ~nargin
    [DCM, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
end
if ~isstruct(DCM)
    load(DCM);
end


%-Get model specification structure (see spm_nlsi)
%--------------------------------------------------------------------------
try
    m  = DCM.M.m;                      % number of inputs
    n  = DCM.M.n;                      % number of hidden states
    l  = DCM.M.l;                      % number of regions (responses)
    U  = 0.9;                          % p-value threshold for display
catch
    m  = size(DCM.U.u,2);
    l  = DCM.n;
end


% experimental input specific reports
%----------------------------------------------------------------------
for i = 1:m
    inputstr{i} = ['    effects of ' DCM.U.name{i}];
end

% connectivity and kernels
%----------------------------------------------------------------------
str = [inputstr,{'fixed connections'}];

if isfield(DCM,'averaged')
    str    = {str{:},'contrast of connections'};
else
    str    = {str{:} ,...
            'contrast of connections',...
            'location of regions',...
            'inputs',...
            'outputs',...
            'kernels'};
    if isfield(DCM,'qP')
        str = {str{:} ,...
            'estimates of states',...
            'estimates of parameters',...
            'estimates of precisions'};
        
        % region specific reports
        %------------------------------------------------------------------
        for i = 1:l
            hiddenstr{i} = ['   hidden states: ' DCM.Y.name{i}];
        end
        str = [str,hiddenstr];
    else
        hiddenstr = {};
    end
end

if ~isfield(DCM,'M')
    str = {'location of regions','inputs','outputs'};
end
str = [str,{'quit'}];

%-Get action
%==========================================================================
try
    action;
catch
    
    %-Get action
    %----------------------------------------------------------------------
    action = str{spm_input('review',1,'m',str)};

end


%-Exponentiate coupling parameters if a two-state model
%--------------------------------------------------------------------------
try
    if DCM.options.two_state
        Ep.A = exp(DCM.Ep.A);
        Ep.B = exp(DCM.Ep.B);
        Ep.D = exp(DCM.Ep.D);
        Ep.C = DCM.Ep.C;
        disp('NB: The (A,B,D) parameters are scale parameters')
    else
        Ep.A = DCM.Ep.A;
        Ep.B = DCM.Ep.B;
        Ep.D = DCM.Ep.D;
        Ep.C = DCM.Ep.C;
    end
end

%-Set up Graphics window
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);
set(Fgraph,'Renderer','zbuffer')

switch action


    %======================================================================
    % Inputs
    %======================================================================
    case inputstr

        %-input effects
        %------------------------------------------------------------------
        subplot(2,1,1)
        i        = find(strcmp(action,str));
        b(:,:,1) = DCM.Pp.B(:,:,i);
        b(:,:,2) =     Ep.B(:,:,i);
        c(:,1)   = DCM.Pp.C(:,i);
        c(:,2)   =     Ep.C(:,i);
        spm_dcm_display(DCM.xY,b,c)
        title(sprintf('%s P(coupling > %0.2f)',action,DCM.T),'FontSize',16)


        %-direct effects - connections
        %------------------------------------------------------------------
        subplot(4,2,5)
        bar(c(:,2))
        title('C - direct effects (Hz)','FontSize',14)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        axis square


        %-direct effects - probabilities
        %------------------------------------------------------------------
        subplot(4,2,6)
        P  = c(:,1);
        bar(P), hold on, plot([0 (length(P) + 1)],[U U],'-.r'), hold off
        title('C - probability','FontSize',14)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        ylabel(sprintf('P(C > %0.2f)',DCM.T))
        axis square


        %-modulatory effects - connections
        %------------------------------------------------------------------
        subplot(4,2,7)
        bar(b(:,:,2))
        title('B - modulatory effects {Hz}','FontSize',14)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        xlabel('target region')
        ylabel('strength (Hz)')


        %-modulatory effects - probabilities
        %------------------------------------------------------------------
        subplot(4,2,8)
        P = b(:,:,1);
        bar(P), hold on, plot([0 (length(P) + 1)],[U U],'-.r'), hold off
        title('B - probability','FontSize',14)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        xlabel('target region')
        ylabel(sprintf('P(B > %0.2f)',DCM.T))
        legend(DCM.Y.name)


    %======================================================================
    % Fixed connections
    %======================================================================
    case {'fixed connections'}

        % Intrinsic effects
        %------------------------------------------------------------------
        subplot(2,1,1)
        a(:,:,1) = DCM.Pp.A;
        a(:,:,2) = Ep.A;
        spm_dcm_display(DCM.xY,a)
        title(sprintf('%s P(coupling > %0.2f)','fixed',DCM.T),'FontSize',16)


        % intrinsic interactions
        %------------------------------------------------------------------
        subplot(2,2,3)
        bar(a(:,:,2))
        title('A - fixed effects','FontSize',16)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        xlabel('target region')
        ylabel('strength (Hz)')
        legend(DCM.Y.name)


        % intrinsic interactions - probabilities
        %------------------------------------------------------------------
        subplot(2,2,4)
        P = a(:,:,1);
        bar(P), hold on, plot([0 (length(P) + 1)],[U U],'-.r'), hold off
        title('A - probability','FontSize',16)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        xlabel('target region')
        ylabel(sprintf('P(A > %0.2f)',DCM.T))


    %======================================================================
    % Contrast
    %======================================================================
    case {'contrast of connections'}

        %-get contrast
        %------------------------------------------------------------------
        T    = spm_input('Threshold (Hz)',1,'r',0);
        D    = spm_input('contrast for',2,'b',{'A','B','C'});
        C    = spm_dcm_contrasts(DCM,D);


        %-posterior density and inference
        %------------------------------------------------------------------
        c    = C'*spm_vec(DCM.Ep);
        v    = C'*DCM.Cp*C;
        x    = c + [-512:512]*sqrt(v)*6/512;
        p    = full(1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v)));
        sw   = warning('off','SPM:negativeVariance');
        PP   = 1 - spm_Ncdf(T,c,v);
        warning(sw);

        figure(Fgraph)
        subplot(2,1,1)
        plot(x,p,[1 1]*T,[0 max(p)],'-.');
        str  = sprintf('%s P(contrast > %0.2f) = %.1f%s','Posterior density',T,PP*100,'%');
        title(str,'FontSize',16)
        xlabel('contrast of parameter estimates')
        ylabel('probability density')

        i    = find(x >= T);
        hold on
        fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
        axis square, grid on
        hold off

        %-contrast
        %------------------------------------------------------------------
        C  = spm_unvec(C,DCM.Ep);
        subplot(4,2,5)
        imagesc(C.A)
        title('contrast {A}','FontSize',12)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        axis image

        subplot(4,2,6)
        imagesc(C.C)
        title('contrast {C}','FontSize',12)
        set(gca,'XTick',[1:m],'XTickLabel',DCM.U.name)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        axis image

        for i = 1:m
            subplot(4,m,i + 3*m)
            imagesc(C.B(:,:,i))
            title(['contrast {B}-' DCM.U.name{i}],'FontSize',12)
            set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
            set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
            axis image
        end


    %======================================================================
    % Location of regions
    %======================================================================
    case {'location of regions'}

        % transverse
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
        text(1.0,y,'Voxels',       'FontSize',14)
        text(2.0,y,'Location (mm)','FontSize',14)
        y = y - 1;
        line([0 4],[y y],'LineWidth',4)
        y = y - 1;
        for i = 1:length(DCM.xY)
            name = DCM.xY(i).name;
            N    = length(DCM.xY(i).s);
            L    = DCM.xY(i).xyz;
            r    = DCM.xY(i).spec;
            text(0,y,name, 'FontWeight','bold',        'FontSize',12)
            text(1,y,sprintf('%0.0f',N),               'FontSize',12)
            text(2,y,sprintf('%-4.0f %-4.0f %-4.0f',L),'FontSize',12)
            y = y - 1;
        end
        line([0 4],[y y])
        axis off square


    %======================================================================
    % Inputs
    %======================================================================
    case {'inputs'}

        % priors
        %------------------------------------------------------------------
        x     = (1:length(DCM.U.u))*DCM.U.dt;
        t     = (1:length(DCM.Y.y))*DCM.Y.dt;
        for i = 1:m

            subplot(m,1,i)
            plot(x,full(DCM.U.u(:,i)))

            % posteriors (if stochastic DCM)
            %--------------------------------------------------------------
            try
                hold on
                plot(t,DCM.qU.v{2}(i,:),':')
                
            end
            hold off
            title(['Input - ' DCM.U.name{i}],'FontSize',16)
            ylabel('event density {Hz}')
            axis tight
            lm = get(gca,'ylim'); set(gca,'ylim',lm + [-1 1]*diff(lm)/16);
        end
        xlabel('time {seconds}')
        if DCM.options.stochastic && isfield(DCM,'M')
            legend({'prior','posterior'})
        end


    %======================================================================
    % Outputs
    %======================================================================
    case {'outputs'}

        % graph
        %------------------------------------------------------------------
        x     = [1:DCM.v]*DCM.Y.dt;
        for i = 1:l
            subplot(l,1,i);
            try
                plot(x,DCM.y(:,i),x,DCM.y(:,i) + DCM.R(:,i),':');
                title([DCM.Y.name{i} ': responses and predictions' ],'FontSize',16);
                xlabel('time {seconds}');
                legend('predicted', 'observed')
            catch
                plot(x,DCM.Y.y(:,i));
                title([DCM.Y.name{i} ': responses' ],'FontSize',16);
                xlabel('time {seconds}');
            end
        end


    %======================================================================
    % Kernels
    %======================================================================
    case {'kernels'}

        % input effects
        %------------------------------------------------------------------
        x     = (1:DCM.M.N)*DCM.M.dt;
        d     = 2/DCM.M.dt;
        for i = 1:m

            % input effects - neuronal
            %--------------------------------------------------------------
            y = DCM.K1(:,:,i);
            subplot(m,2,2*(i - 1) + 1)
            plot(x,y)
            set(gca,'XLim',[0 16])
            axis square
            title(['neuronal responses to ' DCM.U.name{i}],'FontSize',12)
            xlabel('time {seconds}')
            for j = 1:l
                text(x(j*d),y(j*d,j),DCM.Y.name{j},...
                    'FontWeight','bold','FontSize',12,...
                    'HorizontalAlignment','Center')
            end

            % input effects - hemodynamic
            %--------------------------------------------------------------
            y = DCM.K1(:,:,i);
            k = DCM.H1(:,:,i);
            subplot(m,2,2*(i - 1) + 2)
            plot(x,k,x,y,':')
            set(gca,'XLim',[0 16])
            axis square
            title('hemodynamic responses','FontSize',12)
            xlabel('time {seconds}')
            for j = 1:l
                text(x(j*d),k(j*d,j),DCM.Y.name{j},...
                    'FontWeight','bold','FontSize',12,...
                    'HorizontalAlignment','Center')
            end

        end

    %======================================================================
    % DEM estimates (using standard format)
    %======================================================================
    case {'estimates of states'}
        spm_DEM_qU(DCM.qU)

    case {'estimates of parameters'}
        spm_DEM_qP(DCM.qP)

    case {'estimates of precisions'}
        spm_DEM_qH(DCM.qH)

    case hiddenstr
        
        
        % get region
        %------------------------------------------------------------------
        i = find(strcmp(action,str)) - find(strcmp('estimates of precisions',str));
        t = [1:length(DCM.Y.y)]*DCM.Y.dt;


        % hidden states - causes
        %------------------------------------------------------------------
        subplot(4,1,1)
        j = find(DCM.c(i,:));
        if isempty(j)
            x = DCM.qU.v{2};
        else
            x = DCM.qU.v{2}(j,:);
        end
        plot(t,x)
        title(['inputs or causes - ' DCM.Y.name{i}],'FontSize',16)


        % hidden states - neuronal
        %------------------------------------------------------------------
        subplot(4,1,2)
        j = DCM.M.x;
        if DCM.options.two_state
            j(i,[1 2 6]) = 1;
        else
            j(i,[1 2]) = 1;
        end
        j = find(j(:));
        x = DCM.qU.x{1}(j,:);
        plot(t,x)
        title('hidden states - neuronal', 'FontSize',16)
        if DCM.options.two_state
            legend({'excitatory','signal','inhibitory'})
        else
            legend({'excitatory','signal'})
        end


        % hidden states - hemodynamic
        %------------------------------------------------------------------
        subplot(4,1,3)
        j = DCM.M.x;
        j(i,[3 4 5]) = 1;
        j = find(j(:));
        x = DCM.qU.x{1}(j,:);
        plot(t,exp(x))
        title('hidden states - hemodynamic', 'FontSize',16)
        legend({'flow','volume','dHb'})


        % response
        %------------------------------------------------------------------
        subplot(4,1,4)
        y  = DCM.qU.v{1}(i,:);
        x  = DCM.qU.v{1}(i,:) + DCM.qU.z{1}(i,:);
        plot(t,x,t,y)
        title('predicted BOLD signal', 'FontSize',16)
        xlabel('time (seconds)')
        legend({'observed','predicted'})

        
    %======================================================================
    % Quit
    %======================================================================
    case {'quit'}
        return;
        
    %======================================================================
    % Unknown action
    %======================================================================
    otherwise
        disp('unknown option')
end

% retutn to menu
%--------------------------------------------------------------------------
if nargin < 2
    spm_dcm_review(DCM);
end
