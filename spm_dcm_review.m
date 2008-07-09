function spm_dcm_review(DCM)
% Review a previously estimated DCM
% FORMAT spm_dcm_review(DCM)
%
% DCM  - DCM filename
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_review.m 1900 2008-07-09 14:57:34Z guillaume $


% Get figure handles
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


% Get results
%--------------------------------------------------------------------------
if ~nargin
    [DCM, sts] = spm_select(1,'^DCM.*\.mat$','select DCM_???.mat');
    if ~sts, return; end
end
if ~isstruct(DCM)
    P = DCM;
    try
        load(DCM);
    catch
        error('Cannot load DCM file "%s".',P);
    end
else
    % DCM can't be a structure as <Contrasts> requires filename.
end


% Get model specification structure (see spm_nlsi)
%--------------------------------------------------------------------------
try
    m  = DCM.M.m;
    n  = DCM.M.n;
    l  = DCM.M.l;
catch
    error('This model has not been estimated yet.');
end


% Check whether the model is an averaged DCM
%--------------------------------------------------------------------------
try
    averaged = DCM.averaged;
catch
    averaged = false;
end
if averaged
    str = {...
        'This model is an "averaged" DCM, containing the results of a Bayesian fixed effects analysis.', ...
        'It only contains averaged parameter estimates and their posterior probabilities.', ...
        ' ', ...
        '<Review> functions are disabled except for computing contrasts.'};
    spm('alert!',str,'DCM <Review>',spm('Cmdline'));
end


% Get threshold and recompute posterior probabilities
%--------------------------------------------------------------------------
str        = 'Threshold {Hz}';
DCM.T      = spm_input(str,1,'e',0,[1 1]);
sw = warning('off','SPM:negativeVariance');
pp         = 1 - spm_Ncdf(DCM.T,abs(DCM.Ep),diag(DCM.Cp));
warning(sw);
[pA pB pC] = spm_dcm_reshape(pp,m,l,1);
DCM.pA     = pA;
DCM.pB     = pB;
DCM.pC     = pC;


% Menu
%--------------------------------------------------------------------------
str = {};
if averaged
    % averaged DCM: only allow for 'contrast' option
    str   = {'Contrast of connections', 'Quit'};
else
    % normal (non-averaged DCMs)
    for i = 1:m
        str{i} = ['    Effects of ' DCM.U.name{i} '    '];
    end
    str{m + 1} =  '    Intrinsic connections';
    str   = {str{:} ,...
        'Contrast of connections',...
        'Location of regions',...
        'Inputs',...
        'Outputs',...
        '1st order Kernels',...
        'Quit'};
end

while true
    
    %-get options
    %----------------------------------------------------------------------
    OPT = spm_input('display',2,'m',str);
    figure(Fgraph);
    spm_figure('Clear',Fgraph);


    %======================================================================
    % Inputs
    %======================================================================
    if (OPT <= m) && ~averaged

        %-input effects
        %------------------------------------------------------------------
        subplot(2,1,1)
        i        = OPT;
        b        = DCM.pB(:,:,i);
        b(:,:,2) = DCM.B(:,:,i);
        c        = DCM.pC(:,i);
        c(:,2)   = DCM.C(:,i);
        spm_dcm_display(DCM.xY,b,c)
        title({ str{OPT};...
            sprintf('P(|connection| > %0.2f)',DCM.T);...
            'connection strength'},'FontSize',12)


        %-direct effects - connections
        %------------------------------------------------------------------
        subplot(4,2,5)
        bar(DCM.C(:,i))
        title('C - direct effects {Hz}','FontSize',12)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        axis square
        grid on

        %-direct effects - probabilities
        %------------------------------------------------------------------
        subplot(4,2,6)
        bar(DCM.pC(:,i))
        title('C {probabilities}','FontSize',12)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        ylabel(sprintf('P(A > %0.2f)',DCM.T))
        axis square
        grid on

        %-modulatory effects - connections
        %------------------------------------------------------------------
        subplot(4,2,7)
        bar3(DCM.B(:,:,i),0.4)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        title('B - modulatory effects {Hz}','FontSize',12)
        xlabel('from')
        ylabel('to')
        axis square

        %-modulatory effects - probabilities
        %------------------------------------------------------------------
        subplot(4,2,8)
        bar3(DCM.pB(:,:,i),0.4)
        title('B {probabilities}','FontSize',12)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        xlabel('from')
        ylabel('to')
        zlabel(sprintf('P(A > %0.2f)',DCM.T))
        axis square
        grid on


    %======================================================================
    % Intrinsic connections
    %======================================================================
    elseif (OPT == 1 + m) && ~averaged

        % Intrinsic effects
        %------------------------------------------------------------------
        subplot(2,1,1)
        a        = DCM.pA;
        a(:,:,2) = DCM.A;
        spm_dcm_display(DCM.xY,a)
        title({ str{OPT};...
            sprintf('P(|connection| > %0.2f)',DCM.T);...
            'connection strength'},'FontSize',10)

        % intrinsic interactions
        %------------------------------------------------------------------
        subplot(2,2,3)
        bar3(DCM.A,0.4)
        title('A - Intrinsic connections','FontSize',12)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        xlabel('from')
        ylabel('to')
        grid on
        axis square

        % intrinsic interactions - probabilities
        %------------------------------------------------------------------
        subplot(2,2,4)
        bar3(DCM.pA,0.4)
        title('A {probabilities}','FontSize',12)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        xlabel('from')
        ylabel('to')
        zlabel(sprintf('P(A > %0.2f)',DCM.T))
        grid on
        axis square


    %======================================================================
    % Contrast
    %======================================================================
    elseif ((OPT == 2 + m) && ~averaged) || ((OPT == 1) && averaged)

        %-get contrast
        %------------------------------------------------------------------
        D       = spm_input('contrast for',1,'b',{'A','B','C'});

        switch D

            % intrinsic connections
            %--------------------------------------------------------------
            case 'A'
                C       = spm_dcm_contrasts(P,'A');
                i       = find(C); j = 1;
                C       = sparse(i + j,1,C(i),length(DCM.Ep),1);

            % modulatory inputs
            %--------------------------------------------------------------
            case 'B' % modulatory inputs
                C       = spm_dcm_contrasts(P,'B');
                i       = find(C); j = 1 + l*l;
                C       = sparse(i + j,1,C(i),length(DCM.Ep),1);

            % direct (driving) inputs
            %--------------------------------------------------------------
            case 'C'
                C       = spm_dcm_contrasts(P,'C');
                i       = find(C); j = 1 + l*l + l*l*m;
                C       = sparse(i + j,1,C(i),length(DCM.Ep),1);

        end

        %-posterior density and inference
        %------------------------------------------------------------------
        c    = C'*DCM.Ep;
        v    = C'*DCM.Cp*C;
        x    = c + [-512:512]*sqrt(v)*6/512;
        p    = full(1/sqrt(2*pi*v)*exp(-[x - c].^2/(2*v)));
        % conversion to full necessary to account for sparse matrices bug in MATLAB 6.5.0 R13
        sw = warning('off','SPM:negativeVariance');
        PP   = 1 - spm_Ncdf(DCM.T,c,v);
        warning(sw);

        figure(Fgraph)
        subplot(2,1,1)
        plot(x,p,[1 1]*DCM.T,[0 max(p)],'-.');
        title({'Posterior density of contrast',...
            sprintf('P(contrast > %0.2f) = %.1f%s',DCM.T,PP*100,'%')},...
            'FontSize',12)
        xlabel('contrast')
        ylabel('probability density')

        i    = find(x >= DCM.T);
        hold on
        fill([x(i) fliplr(x(i))],[i*0 fliplr(p(i))],[1 1 1]*.8)
        axis square, grid on
        hold off

        %-contrast
        %------------------------------------------------------------------
        [A B C] = spm_dcm_reshape(C,m,l,1);
        subplot(4,2,5)
        imagesc(A)
        title('contrast {A}','FontSize',12)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
        axis image

        subplot(4,2,6)
        imagesc(C)
        title('contrast {C}','FontSize',12)
        set(gca,'XTick',[1:m],'XTickLabel',DCM.U.name)
        set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
        axis image

        for i = 1:m
            subplot(4,m,i + 3*m)
            imagesc(B(:,:,i))
            title(['contrast {B}-' DCM.U.name{i}],'FontSize',12)
            set(gca,'XTick',[1:l],'XTickLabel',DCM.Y.name)
            set(gca,'YTick',[1:l],'YTickLabel',DCM.Y.name)
            axis image
        end


    %======================================================================
    % Location of regions
    %======================================================================
    elseif (OPT == 3 + m) && ~averaged

        % transverse
        %------------------------------------------------------------------
        subplot(2,2,3)
        title('transverse')
        u = [[0 1 0];[-1 0 0]]';
        spm_dcm_display(DCM.xY,[],[],u,32);

        % sagittal
        %------------------------------------------------------------------
        subplot(2,2,1)
        title('sagittal')
        u = [[0 1 0];[0 0 1]]';
        spm_dcm_display(DCM.xY,[],[],u,32);

        % coronal
        %------------------------------------------------------------------
        subplot(2,2,2)
        title('coronal')
        u = [[1 0 0];[0 0 1]]';
        spm_dcm_display(DCM.xY,[],[],u,32);

        % table
        %------------------------------------------------------------------
        subplot(2,2,4)
        title(str{OPT},'FontSize',12)
        y = 0;
        line([0 4],[y y])
        y = y - 1;
        text(0.0,y,'Name',         'FontSize',10)
        text(1.0,y,'Voxels',       'FontSize',10)
        text(2.0,y,'Location (mm)','FontSize',10)
        y = y - 1;
        line([0 4],[y y],'LineWidth',4)
        y = y - 1;
        for i = 1:length(DCM.xY)
            name = DCM.xY(i).name;
            N    = length(DCM.xY(i).s);
            L    = DCM.xY(i).xyz;
            r    = DCM.xY(i).spec;
            text(0,y,name,                             'FontSize',8)
            text(1,y,sprintf('%0.0f',N),               'FontSize',8)
            text(2,y,sprintf('%-4.0f %-4.0f %-4.0f',L),'FontSize',8)
            y = y - 1;
        end
        line([0 4],[y y])
        axis off square

    
    %======================================================================
    % Inputs
    %======================================================================
    elseif (OPT == 4 + m) && ~averaged

        % graph
        %------------------------------------------------------------------
        x     = [1:length(full(DCM.U.u))]*DCM.U.dt;
        for i = 1:m
            subplot(m,1,i)
            plot(x,full(DCM.U.u(:,i)))
            title(['Input - ' DCM.U.name{i}],'FontSize',12)
            %axis tight
            ylabel('event density {Hz}')
            if i==m
                xlabel('time {seconds}')
            end
            %axis square tight
            grid on
            % Change axis so we can see boxcars properly
            u_min=min(full(DCM.U.u(:,i)));
            u_max=max(full(DCM.U.u(:,i)));
            u_range=u_max-u_min;
            axis([0 max(x) u_min-0.1*u_range u_max*1.1]);
        end


    %======================================================================
    % Outputs
    %======================================================================
    elseif (OPT == 5 + m) && ~averaged

        % graph
        %------------------------------------------------------------------
        x  = [1:length(DCM.y)]*DCM.Y.dt;
        for i = 1:l
            subplot(l,1,i);
            plot(x,DCM.Y.y(:,i));
            if i==l
                xlabel('time {seconds}');
            end
            grid on
            hold on

            % Remove zero columns from X0 if there are any
            % (although there should'nt be !)
            % Otherwise we won't be able to evaluate inv(X0'X0)
            X0 = DCM.Y.X0;
            X0 = X0(:,any(X0));

            % Note: adding on X0*beta to DCM prediction y
            % is equivalent to orthogonalising y WRT X0
            % because data have already been orth wrt X0
            % (overall_prediction=DCM.y(:,i)+X0*beta;)
            overall_prediction=DCM.y(:,i)-X0*inv(X0'*X0)*X0'*DCM.y(:,i);
            plot(x,overall_prediction,'r');
            title([DCM.Y.name{i} ': data and model predictions' ],'FontSize',12);
            axis normal
        end


    %======================================================================
    % Kernels
    %======================================================================
    elseif (OPT == 6 + m) && ~averaged

        % input effects
        %------------------------------------------------------------------
        x     = [1:DCM.M.N] * DCM.M.dt;
        d     = 2 / DCM.M.dt;
        for i = 1:m

            % input effects - neuronal
            %--------------------------------------------------------------
            y = DCM.K1(:,:,i);
            subplot(m,2,2*(i - 1) + 1)
            plot(x,y)
            set(gca,'XLim',[0 16])
            axis square
            title(['neuronal responses to ' DCM.U.name{i}])
            grid on
            xlabel('time {seconds}')
            for j = 1:l
                text(x(j*d),y(j*d,j),DCM.Y.name{j},...
                    'FontSize',8,...
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
            title('hemodynamic responses')
            grid on
            xlabel('time {seconds}')
            for j = 1:l
                text(x(j*d),k(j*d,j),DCM.Y.name{j},...
                    'FontSize',8,...
                    'HorizontalAlignment','Center')
            end
            
        end

        
    %======================================================================
    % Quit
    %======================================================================
    elseif ((OPT == 7 + m) && ~averaged) || ((OPT == 2) && averaged)

        return

    end

end