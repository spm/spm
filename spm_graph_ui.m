function varargout = spm_graph_ui(action, varargin)
% Graphical display of adjusted data
% FORMAT [Y,y,beta,Bcov] = spm_graph_ui(xSPM,SPM,hReg)
%
% xSPM   - structure containing SPM, distributional & filtering details
%          about the excursion set
% SPM    - structure containing generic details about the analysis
% hReg   - handle of MIP register or [x y z] coordinates
%
% Y      - fitted   data for the selected voxel
% y      - adjusted data for the selected voxel
% beta   - parameter estimates (ML or MAP)
% Bcov   - Covariance of parameter estimates (ML or conditional)
%
% See spm_getSPM for details.
%__________________________________________________________________________
%
% spm_graph is a Callback function that uses the structures above to:
% (1) send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
%
% a) Contrasts of parameter estimates (e.g. activations) and their
% standard error.
%
% b) Fitted and adjusted responses that can be plotted against time, scan,
% or an indicator variable in the design matrix.
%
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses that would have been seen in the absence of other
% effects. The PSTH (peristimulus-time histogram) option provides a finite
% impulse response (FIR) estimate of the trial-specific evoked response as
% a function of peristimulus time.  This is estimated by refitting a
% convolution model to the selected voxel using an FIR basis set.  This is
% simply a set of small boxes covering successive time bins after trial
% onset.  The width of each bin is usually the TR.  This option provides a
% more time-resolved quantitative characterisation of the evoked
% hemodynamic response.  However, it should not be over-interpreted because
% inference is usually made using a simpler and more efficient basis set
% (e.g., canonical hrf, or canonical plus time derivative).
%
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g. in
% constructing explanatory variables such as in a psychophysiological
% interaction). To remove or correct for specific effects, specify an
% appropriate F contrast and simply plot the fitted (and adjusted)
% responses after selecting that F contrast. The vectors Y (fitted) and y
% (adjusted) in the workspace will now be corrected for the effects in the
% reduced design matrix (X0) specified in the contrast manager with the
% column indices (iX0) of the confounds in this adjustment.
%
% Plotting data:
% All data and graphics use filtered/whitened data and residuals. In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrasts previously defined can be plotted. This
% ensures that the parameters plotted are meaningful even when there is
% collinearity among the design matrix subpartitions.
%
% Selecting contrasts used for PPMs will automatically give plots
% based on conditonal estimates.
%
% The structure     contrast.contrast      = cbeta;
%                   contrast.standarderror = SE;
%                   contrast.interval      = 2*CI;
%
% is assigned in base workspace for plots of contrasts and their error.
%__________________________________________________________________________
% Copyright (C) 1996-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston, Guillaume Flandin
% $Id: spm_graph_ui.m 6314 2015-01-23 17:00:51Z guillaume $


if nargin && ~ischar(action)
    varargin = {action varargin{:}};
    action   = 'Graph';
end

%==========================================================================
switch lower(action), case 'graph'                                  %-Graph
%==========================================================================
% FORMAT [Y,y,beta,Bcov] = spm_graph_ui('Graph',xSPM,SPM,hReg)

    if nargin < 3, error('Not enough input arguments.'); end
    [xSPM,SPM,hReg] = deal(varargin{:});

    [xG, XYZ] = spm_graph_ui('UI',xSPM,SPM,hReg);
    if isempty(xG), varargout = {[],[],[],[]}; return; end

    G = spm_graph_ui('Compute',SPM,XYZ,xG);

    XYZmm = SPM.xVol.M(1:3,:) * [XYZ(:);1];
    spm_graph_ui('Plot',SPM,XYZmm,xG,G,hReg);

    varargout = { G.Y, G.y, G.beta, G.Bcov };


%==========================================================================
case 'ui'                                                              %-UI
%==========================================================================
% FORMAT [xG, XYZ] = spm_graph_ui('UI',xSPM,SPM,hReg)

    if nargin < 3, error('Not enough input arguments.'); end
    [xSPM,SPM,hReg] = deal(varargin{:});

    %-Delete previous axis and their pagination controls (if any)
    %----------------------------------------------------------------------
    Fgraph = spm_figure('GetWin','Graphics');
    spm_results_ui('Clear',Fgraph,2);

    %-Find nearest voxel [Euclidean distance] in point list & update GUI
    %----------------------------------------------------------------------
    if isempty(xSPM.XYZmm)
        spm('alert!','No suprathreshold voxels!',mfilename,0);
        varargout = {[],[]};
        return
    end

    if numel(hReg) == 1
        XYZmm = spm_XYZreg('GetCoords',hReg);
    else
        XYZmm = hReg;
    end
    [XYZmm,i] = spm_XYZreg('NearestXYZ',XYZmm,xSPM.XYZmm);
    if numel(hReg) == 1, spm_XYZreg('SetCoords',XYZmm,hReg); end
    XYZ       = xSPM.XYZ(:,i);


    %-Find out what to plot
    %======================================================================
    Cplot   = {'Contrast estimates and 90% C.I.',...
               'Fitted responses',...
               'Event-related responses',...
               'Parametric responses',...
               'Volterra Kernels'};

    % ensure options are appropriate
    %----------------------------------------------------------------------
    if ~isfield(SPM,'Sess'), Cplot = Cplot(1:2); end
    i       = spm_input('Plot',-1,'m',Cplot);
    xG.def  = Cplot{i};

    switch xG.def

        % select contrast if
        %------------------------------------------------------------------
        case {'Contrast estimates and 90% C.I.','Fitted responses'}

            % determine which contrast
            %--------------------------------------------------------------
            xG.spec.Ic = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});

        % select session and trial if
        %------------------------------------------------------------------
        case {'Event-related responses','Parametric responses','Volterra Kernels'}

            % get session
            %--------------------------------------------------------------
            s     = length(SPM.Sess);
            if  s > 1
                s = spm_input('which session','+1','n1',1,s);
            end
            xG.spec.Sess = s;

            % effect names
            %--------------------------------------------------------------
            switch xG.def
                case 'Volterra Kernels'
                    u = length(SPM.Sess(s).Fc);
                otherwise
                    u = length(SPM.Sess(s).U);
            end
            Uname = {};
            for i = 1:u
                Uname{i} = SPM.Sess(s).Fc(i).name;
            end
            if isempty(Uname), error('No conditions found.'); end

            % get effect
            %--------------------------------------------------------------
            str   = 'which effect';
            u     = spm_input(str,'+1','m',Uname);
            xG.spec.u = u;

    end

    switch xG.def
        case 'Fitted responses'

            % predicted or adjusted response
            %--------------------------------------------------------------
            str   = 'predicted or adjusted response?';
            xG.spec.predicted = spm_input(str,'!+1','b',{'predicted','adjusted'},[1 0]);

            % get ordinates
            %--------------------------------------------------------------
            Xplot = {'an explanatory variable',...
                     'scan or time',...
                     'a user specified ordinate'};
            Cx    = spm_input('plot against','!+1','m',Xplot);
            if     Cx == 1
                str  = 'Which explanatory variable?';
                xG.spec.x.i = spm_input(str,'!+1','m',SPM.xX.name);
            elseif Cx == 2
                xG.spec.x.scan = true;
            elseif Cx == 3
                xG.spec.x.x = spm_input('enter ordinate','!+1','e','',numel(SPM.xY.VY));
            end

        %-Modeling evoked responses based on Sess
        %==================================================================
        case 'Event-related responses'

            % get plot type
            %--------------------------------------------------------------
            Rplot = {'fitted response and PSTH',...
                     'fitted response and 90% C.I.',...
                     'fitted response and adjusted data'};

            Rplot = Rplot{spm_input('plot in terms of','+1','m',Rplot)};
            xG.spec.Rplot = Rplot;

        %-Modeling evoked responses based on Sess
        %==================================================================
        case 'Parametric responses'

            % return gracefully if no parameters
            %--------------------------------------------------------------
            if ~SPM.Sess(s).U(u).P(1).h
                warning('No parametric modulation in condition "%s".',...
                    SPM.Sess(s).Fc(u).name);
                varargout = {[],[]};
                return
            end

            str   = 'which parameter';
            p     = spm_input(str,'+1','m',{SPM.Sess(s).U(u).P.name});
            xG.spec.p = p;
    end
    
    varargout = { xG, XYZ };


%==========================================================================
case 'compute'                                                    %-Compute
%==========================================================================
% FORMAT [Y,y,beta,Bcov,G] = spm_graph_ui('Compute',SPM,XYZ,xG)

spm('Pointer','Watch');
[Y,y,beta,Bcov,G] = spm_graph(varargin{:});
spm('Pointer','Arrow');

G.Y = Y; G.y = y; G.beta = beta; G.Bcov = Bcov;
varargout = { G };


%==========================================================================
case 'plot'                                                          %-Plot
%==========================================================================
% FORMAT spm_graph_ui('Plot',SPM,XYZmm,xG,G,hReg)

    [SPM,XYZmm,xG,G,hReg] = deal(varargin{:});

    %-Get Graphics figure handle
    %----------------------------------------------------------------------
    Fgraph = spm_figure('GetWin','Graphics');
    
    %-Scaled font sizes
    %----------------------------------------------------------------------
    FS     = spm('FontSizes');
    
    %-Voxel coordinates
    %----------------------------------------------------------------------
    XYZstr = sprintf(' at [%g, %g, %g]',XYZmm);

    %-Colour specification
    %----------------------------------------------------------------------
    Col    = [0 0 0; .8 .8 .8; 1 .5 .5];

    switch xG.def

        %-Plot parameter estimates
        %==================================================================
        case 'Contrast estimates and 90% C.I.'

            % bar chart
            %--------------------------------------------------------------
            figure(Fgraph)
            ax = subplot(2,1,2,'Parent',Fgraph);
            cla(ax)
            hold(ax,'on')

            % estimates
            %--------------------------------------------------------------
            cbeta = G.contrast;
            h     = bar(ax,cbeta);
            set(h,'FaceColor',Col(2,:))

            % standard error
            %--------------------------------------------------------------
            CI    = G.interval / 2;
            for j = 1:length(cbeta)
                line([j j],([CI(j) -CI(j)] + cbeta(j)),...
                    'LineWidth',6,'Color',Col(3,:),'Parent',ax)
            end

            TTLstr = {xG.def SPM.xCon(xG.spec.Ic).name};
            if isfield(SPM,'VCbeta')
                TTLstr = [TTLstr {'(conditional estimates)'}];
            end
            title(ax,TTLstr,'FontSize',FS(12))
            xlabel(ax,'contrast','FontSize',FS(12))
            ylabel(ax,['contrast estimate',XYZstr],'FontSize',FS(12))
            set(ax,'XLim',[0.4 (length(cbeta) + 0.6)])
            hold(ax,'off')

            % for compatibility with earlier versions of SPM
            %--------------------------------------------------------------
            assignin('base','contrast',G)

        %-All fitted effects or selected effects
        %==================================================================
        case 'Fitted responses'

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)
            ax = subplot(2,1,2,'Parent',Fgraph);
            cla(ax)
            hold(ax,'on')

            x = G.x;
            y = G.y;
            Y = G.Y;
            [p,q] = sort(x);
            if all(diff(x(q)))
                plot(x(q),Y(q),'LineWidth',4,'Color',Col(2,:));
                plot(x(q),y(q),':','Color',Col(1,:));
                plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(3,:));
            else
                plot(x(q),Y(q),'.','MarkerSize',16,'Color',Col(1,:));
                plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
                xlim = get(ax,'XLim');
                xlim = [-1 1]*diff(xlim)/4 + xlim;
                set(ax,'XLim',xlim)
            end

            TTLstr = {xG.def SPM.xCon(xG.spec.Ic).name};
            if isfield(SPM,'VCbeta')
                TTLstr = [TTLstr {'(conditional estimates)'}];
            end
            title(ax,TTLstr,'FontSize',FS(12))
            switch char(fieldnames(xG.spec.x))
                case 'i' % an explanatory variable
                    XLAB = SPM.xX.name{xG.spec.x.i};
                case 'scan' % scan or time
                    if isfield(SPM.xY,'RT')
                        XLAB = 'time {seconds}';
                    else
                        XLAB = 'scan number';
                    end
                case 'x' % user specified
                    XLAB = 'ordinate';
            end
            xlabel(ax,XLAB,'FontSize',FS(12))
            ylabel(ax,['response',XYZstr],'FontSize',FS(12))
            legend(ax,'fitted','plus error')
            hold(ax,'off')

        %-Modeling evoked responses based on Sess
        %==================================================================
        case 'Event-related responses'

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)
            ax = subplot(2,1,2,'Parent',Fgraph);
            hold(ax,'on')

            x = G.x;
            y = G.y;
            Y = G.Y;
            s = xG.spec.Sess;
            u = xG.spec.u;
            Rplot = xG.spec.Rplot;
            if strcmp(Rplot,'fitted response and PSTH') && ~isfield(G,'PSTH')
                % Data not available
                Rplot = 'fitted response and 90% C.I.';
            end
            switch Rplot

                case 'fitted response and PSTH'
                    %------------------------------------------------------
                    PST  = G.PST;
                    PSTH = G.PSTH;
                    PCI  = G.PCI;
                    errorbar(ax,PST,PSTH,PCI)
                    plot(ax,PST,PSTH,'LineWidth',4,'Color',Col(2,:))
                    plot(ax,x,Y,'-.','Color',Col(3,:))

                case 'fitted response and 90% C.I.'
                    %------------------------------------------------------
                    CI   = G.CI;
                    plot(ax,x,Y,'Color',Col(2,:),'LineWidth',4)
                    plot(ax,x,Y + CI,'-.',x,Y - CI,'-.','Color',Col(1,:))

                case 'fitted response and adjusted data'
                    %------------------------------------------------------
                    pst  = G.pst;
                    plot(ax,x,Y,'Color',Col(2,:),'LineWidth',4)
                    plot(ax,pst,y,'.','Color',Col(3,:))

            end

            % label
            %--------------------------------------------------------------
            [i,j] = max(Y);
            text(ceil(1.1*x(j)),i,SPM.Sess(s).Fc(u).name,'FontSize',FS(8),'Parent',ax);
            title(ax,Rplot,'FontSize',FS(12))
            xlabel(ax,'peristimulus time {secs}','FontSize',FS(12))
            ylabel(ax,['response',XYZstr],'FontSize',FS(12))
            hold(ax,'off')

        %-Parametric responses
        %==================================================================
        case 'Parametric responses'

            % plot
            %--------------------------------------------------------------
            figure(Fgraph)

            pst = G.pst;
            P   = G.P;
            Y   = G.Y;
            s   = xG.spec.Sess;
            u   = xG.spec.u;
            p   = xG.spec.p;

            ax  = subplot(2,2,3,'Parent',Fgraph);
            surf(ax,pst,P,Y')
            shading(ax,'flat')
            title(ax,SPM.Sess(s).U(u).name{1},'FontSize',FS(12))
            xlabel(ax,'PST {secs}','FontSize',FS(12))
            ylabel(ax,SPM.Sess(s).U(u).P(p).name,'FontSize',FS(12))
            zlabel(ax,['responses',XYZstr],'FontSize',FS(12))
            axis(ax,'square')

            ax = subplot(2,2,4,'Parent',Fgraph);
            [j,i] = max(mean(Y,2));
            plot(ax,P,Y(i,:),'LineWidth',4,'Color',Col(2,:))
            str   = sprintf('response at %0.1fs',i*SPM.xBF.dt);
            title(ax,str,'FontSize',FS(12))
            xlabel(ax,SPM.Sess(s).U(u).P(p).name,'FontSize',FS(12))
            axis(ax,'square')
            grid(ax,'on')

        %-Volterra Kernels
        %==================================================================
        case 'Volterra Kernels'

            pst = G.pst;
            Y   = G.Y;
            s   = xG.spec.Sess;
            u   = xG.spec.u;

            % second order kernel
            %--------------------------------------------------------------
            if u > length(SPM.Sess(s).U)

                % plot
                %----------------------------------------------------------
                figure(Fgraph)

                ax = subplot(2,2,3,'Parent',Fgraph);
                imagesc(pst,pst,Y,'Parent',ax)
                axis(ax,'xy')
                axis(ax,'image')
                title(ax,'2nd order Kernel','FontSize',FS(12))
                xlabel(ax,'peristimulus time {secs}','FontSize',FS(12))
                ylabel(ax,'peristimulus time {secs}','FontSize',FS(12))

                ax = subplot(2,2,4,'Parent',Fgraph);
                plot(ax,pst,Y)
                axis(ax,'square')
                grid(ax,'on')
                title(ax,SPM.Sess(s).Fc(u).name,'FontSize',FS(12))
                xlabel(ax,'peristimulus time {secs}','FontSize',FS(12))

            % first  order kernel
            %--------------------------------------------------------------
            else

                % plot
                %----------------------------------------------------------
                figure(Fgraph)
                ax = subplot(2,1,2,'Parent',Fgraph);
                plot(ax,pst,Y)
                grid(ax,'on')
                axis(ax,'square')
                title(ax,{'1st order Volterra Kernel' SPM.Sess(s).Fc(u).name},...
                    'FontSize',FS(12))
                xlabel(ax,'peristimulus time {secs}','FontSize',FS(12))
                ylabel(ax,['impulse response',XYZstr],'FontSize',FS(12))

            end

    end


    %-Turn hold button off - this will alert the user to press it again
    %----------------------------------------------------------------------
    try, set(findobj('Tag','holdButton'),'Value',0); end

    %-Call Plot UI
    %----------------------------------------------------------------------
    spm_results_ui('PlotUi',ax);

    %-Setup registry
    %----------------------------------------------------------------------
    hPB = findobj('Tag','plotButton');
    if ~isempty(hPB)
        setappdata(hPB,'Reg',struct('xG',xG,'Fgraph',Fgraph));
        spm_XYZreg('Add2Reg',hReg,hPB,'spm_graph_ui');
        set(ax,'DeleteFcn',[get(ax,'DeleteFcn') ';' ...
            'spm_XYZreg(''Del2Reg'',hReg,findobj(''Tag'',''plotButton''));']);
    end
    
%==========================================================================
case 'setcoords'                                        %-Coordinate change
%==========================================================================
% FORMAT spm_graph_ui('SetCoords',XYZmm,h,hReg)
    if nargin<3, error('Not enough input arguments.'); end
    XYZmm = varargin{1};
    h     = varargin{2};
    hReg  = varargin{3};
    UD    = getappdata(h,'Reg');

    SPM   = evalin('base','SPM');
    XYZ   = SPM.xVol.iM(1:3,:) * [XYZmm(:);1];
    G     = spm_graph_ui('Compute',SPM,XYZ,UD.xG);

    spm_results_ui('Clear',UD.Fgraph,2);
    
    spm_graph_ui('Plot',SPM,XYZmm,UD.xG,G,hReg);


%==========================================================================
otherwise                                           %-Unknown action string
%==========================================================================
    error('Unknown action string')

end
