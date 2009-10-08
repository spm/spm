function varargout = spm_list(varargin)
% Display and analysis of SPM{.}
% FORMAT TabDat = spm_list('List',SPM,hReg,[Num,Dis,Str])
% Summary list of local maxima for entire volume of interest
% FORMAT TabDat = spm_list('ListCluster',SPM,hReg,[Num,Dis,Str])
% List of local maxima for a single suprathreshold cluster
%
% SPM    - structure containing SPM, distribution & filtering details
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests
% .STAT  - distribution {Z, T, X or F}
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}
% .M     - voxels - > mm matrix
% .VOX   - voxel dimensions {mm}
% .Vspm  - mapped statistic image(s)
% .Ps    - uncorrected P values in searched volume (for voxel FDR)
% .Pp    - uncorrected P values of peaks (for peak FDR)
% .Pc    - uncorrected P values of cluster extents (for cluster FDR)
% .uc    - 0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
% .thresDesc - description of height threshold (string)
%
% (see spm_getSPM for further details of xSPM structures)
%
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% Num    - number of maxima per cluster
% Dis    - distance among clusters (mm)
% Str    - header string
%
% TabDat - Structure containing table data
%        - fields are
% .tit   - table Title (string)
% .hdr   - table header (2x12 cell array)
% .fmt   - fprintf format strings for table data (1x12 cell array)
% .str   - table filtering note (string)
% .ftr   - table footnote information (5x2 cell array)
% .dat   - table data (Nx12 cell array)
%
%                           ----------------
%
% FORMAT spm_list('TxtList',TabDat,c)
% Prints a tab-delimited text version of the table
% TabDat - Structure containing table data (format as above)
% c      - Column of table data to start text table at
%          (E.g. c=3 doesn't print set-level results contained in columns 1 & 2)
%                           ----------------
%
% FORMAT spm_list('SetCoords',xyz,hAx,hC)
% Highlighting of table co-ordinates (used by results section registry)
% xyz    - 3-vector of new co-ordinate
% hAx    - table axis (the registry object for tables)
% hReg   - Handle of caller (not used)
%__________________________________________________________________________
%
% spm_list characterizes SPMs (thresholded at u and k) in terms of
% excursion sets (a collection of face, edge and vertex connected
% subsets or clusters).  The corrected significance of the results are
% based on set, cluster and voxel-level inferences using distributional
% approximations from the Theory of Gaussian Fields.  These
% distributions assume that the SPM is a reasonable lattice
% approximation of a continuous random field with known component field
% smoothness.
%
% The p values are based on the probability of obtaining c, or more,
% clusters of k, or more, resels above u, in the volume S analysed =
% P(u,k,c).  For specified thresholds u, k, the set-level inference is
% based on the observed number of clusters C, = P(u,k,C).  For each
% cluster of size K the cluster-level inference is based on P(u,K,1)
% and for each voxel (or selected maxima) of height U, in that cluster,
% the voxel-level inference is based on P(U,0,1).  All three levels of
% inference are supported with a tabular presentation of the p values
% and the underlying statistic:
%
% Set-level     - c    = number of suprathreshold clusters
%               - P    = prob(c or more clusters in the search volume)
%
% Cluster-level - k    = number of voxels in this cluster
%               - Pc   = prob(k or more voxels in the search volume)
%               - Pu   = prob(k or more voxels in a cluster)
%               - Qc   = lowest FDR bound for which this cluster would be
%                        declared positive
%
% Peak-level    - T/F  = Statistic upon which the SPM is based
%               - Ze   = The equivalent Z score - prob(Z > Ze) = prob(t > T)
%               - Pc   = prob(Ze or higher in the search volume)
%               - Qp   = lowest FDR bound for which this peak would be
%                        declared positive
%               - Pu   = prob(Ze or higher at that voxel)
%
% Voxel-level   - Qu   = Expd(Prop of false positives among voxels >= Ze)
%
% x,y,z (mm)    - Coordinates of the voxel
%
% The table is grouped by regions and sorted on the Ze-variate of the
% primary maxima.  Ze-variates (based on the uncorrected p value) are the
% Z score equivalent of the statistic. Volumes are expressed in voxels.
%
% Clicking on values in the table returns the value to the Matlab
% workspace. In addition, clicking on the co-ordinates jumps the
% results section cursor to that location. The table has a context menu
% (obtained by right-clicking in the background of the table),
% providing options to print the current table as a text table, or to
% extract the table data to the Matlab workspace.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Andrew Holmes
% $Id: spm_list.m 3450 2009-10-08 16:11:13Z guillaume $


% satellite figure global variable
%--------------------------------------------------------------------------
global SatWindow

% Choose between voxel-wise and topological FDR
%--------------------------------------------------------------------------
defaults = spm('GetGlobal','defaults');
try
    topoFDR = defaults.stats.topoFDR;
catch
    topoFDR = true;
end

%==========================================================================
switch lower(varargin{1}), case 'list'                            %-List
%==========================================================================
% FORMAT TabDat = spm_list('list',SPM,hReg)

    %-Tolerance for p-value underflow, when computing equivalent Z's
    %----------------------------------------------------------------------
    tol = eps*10;

    %-Parse arguments and set maxima number and separation
    %----------------------------------------------------------------------
    if nargin < 2,  error('insufficient arguments'),     end
    if nargin < 3,  hReg = []; else  hReg = varargin{3}; end


    %-Get current location (to highlight selected voxel in table)
    %----------------------------------------------------------------------
    xyzmm     = spm_results_ui('GetCoords');

    %-Extract data from xSPM
    %----------------------------------------------------------------------
    S     = varargin{2}.S;
    VOX   = varargin{2}.VOX;
    DIM   = varargin{2}.DIM;
    n     = varargin{2}.n;
    STAT  = varargin{2}.STAT;
    df    = varargin{2}.df;
    u     = varargin{2}.u;
    M     = varargin{2}.M;
    k     = varargin{2}.k;
    try, QPs = varargin{2}.Ps; end
    try, QPp = varargin{2}.Pp; end
    try, QPc = varargin{2}.Pc; end
    try
        thresDesc = sprintf('{%s}', varargin{2}.thresDesc);
    catch
        thresDesc = '';
    end
    
    if STAT~='P'
        R     = varargin{2}.R;
        FWHM  = varargin{2}.FWHM;
    end
    try
        units = varargin{2}.units;
    catch
        units = {'mm' 'mm' 'mm'};
    end
    units{1}  = [units{1} ' '];
    units{2}  = [units{2} ' '];

    DIM       = DIM > 1;              % dimensions
    VOX       = VOX(DIM);             % scaling

    if STAT~='P'
        FWHM  = FWHM(DIM);            % Full width at max/2
        FWmm  = FWHM.*VOX;            % FWHM {units}
        v2r   = 1/prod(FWHM);         % voxels to resels
        k     = k*v2r;                % extent threshold in resels
        R(find(~DIM) + 1) = [];       % eliminate null resel counts
        try, QPs = sort(QPs(:)); end  % Needed for voxel FDR
        try, QPp = sort(QPp(:)); end  % Needed for peak FDR
        try, QPc = sort(QPc(:)); end  % Needed for cluster FDR
    end

    %-get number and separation for maxima to be reported
    %----------------------------------------------------------------------
    if length(varargin) > 3
        Num    = varargin{4};         % number of maxima per cluster
        Dis    = varargin{5};         % distance among clusters (mm)
    else
        Num    = 3;
        Dis    = 8;
    end
    if length(varargin) > 5
        Title  = varargin{6};
    else
        Title  = 'p-values adjusted for search volume';
    end

    %-Setup graphics panel
    %----------------------------------------------------------------------
    spm('Pointer','Watch')
    if SatWindow
        Fgraph = SatWindow;
        figure(Fgraph);
    else
        Fgraph = spm_figure('GetWin','Graphics');
    end
    spm_results_ui('Clear',Fgraph)
    FS    = spm('FontSizes');           %-Scaled font sizes
    PF    = spm_platform('fonts');      %-Font names (for this platform)


    %-Table header & footer
    %======================================================================

    %-Table axes & Title
    %----------------------------------------------------------------------
    if SatWindow, ht = 0.85; bot = 0.14; else ht = 0.4; bot = 0.1; end

    if STAT == 'P'
        Title = 'Posterior Probabilities';
    end

    hAx   = axes('Position',[0.025 bot 0.9 ht],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Units','points',...
                'Visible','off');

    AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])
    dy    = FS(9);
    y     = floor(AxPos(4)) - dy;

    text(0,y,['Statistics:  \it\fontsize{',num2str(FS(9)),'}',Title],...
              'FontSize',FS(11),'FontWeight','Bold');   y = y - dy/2;
    line([0 1],[y y],'LineWidth',3,'Color','r'),        y = y - 9*dy/8;

    %-Construct table header
    %----------------------------------------------------------------------
    set(gca,'DefaultTextFontName',PF.helvetica,'DefaultTextFontSize',FS(8))

    Hc = [];
    Hp = [];
    h  = text(0.01,y,   'set-level','FontSize',FS(9));      Hc = [Hc,h];
    h  = line([0,0.11],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r'); Hc = [Hc,h];
    h  = text(0.08,y-9*dy/8,    '\itc ');                   Hc = [Hc,h];
    h  = text(0.02,y-9*dy/8,    '\itp ');                   Hc = [Hc,h];
    Hp = [Hp,h];
    text(0.22,y,        'cluster-level','FontSize',FS(9));
    line([0.14,0.44],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.15,y-9*dy/8,    '\itp\rm_{FWE-corr}');     Hp = [Hp,h];
    h  = text(0.24,y-9*dy/8,    '\itq\rm_{FDR-corr}');     Hp = [Hp,h];
    h  = text(0.39,y-9*dy/8,    '\itp\rm_{uncorr}');       Hp = [Hp,h];
    h  = text(0.34,y-9*dy/8,    '\itk\rm_E');

    text(0.64,y,        'peak-level','FontSize',FS(9));
    line([0.48,0.88],[1,1]*(y-dy/4),'LineWidth',0.5,'Color','r');
    h  = text(0.49,y-9*dy/8,    '\itp\rm_{FWE-corr}');     Hp = [Hp,h];
    h  = text(0.58,y-9*dy/8,        '\itq\rm_{FDR-corr}'); Hp = [Hp,h];
    h  = text(0.82,y-9*dy/8,    '\itp\rm_{uncorr}');       Hp = [Hp,h];
    h  = text(0.67,y-9*dy/8,     sprintf('\\it%c',STAT));
    h  = text(0.75,y-9*dy/8,    '(\itZ\rm_\equiv)');

    text(0.92,y - dy/2,[units{:}],'Fontsize',FS(8));


    %-Headers for text table...
    %-----------------------------------------------------------------------
    TabDat.tit = Title;
    TabDat.hdr = {  'set',      'c';...
        'set',      'p';...
        'cluster',  'p(FWE-cor)';...
        'cluster',  'p(FDR-cor)';...
        'cluster',  'equivk';...
        'cluster',  'p(unc)';...
        'peak',     'p(FWE-cor)';...
        'peak',     'p(FDR-cor)';...
        'peak',      STAT;...
        'peak',     'equivZ';...
        'peak',     'p(unc)';...
        '',         'x,y,z {mm}'}';...

    TabDat.fmt = {  '%-0.3f','%g',...                          %-Set
        '%0.3f', '%0.3f','%0.0f', '%0.3f',...                  %-Cluster
        '%0.3f', '%0.3f', '%6.2f', '%5.2f', '%0.3f',...        %-Peak
        '%3.0f %3.0f %3.0f'};                                  %-XYZ

    %-Column Locations
    %----------------------------------------------------------------------
    tCol = [ 0.01      0.08 ...                                %-Set
             0.15      0.24      0.33      0.39 ...            %-Cluster
             0.49      0.58      0.65      0.74      0.83 ...  %-Peak
             0.92];                                            %-XYZ

    %-Move to next vertical position marker
    %----------------------------------------------------------------------
    y     = y - 7*dy/4;
    line([0 1],[y y],'LineWidth',1,'Color','r')
    y     = y - 5*dy/4;
    y0    = y;


    %-Table filtering note
    %----------------------------------------------------------------------
    if isinf(Num)
        TabDat.str = sprintf('table shows all local maxima > %.1fmm apart',Dis);
    else
        TabDat.str = sprintf(['table shows %d local maxima ',...
            'more than %.1fmm apart'],Num,Dis);
    end
    text(0.5,4,TabDat.str,'HorizontalAlignment','Center','FontName',PF.helvetica,...
        'FontSize',FS(8),'FontAngle','Italic')


    %-Volume, resels and smoothness (if classical inference)
    %----------------------------------------------------------------------
    line([0 1],[0 0],'LineWidth',1,'Color','r')
    if STAT ~= 'P'
        %------------------------------------------------------------------
        Pz              = spm_P(1,0,u,df,STAT,1,n,S);
        Pu              = spm_P(1,0,u,df,STAT,R,n,S);
        %Qu              = spm_P_FDR(u,df,STAT,n,QPs);
        [P Pn Em En EN] = spm_P(1,k,u,df,STAT,R,n,S);
        
        %-Footnote with SPM parameters
        %------------------------------------------------------------------
        set(gca,'DefaultTextFontName',PF.helvetica,...
            'DefaultTextInterpreter','None','DefaultTextFontSize',FS(8))
        TabDat.ftr    = cell(5,2);
        TabDat.ftr{1} = ...
            sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
            STAT,u,Pz,Pu);
        TabDat.ftr{2} = ...
            sprintf('Extent threshold: k = %0.0f voxels, p = %0.3f (%0.3f)',...
            k/v2r,Pn,P);
        TabDat.ftr{3} = ...
            sprintf('Expected voxels per cluster, <k> = %0.3f',En/v2r);
        TabDat.ftr{4} = ...
            sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
        if any(isnan(varargin{2}.uc))
            TabDat.ftr{5} = ...
            sprintf('FWEp: %0.3f, FDRp: %0.3f',varargin{2}.uc(1:2));
        else
            TabDat.ftr{5} = ...
            sprintf('FWEp: %0.3f, FDRp: %0.3f, FWEc: %0.0f, FDRc: %0.0f',...
            varargin{2}.uc);
        end
        TabDat.ftr{6} = ...
            sprintf('Degrees of freedom = [%0.1f, %0.1f]',df);
        TabDat.ftr{7} = ...
            ['FWHM = ' sprintf('%0.1f ', FWmm) units{:} '; ' ...
            sprintf('%0.1f ', FWHM) '{voxels}'];
        TabDat.ftr{8} = ...
            sprintf('Volume: %0.0f = %0.0f voxels = %0.1f resels', ...
            S*prod(VOX),S,R(end));
        TabDat.ftr{9} = ...
            ['Voxel size: ' sprintf('%0.1f ',VOX) units{:} '; ' ...
            sprintf('(resel = %0.2f voxels)',prod(FWHM))];

        text(0.0,-1*dy,TabDat.ftr{1},...
            'UserData',[u,Pz,Pu],'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-2*dy,TabDat.ftr{2},...
            'UserData',[k/v2r,Pn,P],'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-3*dy,TabDat.ftr{3},...
            'UserData',En/v2r,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-4*dy,TabDat.ftr{4},...
            'UserData',Em*Pn,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.0,-5*dy,TabDat.ftr{5},...
            'UserData',varargin{2}.uc,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-1*dy,TabDat.ftr{6},...
            'UserData',df,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-2*dy,TabDat.ftr{7},...
            'UserData',FWmm,'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-3*dy,TabDat.ftr{8},...
            'UserData',[S*prod(VOX),S,R(end)],...
            'ButtonDownFcn','get(gcbo,''UserData'')')
        text(0.5,-4*dy,TabDat.ftr{9},...
            'UserData',[VOX,prod(FWHM)],...
            'ButtonDownFcn','get(gcbo,''UserData'')')
    else
        TabDat.ftr = {};
    end


    %-Characterize excursion set in terms of maxima
    % (sorted on Z values and grouped by regions)
    %======================================================================
    if isempty(varargin{2}.Z)
        text(0.5,y-6*dy,'no suprathreshold clusters',...
            'HorizontalAlignment','Center',...
            'FontAngle','Italic','FontWeight','Bold',...
            'FontSize',FS(16),'Color',[1,1,1]*.5);
        TabDat.dat = cell(0,12);
        varargout  = {TabDat};
        spm('Pointer','Arrow')
        return
    end

    % Includes Darren Gitelman's code for working around
    % spm_max for conjunctions with negative thresholds
    %----------------------------------------------------------------------
    minz        = abs(min(min(varargin{2}.Z)));
    zscores     = 1 + minz + varargin{2}.Z;
    [N Z XYZ A] = spm_max(zscores,varargin{2}.XYZ);
    Z           = Z - minz - 1;

    %-Convert cluster sizes from voxels to resels
    %----------------------------------------------------------------------
    if STAT~='P'
        if isfield(varargin{2},'VRvp')
            V2R = spm_get_data(varargin{2}.VRvp,XYZ);
        else
            V2R = v2r;
        end
        N       = N.*V2R;
    end

    %-Convert maxima locations from voxels to mm
    %----------------------------------------------------------------------
    XYZmm = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];



    %-Table proper (& note all data in cell array)
    %======================================================================

    %-Pagination variables
    %----------------------------------------------------------------------
    hPage = [];
    set(gca,'DefaultTextFontName',PF.courier,'DefaultTextFontSize',FS(7))


    %-Set-level p values {c} - do not display if reporting a single cluster
    %----------------------------------------------------------------------
    c     = max(A);                                    %-Number of clusters
    if STAT ~= 'P'
        Pc    = spm_P(c,k,u,df,STAT,R,n,S);            %-Set-level p-value
    else
        Pc    = [];
        set(Hp,'Visible','off')
    end

    if c > 1;
        h     = text(tCol(1),y,sprintf(TabDat.fmt{1},Pc),'FontWeight','Bold',...
            'UserData',Pc,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(2),y,sprintf(TabDat.fmt{2},c),'FontWeight','Bold',...
            'UserData',c,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
    else
        set(Hc,'Visible','off')
    end

    TabDat.dat = {Pc,c};            %-Table data
    TabLin     = 1;                 %-Table data line


    %-Local maxima p-values & statistics
    %----------------------------------------------------------------------
    HlistXYZ = [];
    while numel(find(isfinite(Z)))

        % Paginate if necessary
        %------------------------------------------------------------------
        if y < min(Num + 1,3)*dy

            % added Fgraph term to paginate on Satellite window
            %--------------------------------------------------------------
            h     = text(0.5,-5*dy,...
                sprintf('Page %d',spm_figure('#page',Fgraph)),...
                'FontName',PF.helvetica,'FontAngle','Italic',...
                'FontSize',FS(8));

            spm_figure('NewPage',[hPage,h])
            hPage = [];
            y     = y0;
        end

        %-Find largest remaining local maximum
        %------------------------------------------------------------------
        [U,i]   = max(Z);           % largest maxima
        j       = find(A == A(i));  % maxima in cluster


        %-Compute cluster {k} and peak-level {u} p values for this cluster
        %------------------------------------------------------------------
        if STAT ~= 'P'
            Nv      = N(i)/v2r;                       % extent {voxels}
            
            Pz      = spm_P(1,0,   U,df,STAT,1,n,S);  % uncorrected p value
            Pu      = spm_P(1,0,   U,df,STAT,R,n,S);  % FWE-corrected {based on Z}
            [Pk Pn] = spm_P(1,N(i),u,df,STAT,R,n,S);  % [un]corrected {based on k}
            if topoFDR
                Qc  = spm_P_clusterFDR(N(i),df,STAT,R,n,u,QPc); % cluster FDR-corrected {based on k}
                Qp  = spm_P_peakFDR(U,df,STAT,R,n,u,QPp); % peak FDR-corrected {based on Z}
                Qu  = [];
            else
                Qu  = spm_P_FDR(   U,df,STAT,n,QPs);  % voxel FDR-corrected {based on Z}
                Qc  = [];
                Qp  = [];
            end

            if Pz < tol                               % Equivalent Z-variate
                Ze  = Inf;                            % (underflow => can't compute)
            else
                Ze  = spm_invNcdf(1 - Pz);
            end
        else
            Nv      = N(i);
            
            Pz      = [];
            Pu      = [];
            Qu      = [];
            Pk      = [];
            Pn      = [];
            Qc      = [];
            Qp      = [];
            ws      = warning('off','SPM:outOfRangeNormal');
            Ze      = spm_invNcdf(U);
            warning(ws);
        end


        %-Print cluster and maximum peak-level p values {Z}
        %------------------------------------------------------------------
        h     = text(tCol(3),y,sprintf(TabDat.fmt{3},Pk),'FontWeight','Bold',...
            'UserData',Pk,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(4),y,sprintf(TabDat.fmt{4},Qc),'FontWeight','Bold',...
            'UserData',Qc,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(5),y,sprintf(TabDat.fmt{5},Nv),'FontWeight','Bold',...
            'UserData',Nv,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(6),y,sprintf(TabDat.fmt{6},Pn),'FontWeight','Bold',...
            'UserData',Pn,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];

        h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu),'FontWeight','Bold',...
            'UserData',Pu,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        if topoFDR
        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),'FontWeight','Bold',...
            'UserData',Qp,'ButtonDownFcn','get(gcbo,''UserData'')');
        else
        h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu),'FontWeight','Bold',...
            'UserData',Qu,'ButtonDownFcn','get(gcbo,''UserData'')');
        end
        hPage = [hPage, h];
        h     = text(tCol(9),y,sprintf(TabDat.fmt{9},U),'FontWeight','Bold',...
            'UserData',U,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze),'FontWeight','Bold',...
            'UserData',Ze,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];
        h     = ...
            text(tCol(11),y,sprintf(TabDat.fmt{11},Pz),'FontWeight','Bold',...
            'UserData',Pz,'ButtonDownFcn','get(gcbo,''UserData'')');
        hPage = [hPage, h];

        % Specifically changed so it properly finds hMIPax
        %------------------------------------------------------------------
        tXYZmm = XYZmm(DIM,i);
        h     = text(tCol(12),y,sprintf(TabDat.fmt{12},tXYZmm),...
            'FontWeight','Bold',...
            'Tag','ListXYZ',...
            'ButtonDownFcn',[...
            'hMIPax = findobj(''tag'',''hMIPax'');',...
            'spm_mip_ui(''SetCoords'',get(gcbo,''UserData''),hMIPax);'],...
            'Interruptible','off','BusyAction','Cancel',...
            'UserData',XYZmm(:,i));

        HlistXYZ = [HlistXYZ, h];
        if spm_XYZreg('Edist',xyzmm,XYZmm(:,i))<tol && ~isempty(hReg)
            set(h,'Color','r')
        end
        hPage  = [hPage, h];

        y      = y - dy;

        if topoFDR
        [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,Nv,Pn,Pu,Qp,U,Ze,Pz,XYZmm(:,i));
        else
        [TabDat.dat{TabLin,3:12}] = deal(Pk,Qc,Nv,Pn,Pu,Qu,U,Ze,Pz,XYZmm(:,i));
        end
        TabLin = TabLin + 1;

        %-Print Num secondary maxima (> Dis mm apart)
        %------------------------------------------------------------------
        [l q] = sort(-Z(j));                % sort on Z value
        D     = i;
        for i = 1:length(q)
            d = j(q(i));
            if min(sqrt(sum((XYZmm(:,D)-XYZmm(:,d)*ones(1,size(D,2))).^2)))>Dis;

                if length(D) < Num

                    % Paginate if necessary
                    %------------------------------------------------------
                    if y < dy
                        h = text(0.5,-5*dy,sprintf('Page %d',...
                            spm_figure('#page',Fgraph)),...
                            'FontName',PF.helvetica,...
                            'FontAngle','Italic',...
                            'FontSize',FS(8));

                        spm_figure('NewPage',[hPage,h])
                        hPage = [];
                        y     = y0;
                    end

                    % voxel-level p values {Z}
                    %------------------------------------------------------
                    if STAT ~= 'P'
                        Pz    = spm_P(1,0,Z(d),df,STAT,1,n,S);
                        Pu    = spm_P(1,0,Z(d),df,STAT,R,n,S);
                        if topoFDR
                            Qp = spm_P_peakFDR(Z(d),df,STAT,R,n,u,QPp);
                            Qu = [];
                        else
                            Qu = spm_P_FDR(Z(d),df,STAT,n,QPs);
                            Qp = [];
                        end
                        if Pz < tol
                            Ze = Inf;
                        else
                            Ze = spm_invNcdf(1 - Pz); 
                        end
                    else
                        Pz    = [];
                        Pu    = [];
                        Qu    = [];
                        Qp    = [];
                        ws      = warning('off','SPM:outOfRangeNormal');
                        Ze    = spm_invNcdf(Z(d));
                        warning(ws);
                    end

                    h     = text(tCol(7),y,sprintf(TabDat.fmt{7},Pu),...
                        'UserData',Pu,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    if topoFDR
                    h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qp),...
                        'UserData',Qp,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    else
                    h     = text(tCol(8),y,sprintf(TabDat.fmt{8},Qu),...
                        'UserData',Qu,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    end
                    hPage = [hPage, h];
                    h     = text(tCol(9),y,sprintf(TabDat.fmt{9},Z(d)),...
                        'UserData',Z(d),...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];
                    h     = text(tCol(10),y,sprintf(TabDat.fmt{10},Ze),...
                        'UserData',Ze,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];
                    h     = text(tCol(11),y,sprintf(TabDat.fmt{11},Pz),...
                        'UserData',Pz,...
                        'ButtonDownFcn','get(gcbo,''UserData'')');
                    hPage = [hPage, h];

                    % specifically modified line to use hMIPax
                    %------------------------------------------------------
                    tXYZmm = XYZmm(DIM,d);
                    h     = text(tCol(12),y,...
                        sprintf(TabDat.fmt{12},tXYZmm),...
                        'Tag','ListXYZ',...
                        'ButtonDownFcn',[...
                        'hMIPax = findobj(''tag'',''hMIPax'');',...
                        'spm_mip_ui(''SetCoords'',',...
                        'get(gcbo,''UserData''),hMIPax);'],...
                        'Interruptible','off','BusyAction','Cancel',...
                        'UserData',XYZmm(:,d));

                    HlistXYZ = [HlistXYZ, h];
                    if spm_XYZreg('Edist',xyzmm,XYZmm(:,d))<tol && ...
                            ~isempty(hReg)
                        set(h,'Color','r')
                    end
                    hPage = [hPage, h];
                    D     = [D d];
                    y     = y - dy;
                    if topoFDR
                    [TabDat.dat{TabLin,7:12}] = ...
                        deal(Pu,Qp,Z(d),Ze,Pz,XYZmm(:,d));
                    else
                    [TabDat.dat{TabLin,7:12}] = ...
                        deal(Pu,Qu,Z(d),Ze,Pz,XYZmm(:,d));
                    end
                    TabLin = TabLin+1;
                end
            end
        end
        Z(j) = NaN;     % Set local maxima to NaN
    end             % end region


    %-Number and register last page (if paginated)
    %-Changed to use Fgraph for numbering
    %----------------------------------------------------------------------
    if spm_figure('#page',Fgraph)>1
        h = text(0.5,-5*dy,sprintf('Page %d/%d',spm_figure('#page',Fgraph)*[1,1]),...
            'FontName',PF.helvetica,'FontSize',FS(8),'FontAngle','Italic');
        spm_figure('NewPage',[hPage,h])
    end

    %-End: Store TabDat in UserData of axes & reset pointer
    %======================================================================
    h      = uicontextmenu('Tag','TabDat',...
        'UserData',TabDat);
    set(gca,'UIContextMenu',h,...
        'Visible','on',...
        'XColor','w','YColor','w')
    uimenu(h,'Label','Table')
    uimenu(h,'Separator','on','Label','Print text table',...
        'Tag','TD_TxtTab',...
        'CallBack',...
        'spm_list(''txtlist'',get(get(gcbo,''Parent''),''UserData''),3)',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Separator','off','Label','Extract table data structure',...
        'Tag','TD_Xdat',...
        'CallBack','get(get(gcbo,''Parent''),''UserData'')',...
        'Interruptible','off','BusyAction','Cancel');
    uimenu(h,'Separator','on','Label','help',...
        'Tag','TD_Xdat',...
        'CallBack','spm_help(''spm_list'')',...
        'Interruptible','off','BusyAction','Cancel');

    %-Setup registry
    %----------------------------------------------------------------------
    set(hAx,'UserData',struct('hReg',hReg,'HlistXYZ',HlistXYZ))
    spm_XYZreg('Add2Reg',hReg,hAx,'spm_list');

    %-Return TabDat structure & reset pointer
    %----------------------------------------------------------------------
    varargout = {TabDat};
    spm('Pointer','Arrow')


    %======================================================================
    case 'listcluster'                      %-List for current cluster only
    %======================================================================
    % FORMAT TabDat = spm_list('listcluster',SPM,hReg)

        spm('Pointer','Watch')

        %-Parse arguments
        %------------------------------------------------------------------
        if nargin < 2,  error('insufficient arguments'),     end
        if nargin < 3,  hReg = []; else hReg = varargin{3}; end
        SPM    = varargin{2};

        %-get number and separation for maxima to be reported
        %------------------------------------------------------------------
        if length(varargin) > 3

            Num    = varargin{4};       % number of maxima per cluster
            Dis    = varargin{5};       % distance among clusters (mm)
        else
            Num    = 32;
            Dis    = 4;
        end


        %-if there are suprathreshold voxels, filter out all but current cluster
        %------------------------------------------------------------------
        if ~isempty(SPM.Z)

            %-Jump to voxel nearest current location
            %--------------------------------------------------------------
            [xyzmm,i] = spm_XYZreg('NearestXYZ',...
                spm_results_ui('GetCoords'),SPM.XYZmm);
            spm_results_ui('SetCoords',SPM.XYZmm(:,i));

            %-Find selected cluster
            %--------------------------------------------------------------
            A         = spm_clusters(SPM.XYZ);
            j         = find(A == A(i));
            SPM.Z     = SPM.Z(j);
            SPM.XYZ   = SPM.XYZ(:,j);
            SPM.XYZmm = SPM.XYZmm(:,j);
            if isfield(SPM,'Rd'), SPM.Rd = SPM.Rd(:,j); end
        end

        %-Call 'list' functionality to produce table
        %------------------------------------------------------------------
        varargout = {spm_list('list',SPM,hReg,Num,Dis)};


    %======================================================================
    case 'txtlist'                                 %-Print ASCII text table
    %======================================================================
    % FORMAT spm_list('TxtList',TabDat,c)

        if nargin<2, error('Insufficient arguments'), end
        if nargin<3, c=1; else c=varargin{3}; end
        TabDat = varargin{2};

        %-Table Title
        %------------------------------------------------------------------
        fprintf('\n\nSTATISTICS: %s\n',TabDat.tit)
        fprintf('%c',repmat('=',1,80)), fprintf('\n')

        %-Table header
        %------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1,c:end-1}), fprintf('%s\n',TabDat.hdr{1,end})
        fprintf('%s\t',TabDat.hdr{2,c:end-1}), fprintf('%s\n',TabDat.hdr{2,end})
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table data
        %------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for j=c:size(TabDat.dat,2)
                fprintf(TabDat.fmt{j},TabDat.dat{i,j})
                fprintf('\t')
            end
            fprintf('\n')
        end
        for i=1:max(1,12-size(TabDat.dat,1)), fprintf('\n'), end
        fprintf('%s\n',TabDat.str)
        fprintf('%c',repmat('-',1,80)), fprintf('\n')

        %-Table footer
        %------------------------------------------------------------------
        fprintf('%s\n',TabDat.ftr{:})
        fprintf('%c',repmat('=',1,80)), fprintf('\n\n')



        %==================================================================
    case 'setcoords'                                   %-Co-ordinate change
        %==================================================================
        % FORMAT spm_list('SetCoords',xyz,hAx,hReg)
        if nargin<3, error('Insufficient arguments'), end
        hAx      = varargin{3};
        xyz      = varargin{2};
        UD       = get(hAx,'UserData');
        HlistXYZ = UD.HlistXYZ(ishandle(UD.HlistXYZ));

        %-Set all co-ord strings to black
        %------------------------------------------------------------------
        set(HlistXYZ,'Color','k')

        %-If co-ord matches a string, highlight it in red
        %------------------------------------------------------------------
        XYZ      = get(HlistXYZ,'UserData');
        if iscell(XYZ), XYZ = cat(2,XYZ{:}); end
        [null,i,d] = spm_XYZreg('NearestXYZ',xyz,XYZ);
        if d<eps
            set(HlistXYZ(i),'Color','r')
        end

        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string')
end
%==========================================================================
