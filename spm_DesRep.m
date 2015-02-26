function varargout = spm_DesRep(varargin)
% Design reporting utilities
% FORMAT varargout = spm_DesRep(action,varargin)
%
% spm_DesRep (design reporting) is a suite of utility functions for various
% graphical reports on a given experimental design, embodied in the design
% matrix structure and other associated data structures.
% For detailed programmers comments, see format specifications in main body
% of code.
%
%                           ----------------
%
% By default, spm_DesRep prompts for selection of a SPM.mat file.
%
% Given details of a design spm_DesRep sets up a "Design" menu in the
% SPM 'Interactive' window.  The menu options launch various graphical
% summaries of the current SPM design in the SPM 'Graphics' window.
%
% * Design Matrix  - Displays graphical summary of the design matrix
%
%     The design is labelled with the corresponding parameter and
%     file names, and is displayed as an image scaled (using
%     spm_DesMtx('sca',...) such that zero is mid-grey, -1 is black, and +1
%     is white. Covariates exceeding this randge are scaled to fit.
%
%     The design matrix is "surfable": Clicking (and holding or dragging)
%     around the design matrix image reports the corresponding value of the
%     design matrix ('normal' click - "left" mouse button usually), the
%     image filename ('extend' mouse click - "middle" mouse), or parameter
%     name ('alt' click - "right" mouse). Double clicking the design matrix
%     image extracts the design matrix into the base MATLAB workspace.
%
%     Under the design matrix the parameter estimability is displayed as a
%     1xp matrix of grey and white squares. Parameters that are not
%     uniquely specified by the model are shown with a grey patch. Surfing
%     the estimability image reports the parameter names and their
%     estimability. Double clicking extracts the estimability vector into
%     the base MATLAB workspace.
%
% * Design orthogonality - Displays orthogonality matrix for this design
%
%     The design matrix is displayed as in "Design Matrix" view above,
%     labelled with the parameter names.
%
%     Under the design matrix the design orthogonality matrix is
%     displayed. For each pair of columns of the design matrix, the
%     orthogonality matrix depicts the magnitude of the cosine of the
%     angle between them, with the range 0 to 1 mapped to white to
%     black. Orthogonal vectors (shown in white), have cosine of zero.
%     Colinear vectors (shown in black), have cosine of 1 or -1.
%
%     The cosine of the angle between two vectors a & b is obtained by
%     dividing the dot product of the two vectors by the product of
%     their lengths:
%
%                         a'*b
%                 ------------------------
%                 sqrt(sum(a.^2)*sum(b.^2)
%
%     If (and only if) both vectors have zero mean, i.e.
%     sum(a)==sum(b)==0, then the cosine of the angle between the
%     vectors is the same as the correlation between the two variates.
%     
%     The design orthogonality matrix is "surfable": Clicking (and
%     holding or dragging) the cursor around the design orthogonality
%     image reports the orthogonality of the correponding pair of
%     columns. Double clicking on the orthogonality matrix extracts
%     the contrast orthogonality matrix into the base MATLAB
%     workspace.
%
% * Explore design - Sub-menu's for detailed design exploration.
%
%     If this is an fMRI design, then the session & trial/condition
%     structure of the design is reflected in the sub-menu structure.
%     Selecting a given session, and then trial/condition within the
%     session, launches a comprehensive display of the parameters of
%     that design.
%
%     If not an fMRI design, then the Explore sub-menu has two options:
%     "Files and factors" & "Covariates".
%
% * Explore: Files and factors - Multi-page listing of filenames, 
%                                factor indicies and covariates.
%
%     The covariates printed are the raw covariates as entered into
%     SPM, with the exception of the global value, which is printed
%     after any grand mean scaling.
%
% * Explore: Covariates - Plots of the covariates, showing how they are
%                         included into the model. 
%                  
%     Covariates are plotted, one per page, overlaid on the design
%     matrix. The description strings in the xC covariate structure
%     array are displayed. The corresponding design matrix column(s)
%     is(are) highlighted.
%
% * Clear    - clears Graphics window, re-instating Results section MIP
%              & design matrix graphics (if in the results section).
% * Help     - displays this text!
%
%                           ----------------
%
% spm_DesRep also handles "surfing" of contrast depictions, which are
% bar-graphs for T-contrasts and images for F-contrasts. Clicking
% ('normal' click - "left" mouse button usually) with the on the bars
% of the bar-graphs, or anywhere in an image, and dragging, dynamically
% reports the contrast weight depicted under the cursor. The format of
% the report string is:
%   #{T/F}: <name> (ij) = <weight>
% ...where # is the contrast number, T/F indicates the type of contrast,
% <name> the name given to the contrast, ij the index into the contrast
% vector/matrix weight under the cursor, and <weight> the corresponding
% contrast weight.
%
% Double clicking on a contrast depiction extracts the contrast weights
% into the base workspace.
%__________________________________________________________________________
% Copyright (C) 1999-2015 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_DesRep.m 6351 2015-02-26 16:37:18Z guillaume $


%==========================================================================
% - FORMAT specifications for embedded functions
%==========================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
% FORMAT h = spm_DesRep('DesRepUI',SPM)
% Setup "design" menu for design reporting.
% SPM    - structure containing design details. Required fields are:
%
% .xX    - design matrix structure
%          (See spm_fmri_spm_ui.m & spm_spm.m for formats)
% .xY.P  - (char or cellstr) array of filenames
% .xY.VY - array of file handles (from spm_fmri_spm_ui.m)
% .xM    - Masking structure
% .xC    - Covariate definition structure (not fMRI)
%          (see spm_spm_ui.m for format)
% .Sess  - fMRI session structure (not needed if not fMRI)
%          (see spm_fmri_spm_ui.m for format)
% .xsDes - Design description structure
%          (see spm_fmri_spm_ui.m for details)
% .swd   - SPM working directory - directory where configuration file resides
%          [defaults to empty]
% .SPMid - (recommended) ID string of creator program.
% h      - handle of menu created ('Tag'ged as 'DesRepUI') 
%
%
% FORMAT spm_DesRep('Files&Factors',P,I,xC,sF,xs)
% Produces multi-page listing of files, factor indices, and covariates.
% P   - nxv CellStr of filenames (i.e. reshape(cellstr(SPM.xY.P),size(V)))
% I   - nx4 matrix of factor indices
% xC  - Covariate structure array (see spm_spm_ui.m for definitions)
%       ('rc' & 'cname' fields used)
% sF  - 1x4 CellStr of factor names (i.e. D.sF of spm_spm_ui)
% xs  - [optional] structure of extra strings containing descriptive
%       information which is printed out after the files & variables listing.
%       The field names are used as sub-headings, the field values (which
%       must be strings or CellStr) printed alongside.
%
% FORMAT spm_DesRep('DesMtx',xX,fnames,xs)
% Produces a one-page graphical summary of the design matrix
% FORMAT spm_DesRep('DesOrth',xX)
% Produces a one-page graphical summary of the design orthogonality
% xX      - Structure containing design matrix information
%         - the first of {xX.nX, xX.xKXs.X, xX.X} is used for display
% .nX     - Desgin matrix already scaled for display
% .xKXs.X - temporally filtered design matrix (within space structure)
% .X      - "raw" design matrix (as setup by spm_fmri_spm_ui)
% .name   - [optional] px1 CellStr of parameter names
% fnames  - [optional] nxv CellStr of filenames (i.e. reshape(cellstr(SPM.xY.P),size(V)))
% xs      - [optional] structure of extra strings containing descriptive
%           information which is printed at the foot of the page ('DesMtx' usage)
%           The field names are used as sub-headings, the field values
%           (which must be strings or CellStr) printed alongside.
%
% FORMAT spm_DesRep('fMRIDesMtx',SPM,s,i)
% Interactive review of fMRI design matrix.
% Sess(s).U(i)  -  see spm_fMRI_design for session s, trial i
%
% FORMAT spm_DesRep('Covs',xC,X,Xnames)
% Plots the covariates and describes how they are included into the model.
% xC     - Covariate structure array (see spm_spm_ui.m for details)
%          ('rcname','rc','descrip','cname' & 'cols' fields used)
% X      - nxp Design matrix
% Xnames - px1 CellStr of parameter names
%
% =========================================================================
% Utility functions and CallBack handlers:
%
% FORMAT s = spm_DesRep('ScanTick',nScan,lim)
% Pares down 1:nScan to at most lim items, showing every 2nd/3rd/... as
% necessary to pair  down to <lim items. Always ends with nScan so
% #images is indicated.
% nScan  - number of scans
% lim    - limit to number of elements of s
% s      - 1:nScan pared down accordingly
%
% FORMAT spm_DesRep('SurfDesMtx_CB')
% 'ButtonDownFcn' CallBack for surfing clickable design matrix images
% FORMAT spm_DesRep('SurfDesMtxMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing clickable design matrix images
% FORMAT spm_DesRep('SurfDesMtxUp_CB')
% 'ButtonUpFcn' CallBack for ending surfing of design matrix images
%
% The design matrix, parameter names and image file names should be
% saved in the UserData of the image object as a structure with fields
% 'X','Xnames' & 'fnames' respectively. The code is robust to any of
% these being empty or mis-specified - surfing simply reports "no
% cached data".
%
% FORMAT spm_DesRep('SurfEstIm_CB')
% 'ButtonDownFcn' CallBack for surfing clickable parameter estimability images
% FORMAT spm_DesRep('SurfEstImMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing parameter estimability images
% FORMAT spm_DesRep('SurfEstImUp_CB')
% 'ButtonUpFcn' CallBack for ending surfing of parameter estimability images
%
% The binary parameter estimability matrix and parameter names should
% be saved in the UserData of the image object as a structure with
% fields 'est' & 'Xnames' respectively. The code is robust to any of
% these being empty or mis-specified - surfing simply reports "no
% cached data".
%
% FORMAT spm_DesRep('SurfDesO_CB')
% 'ButtonDownFcn' CallBack for surfing clickable design orthogonality images
% FORMAT spm_DesRep('SurfDesOMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing design orthogonality images
% FORMAT spm_DesRep('SurfDesOUp_CB')
% 'ButtonUpFcn' CallBack for ending surfing of design orthogonality images
%
% The design orthogonality matrix (cosine of angle between vectors),
% correlation index matrix (1 when both columns have zero mean, when
% cos(\theta) equals the correlation), and parameter names should be
% saved in the UserData of the image object as a structure with fields
% 'O', 'bC' & 'Xnames' respectively.
%
% FORMAT spm_DesRep('SurfCon_CB')
% 'ButtonDownFcn' CallBack for surfing clickable contrast depictions
% FORMAT spm_DesRep('SurfConOMo_CB')
% 'WindowButtonMotionFcn' CallBack for surfing contrast depictions
% FORMAT spm_DesRep('SurfConOUp_CB')
% 'ButtonUpFcn' CallBack for ending surfing of contrast depictions
%
% The contrast number, handle of text object to use for reporting and
% contrast structure (for this contrast) should be saved in the UserData
% of the image object as a structure with fields 'i', 'h' & 'xCon'
% respectively.
%__________________________________________________________________________


SVNid = '$Rev: 6351 $'; 

%-Format arguments
%--------------------------------------------------------------------------
if ~nargin
    SPMid = spm('FnBanner',mfilename,SVNid);
    hC    = spm_DesRep('DesRepUI'); 
    SPM   = get(hC,'UserData');
    if ~isempty(SPM)
        cb_menu([],[],'DesMtx',SPM);
        drawnow;
    end
    return;
end

switch lower(varargin{1})

%==========================================================================
case 'desrepui'                                       %-Design reporting UI
%==========================================================================
% h = spm_DesRep('DesRepUI')
% h = spm_DesRep('DesRepUI',SPM)

%-Load design data from file if not passed as argument
%--------------------------------------------------------------------------
if nargin < 2
    [spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
    if ~sts, varargout = {[]}; return; end
    swd = spm_file(spmmatfile,'fpath');
    try
        load(fullfile(swd,'SPM.mat'));
    catch
        error(['Cannot read ' fullfile(swd,'SPM.mat')]);
    end
    SPM.swd = swd;
else
    SPM     = varargin{2};
end

%-Canonicalise data
%--------------------------------------------------------------------------
%-Work out where design configuration has come from
if ~isfield(SPM,'cfg')
    if     isfield(SPM.xX,'V'),  cfg = 'SPMest';
    elseif isfield(SPM.xY,'VY'), cfg = 'SPMdata';
    elseif isfield(SPM,'Sess'),  cfg = 'SPMcfg';
    else   error('Can''t fathom origin!')
    end
    SPM.cfg = cfg;
end

%-Work out what modality this is
%--------------------------------------------------------------------------
try
    SPM.Sess(1);
    SPM.modality = 'fMRI';
catch
    try
        SPM.xC;
    catch
        SPM.xC = {};
    end
    SPM.modality = 'PET';
end

%-Add a scaled design matrix to the design data structure
%--------------------------------------------------------------------------
if ~isfield(SPM.xX,'nKX')
    SPM.xX.nKX = spm_DesMtx('Sca',SPM.xX.X,SPM.xX.name);
end


%-Draw menu
%==========================================================================

%-Get Interactive window and delete any previous DesRepUI menu
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
delete(findobj(get(Finter,'Children'),'flat','Tag','DesRepUI'))

%-Draw top level menu
%--------------------------------------------------------------------------
hC      = uimenu(Finter,'Label','Design',...
        'Separator','on',...
        'Tag','DesRepUI',...
        'UserData',SPM,...
        'HandleVisibility','on');

%-DesMtx
%--------------------------------------------------------------------------
hDesMtx = uimenu(hC,'Label','Design Matrix','Accelerator','D',...
        'CallBack',{@cb_menu,'DesMtx',SPM},...
        'UserData',hC,...
        'HandleVisibility','off');

%-Design matrix orthogonality
%--------------------------------------------------------------------------
h = uimenu(hC,'Label','Design orthogonality','Accelerator','O',...
        'CallBack',{@cb_menu,'DesOrth',SPM},...
        'UserData',hC,...
        'HandleVisibility','off');

%-Explore design
%--------------------------------------------------------------------------
hExplore = uimenu(hC,'Label','Explore','HandleVisibility','off');

switch SPM.modality
case 'PET'
    hFnF = uimenu(hExplore,'Label','Files and factors','Accelerator','F',...
        'CallBack',{@cb_menu,'Files&Factors',SPM},...
        'UserData',hC,...
        'HandleVisibility','off');
    hCovs = uimenu(hExplore,'Label','Covariates','Accelerator','C',...
        'CallBack',{@cb_menu,'Covs',SPM},...
        'UserData',hC,...
        'HandleVisibility','off');
    if isempty(SPM.xC), set(hCovs,'Enable','off'), end
case 'fMRI'
    for j = 1:length(SPM.Sess)
        h = uimenu(hExplore,'Label',sprintf('Session %.0f ',j),...
            'HandleVisibility','off');
        for k = 1:length(SPM.Sess(j).Fc)
            uimenu(h,'Label',SPM.Sess(j).Fc(k).name,...
                 'CallBack',{@cb_menu,'fMRIDesMtx',SPM,j,k},...
                 'UserData',hC,...
                 'HandleVisibility','off')
        end
    end
end

hxvi = uimenu(hExplore, 'Label','Covariance structure', ...
        'Callback',{@cb_menu,'xVi',SPM},...
        'UserData',hC, 'HandleVisibility','off');
if ~isfield(SPM,'xVi') || (isfield(SPM.xVi,'iid') && SPM.xVi.iid) || ...
    ~(isfield(SPM.xVi,'V') || isfield(SPM.xVi,'Vi'))
    set(hxvi,'Enable','off');
end

%-Clear, Quit
%--------------------------------------------------------------------------
uimenu(hC,'Label','Clear','Accelerator','L','Separator','on',...
    'CallBack','spm_results_ui(''Clear'')',...
    'HandleVisibility','off');
uimenu(hC,'Label','Quit','Accelerator','Q','Separator','on',...
    'CallBack','spm_results_ui(''Close'');',...
    'HandleVisibility','off');

%-Pop open 'Interactive' window
%--------------------------------------------------------------------------
figure(Finter)

%-Return handle of menu
%--------------------------------------------------------------------------
varargout = {hC};


%==========================================================================
case 'files&factors'                            %-Summarise files & factors
%==========================================================================
% spm_DesRep('Files&Factors',fnames,I,xC,sF,xs)
if nargin<4, error('insufficient arguments'), end
fnames  = varargin{2};
I       = varargin{3};
xC      = varargin{4};
sF      = varargin{5};
if nargin<6, xs=[]; else xs = varargin{6}; end %-Structure of strings

[fnames,CPath] = spm_str_manip(fnames,'c'); %-extract common path component
nScan          = size(I,1);                 %-#images
bL             = any(diff(I,1),1);          %-Multiple factor levels?
for i=(numel(sF)+1):size(I,2)
    sF{i} = sprintf('sF%d',i);              %-If more factors than expected
end

%-Get Graphics window & window scaling
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)
FS     = spm('FontSizes');

%-Display header information
%--------------------------------------------------------------------------
hTax = axes('Position',[0.03,0.85,0.94,0.1],...
    'DefaultTextFontSize',FS(9),...
    'XLim',[0,1],'YLim',[0,1],...
    'Visible','off');

text(0.5,1,'Statistical analysis: Image files & covariates...',...
    'Fontsize',FS(14),'Fontweight','Bold',...
    'HorizontalAlignment','center');

dx1 = 0.05;
dx2 = 0.08;

x = 0; text(x+.02,.1,'image #','Rotation',90);
for i=numel(bL):-1:1
    if bL(i), x=x+dx1; text(x+.01,.1,sF{i},'Rotation',90); end
end

for j = 1:length(xC)
    n = size(xC(j).rc,2);
    if n>1, tmp=xC(j).cname; else tmp={xC(j).rcname}; end
    for k=1:n
        x=x+dx2;
        text(x,.1,tmp{k},'Rotation',90,'Interpreter','TeX');
    end
end

x = x + dx2;
text(x,0.65,'Base directory:','FontWeight','Bold');
text(x,0.5,CPath,'FontSize',FS(8));
text(x,0.2,'filename tails...');

line('XData',[0 1],'YData',[0 0],'LineWidth',3,'Color','r');

%-Tabulate file & covariate information
%--------------------------------------------------------------------------
hAx = axes('Position',[0.03,0.05,0.94,0.8],...
    'DefaultTextFontSize',FS(8),...
    'Units','points',...
    'Visible','off');
AxPos = get(hAx,'Position'); set(hAx,'YLim',[0,AxPos(4)])

dy = FS(9); y0 = floor(AxPos(4)) -dy; y  = y0;

for i = 1:nScan

    %-Scan indices
    x = 0; text(x,y,sprintf('%03d',i));
    for j=numel(bL):-1:1
        if bL(j), x=x+dx1; text(x,y,sprintf('%02d',I(i,j))); end
    end

    %-Covariates
    for j = 1:length(xC)
        for k=1:size(xC(j).rc,2)
            x=x+dx2;
            text(x,y,sprintf('%6g',xC(j).rc(i,k)),...
                'HorizontalAlignment','Center');
        end
    end

    %-Filename tail(s) - could be multivariate
    x=x+dx2;
    text(x,y,fnames{i});
    y=y-dy;

    %-Paginate if necessary
    if y<dy
        text(0.5,0,sprintf('Page %d',spm_figure('#page')),...
            'FontSize',FS(8),'FontAngle','italic');
        spm_figure('NewPage',[hAx;get(hAx,'Children')])
        hAx = axes('Units','points','Position',AxPos,...
            'DefaultTextFontSize',FS(8),'YLim',[0,AxPos(4)],...
            'Visible','off');
        y = y0;
        text(y,0,'continued...','FontAngle','Italic');
    end
end

line('XData',[0 1],'YData',[y y],'LineWidth',3,'Color','r');


%-Display description strings at bottom of current page
%--------------------------------------------------------------------------
if ~isempty(xs)
    y = y - 2*dy;
    for sf = fieldnames(xs)'
        text(0.3,y,[strrep(sf{1},'_',' '),' :'],...
            'HorizontalAlignment','Right','FontWeight','Bold',...
            'FontSize',FS(9));
        s = xs.(sf{1});
        if ~iscellstr(s), s={s}; end
        for i=1:numel(s)
            text(0.31,y,s{i},'FontSize',FS(9));
            y = y - dy;
        end
    end
end

%-Register last page if paginated
if spm_figure('#page')>1
    text(0.5,0,sprintf('Page %d/%d',spm_figure('#page')*[1,1]),...
        'FontSize',FS(8),'FontAngle','italic');
    spm_figure('NewPage',[hAx;get(hAx,'Children')]);
end

%-Pop up the Graphics window
%--------------------------------------------------------------------------
figure(Fgraph)


%==========================================================================
case 'xvi'
%==========================================================================
% spm_DesRep('xVi',xVi)
if nargin<2, error('insufficient arguments'), end
if ~isstruct(varargin{2}), error('covariance matrix structure required'), end

%-Display
%==========================================================================

%-Get graphics window & FontSizes
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)
FS     = spm('FontSizes');

%-Title
%--------------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
      'DefaultTextFontSize',FS(9),...
      'XLim',[0,1],'YLim',[0,1],...
      'Visible','off');

str  = 'Statistical analysis: Covariance structure';
text(0.5,0.95,str,'Fontsize',FS(14),'Fontweight','Bold',...
      'HorizontalAlignment','center');

line('Parent',hTax,...
      'XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r');


%-Display covariance matrix
%--------------------------------------------------------------------------
hCovMtx(1) = axes('Position',[.07 .4 .6 .4]);
hCovMtxSc  = [];
if isfield(varargin{2},'V')
  clim = full([-max(varargin{2}.V(:))/2 max(varargin{2}.V(:))]); % scale 0 to gray
  hCovMtxIm(1) = imagesc(varargin{2}.V, clim);
  hCovMtxSc = colorbar;
  set(hCovMtxSc,'Ylim',[0 clim(2)]); % cut colorbar at 0
  %-Setup callbacks to allow interrogation of covariance matrix
  %------------------------------------------------------------------------
  set(hCovMtxIm(1),'UserData',...
      varargin{2})
  set(hCovMtxIm(1),'ButtonDownFcn','spm_DesRep(''SurfxVi_CB'')')
  if isfield(varargin{2},'form')
    str = sprintf('Covariance structure V: %s',varargin{2}.form);
  else
    str = 'Covariance structure V';
  end
  if isfield(varargin{2},'var') && isfield(varargin{2},'dep') && ...
     isfield(varargin{2},'sF') && isfield(varargin{2},'I')
    r = find((max(varargin{2}.I) > 1) & ~varargin{2}.var & ~varargin{2}.dep);
    if any(varargin{2}.dep)
      cmstr = 'yes';
    else
      cmstr = 'no';
    end
    str = sprintf(['Covariance structure V - replications over ''%s''\n'...
        'Correlated repeated measures: %s'], ...
      varargin{2}.sF{r},cmstr);
  end
  xlabel(str);
else
  text(.5,.5, 'Covariance not (yet) estimated.', 'HorizontalAlignment','center');
end

if isfield(varargin{2},'h')
  hPEstAx   = axes('Position',[.07 .315 .6 .025],...
      'DefaultTextInterpreter','TeX');
  hParEstIm   = imagesc(varargin{2}.h',clim);
  set(hPEstAx,...
      'XLim',[0,length(varargin{2}.h)]+.5,'XTick',[1:length(varargin{2}.h)-1]+.5,...
      'XTickLabel','',...
      'YLim',[0,1]+.5,'YDir','reverse','YTick',[],...
      'Box','on','TickDir','in','XGrid','on','GridLineStyle','-');
  xlabel('hyperparameter estimates')
  set(hParEstIm,'UserData',varargin{2}.h)
  set(hParEstIm,'ButtonDownFcn','spm_DesRep(''SurfHPEstIm_CB'')')
else
  hPEstAx = [];
end
spm_figure('NewPage',[hCovMtx;get(hCovMtx,'Children');hPEstAx;get(hPEstAx,'Children');hCovMtxSc])
 
%-Show components of covariance matrix
%--------------------------------------------------------------------------
if isfield(varargin{2},'Vi')
    for k = 1:length(varargin{2}.Vi)
        %-Display covariance component xVi.Vi{k}
        %------------------------------------------------------------------
        hCovMtx(k+1) = axes('Position',[.07 .4 .6 .4]);
        hCovMtxIm(k+1) = imagesc(varargin{2}.Vi{k});
        xlabel(sprintf('Covariance component Vi{%d}', k));
        if k<length(varargin{2}.Vi)
            spm_figure('NewPage',[hCovMtx(k+1);get(hCovMtx(k+1),'Children')])
        end
    end
    if spm_figure('#page')>1
        spm_figure('NewPage',[hCovMtx(k+1);get(hCovMtx(k+1),'Children')])
    end
end


%==========================================================================
case {'desmtx','desorth'}    %-Display design matrix / design orthogonality
%==========================================================================
% spm_DesRep('DesMtx',xX,fnames,xs)
% spm_DesRep('DesOrth',xX)
if nargin<2, error('insufficient arguments'), end
if ~isstruct(varargin{2}), error('design matrix structure required'), end
if nargin<3, fnames={}; else fnames=varargin{3}; end
if nargin<4, xs=[]; else xs=varargin{4}; end

desmtx = strcmpi(varargin{1},'desmtx');

%-Locate DesMtx (X), scaled DesMtx (nX) & get parameter names (Xnames)
%--------------------------------------------------------------------------
if isfield(varargin{2},'xKXs') && ...
        ~isempty(varargin{2}.xKXs) && isstruct(varargin{2}.xKXs)
    iX = 1;
    [nScan,nPar] = size(varargin{2}.xKXs.X);
elseif isfield(varargin{2},'X') && ~isempty(varargin{2}.X)
    iX = 0;
    [nScan,nPar] = size(varargin{2}.X);
else
    error('Can''t find DesMtx in this structure!')
end

if isfield(varargin{2},'nKX') && ~isempty(varargin{2}.nKX)
    inX = 1; else inX = 0; end

if isfield(varargin{2},'name') && ~isempty(varargin{2}.name)
    Xnames = varargin{2}.name; else Xnames = {}; end


%-Compute design orthogonality matrix if DesOrth
%--------------------------------------------------------------------------
if ~desmtx
    if iX
        tmp = sqrt(sum(varargin{2}.xKXs.X.^2));
        O   = varargin{2}.xKXs.X'*varargin{2}.xKXs.X./kron(tmp',tmp);
        tmp = sum(varargin{2}.xKXs.X);
    else
        tmp = sqrt(sum(varargin{2}.X.^2));
        O   = varargin{2}.X'*varargin{2}.X./kron(tmp',tmp);
        tmp = sum(varargin{2}.X);
    end
    tmp     = abs(tmp)<eps*1e5;
    bC      = kron(tmp',tmp);
end


%-Display
%==========================================================================

%-Get graphics window & FontSizes
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)
FS = spm('FontSizes');

%-Title
%--------------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
    'DefaultTextFontSize',FS(9),...
    'XLim',[0,1],'YLim',[0,1],...
    'Visible','off');

str='Statistical analysis: Design'; if ~desmtx, str=[str,' orthogonality']; end
text(0.5,0.95,str,'Fontsize',FS(14),'Fontweight','Bold',...
    'HorizontalAlignment','center');

line('Parent',hTax,...
    'XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r');


%-Display design matrix
%--------------------------------------------------------------------------
hDesMtx = axes('Position',[.07 .4 .6 .4]);
if inX      %-Got a scaled DesMtx
    hDesMtxIm = image((varargin{2}.nKX + 1)*32);
elseif iX   %-No scaled DesMtx, DesMtx in .xKXs structure
    hDesMtxIm = image((spm_DesMtx('sca',varargin{2}.xKXs.X,Xnames) + 1)*32);
else        %-No scaled DesMtx, no .xKXs, DesMtx in .X
    hDesMtxIm = image((spm_DesMtx('sca',varargin{2}.X,     Xnames) + 1)*32);
end

STick = spm_DesRep('ScanTick',nScan,32);
PTick = spm_DesRep('ScanTick',nPar,32);

set(hDesMtx,'TickDir','out',...
    'XTick',PTick,'XTickLabel','',...
    'YTick',STick,'YTickLabel','')
if desmtx
    xlabel('parameters'), ylabel('images')
else
    ylabel('design matrix')
end

%-Parameter names
if ~isempty(Xnames)
    axes('Position',[.07 .8 .6 .1],'Visible','off',...
        'DefaultTextFontSize',FS(8),'DefaultTextInterpreter','TeX',...
        'XLim',[0,nPar]+0.5)
    for i=PTick, text(i,.05,Xnames{i},'Rotation',90), end
end

%-Filenames
% ( Show at most 32, showing every 2nd/3rd/4th/... as necessary to pair )
% ( down to <32 items. Always show last item so #images is indicated.   )     
if desmtx && ~isempty(fnames)
    axes('Position',[.68 .4 .3 .4],'Visible','off',...
        'DefaultTextFontSize',FS(8),...
        'YLim',[0,nScan]+0.5,'YDir','Reverse')
    for i = STick
        try
            str  = strrep(fnames(i,:),'\','/');
        catch
            str  = strrep(fnames{i},'\','/');
        end
        text(0,i,spm_file(str,'short35'));
    end
end

%-Setup callbacks to allow interrogation of design matrix
%--------------------------------------------------------------------------
if iX,  set(hDesMtxIm,'UserData',...
    struct('X',varargin{2}.xKXs.X,'Xnames',{Xnames},'fnames',{fnames}))
else   set(hDesMtxIm,'UserData',...
    struct('X',varargin{2}.X,     'Xnames',{Xnames},'fnames',{fnames}))
end
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')

if desmtx
    %-Parameter estimability/uniqueness
    %----------------------------------------------------------------------
    hPEstAx     = axes('Position',[.07 .315 .6 .025],...
            'DefaultTextInterpreter','TeX');
    if iX,  est = spm_SpUtil('IsCon',varargin{2}.xKXs);
    else    est = spm_SpUtil('IsCon',varargin{2}.X); end
    hParEstIm = image((est+1)*32);
    set(hPEstAx,...
        'XLim',[0,nPar]+.5,'XTick',[1:nPar-1]+.5,'XTickLabel','',...
        'YLim',[0,1]+.5,'YDir','reverse','YTick',[],...
        'Box','on','TickDir','in','XGrid','on','GridLineStyle','-');
    xlabel('parameter estimability')
    text((nPar+0.5 + nPar/30),1,...
        '(gray \rightarrow \beta not uniquely specified)',...
        'Interpreter','TeX','FontSize',FS(8))
    set(hParEstIm,'UserData',struct('est',est,'Xnames',{Xnames}))
    set(hParEstIm,'ButtonDownFcn','spm_DesRep(''SurfEstIm_CB'')')
else
    %-Design orthogonality
    %----------------------------------------------------------------------
    hDesO   = axes('Position',[.07 .18 .6 .2]);
    tmp = 1-abs(O); tmp(logical(tril(ones(nPar),-1))) = 1;
    hDesOIm = image(tmp*64);
    
    set(hDesO,'Box','off','TickDir','out',...
        'XaxisLocation','top','XTick',PTick,'XTickLabel','',...
        'YaxisLocation','right','YTick',PTick,'YTickLabel','',...
        'YDir','reverse')
    tmp = [1,1]'*[[0:nPar]+0.5];
    line('Xdata',tmp(1:end-1)','Ydata',tmp(2:end)')

    xlabel('design orthogonality')
    set(get(hDesO,'Xlabel'),'Position',[0.5,nPar,0],...
        'HorizontalAlignment','left',...
        'VerticalAlignment','top')
    set(hDesOIm,...
        'UserData',struct('O',O,'bC',bC,'Xnames',{Xnames}),...
        'ButtonDownFcn','spm_DesRep(''SurfDesO_CB'')')

    if ~isempty(Xnames)
        axes('Position',[.69 .18 0.01 .2],'Visible','off',...
            'DefaultTextFontSize',FS(10),...
            'DefaultTextInterpreter','TeX',...
            'YDir','reverse','YLim',[0,nPar]+0.5)
        for i=PTick
            text(0,i,Xnames{i},'HorizontalAlignment','left')
        end
    end

end

%-Design descriptions
%--------------------------------------------------------------------------
if desmtx
    str = 'Design description...';
    line('Parent',hTax,...
        'XData',[0.3 0.7],'YData',[0.28 0.28],'LineWidth',3,'Color','r')
    hAx = axes('Position',[0.03,0.05,0.94,0.22],'Visible','off');
else
    str = '';
    line('Parent',hTax,...
        'XData',[0.3 0.7],'YData',[0.14 0.14],'LineWidth',3,'Color','r')
    hAx = axes('Position',[0.03,0.05,0.94,0.08],'Visible','off');
    xs = struct('Measure',  ['abs. value of cosine of angle between ',...
                 'columns of design matrix'],...
            'Scale',    {{  'black - colinear (cos=+1/-1)';...
                    'white - orthogonal (cos=0)';...
                    'gray  - not orthogonal or colinear'}});
end

if ~isempty(xs)
    set(hAx,'Units','points');
    AxPos = get(hAx,'Position');
    set(hAx,'YLim',[0,AxPos(4)])
    
    dy = FS(9); y0 = floor(AxPos(4)) -dy; y = y0;

    text(0.3,y,str,...
        'HorizontalAlignment','Center',...
        'FontWeight','Bold','FontSize',FS(11))
    y=y-2*dy;
    
    for sf = fieldnames(xs)'
        text(0.3,y,[strrep(sf{1},'_',' '),' :'],...
            'HorizontalAlignment','Right','FontWeight','Bold',...
            'FontSize',FS(9))
        s = xs.(sf{1});
        if ~iscellstr(s), s={s}; end
        for i=1:numel(s)
            text(0.31,y,s{i},'FontSize',FS(9))
            y=y-dy;
        end
    end
end

%-Pop up the Graphics window
%--------------------------------------------------------------------------
figure(Fgraph)


%==========================================================================
case 'fmridesmtx'                %-Interactive review of fMRI design matrix
%==========================================================================
%spm_DesRep('fMRIDesMtx',SPM,s,i)
SPM  = varargin{2};
Sess = SPM.Sess;
if nargin < 4, i = 1; else i = varargin{4}; end
if nargin < 3, s = 1; else s = varargin{3}; end

%-Get Graphics window
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)


%-Trial-specific regressors - time domain
%--------------------------------------------------------------------------
sX    = SPM.xX.X(Sess(s).row,Sess(s).col);
rX    = sX(:,Sess(s).Fc(i).i);
subplot(2,2,1)
plot(Sess(s).row,rX)
xlabel('scan')
ylabel('regressor[s]')
title({'Time domain',['Regressors for ' Sess(s).Fc(i).name]})
grid on
axis tight

%-Trial-specific regressors - frequency domain
%--------------------------------------------------------------------------
subplot(2,2,2)
gX    = abs(fft(rX)).^2;
gX    = gX*diag(1./sum(gX));
q     = size(gX,1);
Hz    = [0:(q - 1)]/(q*SPM.xY.RT);
q     = 2:fix(q/2);
plot(Hz(q),gX(q,:))
HPF   = SPM.xX.K(s).HParam;
patch([0 1 1 0]/HPF,[0 0 1 1]*max(max(gX)),[1 1 1]*.9,'facealpha',.5);
xlabel('Frequency (Hz)')
ylabel('relative spectral density')
h=title(['Frequency domain',sprintf('\n'), ' {\bf',num2str(HPF),'}', ...
    ' second High-pass filter'],'Interpreter','Tex');
grid on
axis tight

% if trial (as opposed to trial x trial interaction)
%--------------------------------------------------------------------------
if length(Sess(s).U) >= i

    % Basis set and peristimulus sampling
    %----------------------------------------------------------------------
    subplot(2,2,3)
    dt   = Sess(s).U(i).dt;
    RT   = SPM.xY.RT;
    t    = [1:size(SPM.xBF.bf,1)]*dt;
    pst  = Sess(s).U(i).pst;
    plot(t,SPM.xBF.bf,pst,0*pst,'.','MarkerSize',16)
    str  = sprintf('TR = %0.2fs',RT);
    xlabel({'time {secs}' str sprintf('%0.0fms time bins',1000*dt)})
    title({'Basis set and peristimulus sampling' SPM.xBF.name})
    axis tight
    grid on

    % if a paramteric variate is specified
    %----------------------------------------------------------------------
    for p = 1:length(Sess(s).U(i).P)

        if Sess(s).U(i).P(p).h

        % onsets and parametric modulation
        %------------------------------------------------------------------
        subplot(2,2,4)
        ons = Sess(s).U(i).ons;
        plot(ons,Sess(s).U(i).P(p).P,'.','MarkerSize',8)
        xlabel('time {secs}')
        title('Parameters')
        grid on
        hold on

        end
    end
end

%-Pop up Graphics figure window
%--------------------------------------------------------------------------
figure(Fgraph);


%==========================================================================
case 'covs'                   %-Plot and describe covariates (one per page)
%==========================================================================
% spm_DesRep('Covs',xX,xC)
if nargin<3, error('insufficient arguments'), end
xC     = varargin{3};   %-Struct array of covariate information
if ~isstruct(varargin{2}), error('design matrix structure required'), end

if isempty(xC), spm('alert!','No covariates!',mfilename), return, end

%-Get graphics window & window scaling
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph,0)
FS     = spm('FontSizes');

%-Title
%--------------------------------------------------------------------------
hTax = axes('Position',[0.03,0,0.94,1],...
    'DefaultTextFontSize',FS(9),...
    'XLim',[0,1],'YLim',[0,1],...
    'Visible','off');

text(0.5,0.95,'Statistical analysis: Covariates',...
    'Fontsize',FS(14),'Fontweight','Bold',...
    'HorizontalAlignment','center')

text(0.5,0.82,'(covariates plotted over transposed design matrix)',...
    'FontSize',FS(8),'HorizontalAlignment','center')

line('XData',[0.3 0.7],'YData',[0.92 0.92],'LineWidth',3,'Color','r')
line('XData',[0.3 0.7],'YData',[0.44 0.44],'LineWidth',3,'Color','r')


%-Design matrix (as underlay for plots) and parameter names
%--------------------------------------------------------------------------
[nScan,nPar] = size(varargin{2}.X);
if isfield(varargin{2},'name') && ~isempty(varargin{2}.name)
    Xnames = varargin{2}.name; else Xnames = {}; end

%-Design matrix
hDesMtx = axes('Position',[.1 .5 .7 .3]);
if isfield(varargin{2},'nKX') && ~isempty(varargin{2}.nKX)
    image(varargin{2}.nKX'*32+32)
elseif isfield(varargin{2},'xKXs') && ~isempty(varargin{2}.xKXs)
    image(spm_DesMtx('sca',varargin{2}.xKXs.X,Xnames)*32+32)
else
    image(spm_DesMtx('sca',varargin{2}.X,Xnames)*32+32)
end
set(hDesMtx,'Visible','off')

%-Parameter names
hParAx = axes('Position',[.8 .5 .2 .3],'Visible','off',...
    'DefaultTextFontSize',FS(8),'DefaultTextInterpreter','TeX',...
    'YLim',[0.5,nPar+0.5],'YDir','Reverse');
for i = 1:nPar, hPNames(i) = text(.05,i,Xnames{i}); end


%-Covariates - one page each
%--------------------------------------------------------------------------
for i = 1:length(xC)

    %-Title
    %----------------------------------------------------------------------
    hSTitle = text(0.5,0.87,sprintf('%d : %s',i,xC(i).rcname),...
            'Parent',hTax,...
            'HorizontalAlignment','center',...
            'FontSize',FS(13),'FontWeight','Bold');

    %-Plot
    %----------------------------------------------------------------------
    hAx = axes('Position',[.1 .5 .7 .3],...
            'TickDir','out','Box','off','Color','none',...
            'NextPlot','add',...
            'XLim',[0,nScan]+0.5);
    plot(xC(i).rc,'LineWidth',2)
    if nScan<48, plot(xC(i).rc,'.k','MarkerSize',20); end
    xlabel('image #')
    ylabel('covariate value')


    %-Descriptions
    %----------------------------------------------------------------------
    hDAx = axes('Position',[0.03,0.1,0.94,0.30],'Visible','off');
    
    set(hDAx,'Units','points');
    tmp = get(hDAx,'Position');
    set(hDAx,'YLim',[0,tmp(4)])
    
    dy = FS(9); y0 = floor(tmp(4)) -dy; y = y0;

    %-Description strings from xC(i).descrip
    text(0.3,y,'Details :',...
        'HorizontalAlignment','Right',...
        'FontWeight','Bold','FontSize',FS(9))
    s = xC(i).descrip;
    if ~iscellstr(s), s={s}; end
    for j=1:numel(s)
        text(0.31,y,s{j},'FontSize',FS(9))
        y=y-dy;
    end
    y=y-dy;

    %-Key (if block of covariates entered)
    %----------------------------------------------------------------------
    if size(xC(i).rc,2)>1
        ColorOrder = get(hAx,'ColorOrder');
        text(0.3,y,'Key :',...
            'HorizontalAlignment','Right',...
            'FontWeight','Bold','FontSize',FS(9))
        for j = 1:size(xC(i).rc,2)
            color = ColorOrder(mod(j-1,size(ColorOrder,1))+1,:);
            if size(xC(i).rc,2)==length(xC(i).cname)
                str = xC(i).cname{j};
            else
                str = sprintf('column %d',j);
            end
            text(0.31,y,str,'FontSize',FS(9),...
                'Color',color)
            text(0.5,xC(i).rc(1,j),[str,' \rightarrow'],...
                'Parent',hAx,...
                'FontSize',FS(8),'FontWeight','Bold',...
                'HorizontalAlignment','Right',...
                'Interpreter','TeX',...
                'Color',color)
            y=y-dy;
        end
        y=y-dy;
    end


    %-Associated parameters
    %----------------------------------------------------------------------
    text(0.3,y,'Design matrix columns :',...
        'HorizontalAlignment','Right',...
        'FontWeight','Bold','FontSize',FS(9))
    if isempty(xC(i).cols)
        text(0.31,y,'(none)','FontSize',FS(9))
    else
        for j = xC(i).cols
            text(0.31,y,sprintf('%d : %s',j,Xnames{j}),...
                'FontSize',FS(9),'Interpreter','TeX')
            y=y-dy;
        end
    end
    y=y-dy;


    %-Highlight parameter names
    %----------------------------------------------------------------------
    hCurPNames = hPNames(xC(i).cols);
    set(hCurPNames,'Color','r','FontWeight','Bold','FontSize',FS(8))


    %-Paginate (if more than one covariate)
    %----------------------------------------------------------------------
    if length(xC)>1
        spm_figure('NewPage',[hSTitle; hAx; get(hAx,'Children');...
            hCurPNames; hDAx; get(hDAx,'Children')]);
    end

end

%-Pop up the Graphics window
%--------------------------------------------------------------------------
figure(Fgraph)


%==========================================================================
case 'scantick'
%==========================================================================
% spm_DesRep('ScanTick',nScan,lim)
% ( Show at most 32, showing every 2nd/3rd/4th/... as necessary to pair )
% ( down to <32 items. Always show last item so #images is indicated.   )     
if nargin<3, lim=32; else lim=varargin{3}; end
if nargin<2, error('insufficient arguments'), end
nScan = varargin{2};

p = max(1,ceil(nScan/lim));
s = 1:p:nScan; s(end)=nScan;

varargout = {s,lim};


%==========================================================================
case {'surfdesmtx_cb','surfdesmtxmo_cb','surfdesmtxup_cb'}    %-Surf DesMtx
%==========================================================================
% spm_DesRep('SurfDesMtx_CB')
% spm_DesRep('SurfDesMtxMo_CB')
% spm_DesRep('SurfDesMtxUp_CB')

h = get(gca,'Xlabel');

if strcmpi(varargin{1},'surfdesmtxup_cb')
    UD = get(h,'UserData');
    set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
        'UserData',UD.UserData)
    set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
    return
end

if strcmpi(varargin{1},'surfdesmtx_cb')
    UD = struct('String',      get(h,'String'),...
                'Interpreter', get(h,'Interpreter'),...
                'UserData',    get(h,'UserData'));
    set(h,'UserData',UD)
    set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfDesMtxMo_CB'')',...
         'WindowButtonUpFcn','spm_DesRep(''SurfDesMtxUp_CB'')')
end

mm  = [get(gca,'YLim')',get(gca,'XLim')']+[.5,.5;-.5,-.5];
ij  = get(gca,'CurrentPoint');
ij  = round(min(max(ij(1,[2,1]),mm(1,:)),mm(2,:)));

istr = 'none';
switch get(gcbf,'SelectionType')
    case 'normal'
        try, str = sprintf('X(%d,%d) = %g',ij(1),ij(2),...
                subsref(get(gco,'UserData'),...
                struct('type',{'.','()'},'subs',{'X',{ij(1),ij(2)}})));
        catch, str='(no cached design matrix to surf)'; end
    case 'extend'
        try, str = sprintf('Image %d: %s',ij(1),...
                char(spm_file(...
                subsref(get(gco,'UserData'),...
                struct('type',{'.','()'},...
                'subs',{'fnames',{ij(1),':'}})),'short40')));
        catch, str='(no cached image filenames to surf)'; end
    case 'alt'
        try, str = sprintf('Parameter %d: %s',ij(2),...
                subsref(get(gco,'UserData'),...
                struct('type',{'.','{}'},'subs',{'Xnames',{ij(2)}})));
            istr = 'tex';
        catch, str='(no cached parameter names to surf)'; end
    case 'open'
        try,    assignin('base','ans',subsref(get(gco,'UserData'),...
                struct('type',{'.'},'subs',{'X'})))
            evalin('base','ans')
        catch,  fprintf('%s GUI: can''t find design matrix\n',mfilename)
        end
        return
end

set(h,'String',str,'Interpreter',istr)


%==========================================================================
case {'surfestim_cb','surfestimmo_cb','surfestimup_cb'}     %-Surf ParEstIm
%==========================================================================
% spm_DesRep('SurfEstIm_CB')
% spm_DesRep('SurfEstImMo_CB')
% spm_DesRep('SurfEstImUp_CB')

h = get(gca,'Xlabel');

if strcmpi(varargin{1},'surfestimup_cb')
    UD = get(h,'UserData');
    set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
        'UserData',UD.UserData)
    set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
    return
end

if strcmpi(varargin{1},'surfestim_cb')
    UD = struct(    'String',   get(h,'String'),...
            'Interpreter',  get(h,'Interpreter'),...
            'UserData', get(h,'UserData'));
    set(h,'UserData',UD)
    set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfEstImMo_CB'')',...
         'WindowButtonUpFcn',    'spm_DesRep(''SurfEstImUp_CB'')')
end

mm  = [get(gca,'XLim')]+[.5,-.5];
i   = get(gca,'CurrentPoint');
i   = round(min(max(i(1,1),mm(1)),mm(2)));

istr = 'none';
switch get(gcbf,'SelectionType')
case 'normal'
    try, tmp = {' (not unique)',' (unique)'};
    str = sprintf('Parameter %d : %s%s',...
        i,...
        subsref(get(gco,'UserData'),...
            struct('type',{'.','{}'},'subs',{'Xnames',{i}})),...
        tmp{subsref(get(gco,'UserData'),...
            struct('type',{'.','()'},'subs',{'est',{i}}))+1});
        istr = 'tex';
    catch, str='(no cached data to surf)'; end
case {'extend','alt'}
    return
case 'open'
    try,  UD = get(gco,'UserData');
        assignin('base','ans',...
            subsref(get(gco,'UserData'),...
                struct('type',{'.'},'subs',{'est'})))
        evalin('base','ans')
    catch,  fprintf('%s GUI: can''t find design orthogonality\n',mfilename)
    end
    return
end

set(h,'String',str,'Interpreter',istr);


%==========================================================================
case {'surfdeso_cb','surfdesomo_cb','surfdesoup_cb'}       %-Surf DesOrthIm
%==========================================================================
% spm_DesRep('SurfDesO_CB')
% spm_DesRep('SurfDesOMo_CB')
% spm_DesRep('SurfDesOUp_CB')

h = get(gca,'Xlabel');

if strcmpi(varargin{1},'surfdesoup_cb')
    UD = get(h,'UserData');
    set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
        'UserData',UD.UserData)
    set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
    return
end

if strcmpi(varargin{1},'surfdeso_cb')
    UD = struct(    'String',   get(h,'String'),...
            'Interpreter',  get(h,'Interpreter'),...
            'UserData', get(h,'UserData'));
    set(h,'UserData',UD)
    set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfDesOMo_CB'')',...
         'WindowButtonUpFcn',    'spm_DesRep(''SurfDesOUp_CB'')')
end

mm  = [get(gca,'YLim')',get(gca,'XLim')']+[.5,.5;-.5,-.5];
ij  = get(gca,'CurrentPoint');
ij  = round(min(max(ij(1,[2,1]),mm(1,:)),mm(2,:)));
if ij(1)>ij(2), return, end

istr = 'none';
switch get(gcbf,'SelectionType')
case 'normal'
    try
        UD = get(gco,'UserData');
        if abs(abs(UD.O(ij(1),ij(2)))-1) < eps*1e1
            str = '{\bf colinear}';
        elseif abs(UD.O(ij(1),ij(2))) < eps*1e1
            str = '{\bf orthogonal}';
        else
            str = '{\bf not orthogonal}';
        end
        if ~diff(ij), str=[str,' {\it(same column)}']; end
        if UD.bC(ij(1),ij(2)), tmp=' ={\it r}'; else tmp=''; end
        str = { sprintf('{\\bf %s} (col %d)',...
                UD.Xnames{ij(1)},ij(1)),...
                sprintf('& {\\bf %s} (col %d)',...
                UD.Xnames{ij(2)},ij(2)),...
                sprintf('cos(\\theta)%s = %1.2f',...
                tmp,UD.O(ij(1),ij(2))),...
            ['\rightarrow ',str]};
        istr = 'tex';
    catch, str='(no cached data to surf)'; end
case {'extend','alt'}
    return
case 'open'
    try,    UD = get(gco,'UserData');
        assignin('base','ans',UD.O)
        evalin('base','ans')
    catch,  fprintf('%s GUI: can''t find design orthogonality\n',mfilename)
    end
    return
end

set(h,'String',str,'Interpreter',istr)


%==========================================================================
case {'surfcon_cb','surfconmo_cb','surfconup_cb'}           %-Surf Contrast
%==========================================================================
% spm_DesRep('SurfCon_CB')
% spm_DesRep('SurfConOMo_CB')
% spm_DesRep('SurfConOUp_CB')

cUD = get(gco,'UserData');
if ~isstruct(cUD) || ~isfield(cUD,'h')
    warning('contrast GUI objects setup incorrectly'), return
end
h    = cUD.h;

if strcmpi(varargin{1},'surfconup_cb')
    UD = get(h,'UserData');
    set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
        'UserData',UD.UserData)
    set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
    return
end

if strcmpi(varargin{1},'surfcon_cb')
    UD = struct(    'String',   get(h,'String'),...
            'Interpreter',  get(h,'Interpreter'),...
            'UserData', get(h,'UserData'));
    set(h,'UserData',UD)
    set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfConMo_CB'')',...
         'WindowButtonUpFcn',    'spm_DesRep(''SurfConUp_CB'')')
end

mm  = [get(gca,'YLim')',get(gca,'XLim')']+[.5,.5;-.5,-.5];
ij  = get(gca,'CurrentPoint');
ij  = round(min(max(ij(1,[2,1]),mm(1,:)),mm(2,:)));

istr = 'none';
switch get(gcbf,'SelectionType')
case 'normal'
    try
        if cUD.i>0, str = sprintf('%d',cUD.i); else str = ''; end
        switch get(gco,'Type')
        case 'image'
            str = sprintf('%s\\{F\\}: {\\bf%s} (%d,%d) = %.2f',...
                str,cUD.xCon.name,ij(2),ij(1),...
                cUD.xCon.c(ij(2),ij(1)));
        case 'patch'
            str = sprintf('%s\\{T\\}: {\\bf%s} (%d) = %.2f',...
                str,cUD.xCon.name,ij(2),...
                cUD.xCon.c(ij(2)));
        otherwise, error('unexpected object type')
        end
        istr = 'TeX';
    catch, str='(no cached data to surf)'; end
case {'alt','extend'}
    return
case 'open'
    try,    assignin('base','ans',cUD.xCon.c')
        evalin('base','ans')
    catch,  fprintf('%s GUI: can''t find contrast\n',mfilename)
    end
    return
end

set(h,'String',str,'Interpreter',istr);


%==========================================================================
case {'surfxvi_cb','surfxvimo_cb','surfxviup_cb'}                %-Surf Xvi
%==========================================================================
% spm_DesRep('SurfxVi_CB')
% spm_DesRep('SurfxViMo_CB')
% spm_DesRep('SurfxViUp_CB')

h = get(gca,'Xlabel');

if strcmpi(varargin{1},'surfxviup_cb')
      UD = get(h,'UserData');
      set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
            'UserData',UD.UserData)
      set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
      return
end


if strcmpi(varargin{1},'surfxvi_cb')
      UD = struct(    'String',       get(h,'String'),...
                      'Interpreter',  get(h,'Interpreter'),...
                      'UserData',     get(h,'UserData'));
      set(h,'UserData',UD)
      set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfxViMo_CB'')',...
               'WindowButtonUpFcn',    'spm_DesRep(''SurfxViUp_CB'')')
end

mm  = [get(gca,'YLim')',get(gca,'XLim')']+[.5,.5;-.5,-.5];
ij  = get(gca,'CurrentPoint');
ij  = round(min(max(ij(1,[2,1]),mm(1,:)),mm(2,:)));
istr = 'none';

switch get(gcbf,'SelectionType')
case 'normal'
    try
        str = sprintf('V(%d,%d) = %g',ij(1),ij(2),...
              full(subsref(get(gco,'UserData'),...
              struct('type',{'.','()'},'subs',{'V',{ij(1),ij(2)}}))));
    catch
        str = '(no cached covariance matrix to surf)';
    end
case 'extend'
    try
        ind = 1:length(subsref(get(gco,'Userdata'),...
            struct('type','.','subs','h')));
        isel = false(size(ind));
        for k = 1:length(ind)
            isel(k) = subsref(get(gco,'UserData'),...
                struct('type',{'.','{}','()'},'subs',{'Vi',{k},{ij(1),ij(2)}})) ~= 0;
        end
        if any(isel)
            str = [sprintf('V(%d,%d): ',ij(1),ij(2)) sprintf('Vi{%d}',ind(isel))];
        else
            str = sprintf('no Vi at (%d,%d)',ij(1),ij(2));
        end
    catch
        return
    end
case 'alt'
    return
case 'open'
    try
        assignin('base','ans',subsref(get(gco,'UserData'),...
            struct('type',{'.'},'subs',{'V'})));
        evalin('base','ans');
    catch
        fprintf('%s GUI: can''t find covariance matrix\n',mfilename)
    end
    return
end

set(h,'String',str,'Interpreter',istr)


%==========================================================================
case {'surfhpestim_cb','surfhpestimmo_cb','surfhpestimup_cb'}  %-Surf ParHpestim
%==========================================================================
% spm_DesRep('SurfHPEstim_CB')
% spm_DesRep('SurfHPEstimMo_CB')
% spm_DesRep('SurfHPEstimUp_CB')

h = get(gca,'Xlabel');

if strcmpi(varargin{1},'surfhpestimup_cb')
      UD = get(h,'UserData');
      set(h,'String',UD.String,'Interpreter',UD.Interpreter,...
              'UserData',UD.UserData)
      set(gcbf,'WindowButtonMotionFcn','','WindowButtonUpFcn','')
      return
end

if strcmpi(varargin{1},'surfhpestim_cb')
      UD = struct(    'String',       get(h,'String'),...
                      'Interpreter',  get(h,'Interpreter'),...
                      'UserData',     get(h,'UserData'));
      set(h,'UserData',UD)
      set(gcbf,'WindowButtonMotionFcn','spm_DesRep(''SurfHPEstimMo_CB'')',...
               'WindowButtonUpFcn',    'spm_DesRep(''SurfHPEstimUp_CB'')')
end

mm  = [get(gca,'XLim')]+[.5,-.5];
i   = get(gca,'CurrentPoint');
i   = round(min(max(i(1,1),mm(1)),mm(2)));

istr = 'none';
switch get(gcbf,'SelectionType')
case 'normal'
      try,
      str = sprintf('Hyperparameter %d : %f',...
              i,...
              subsref(get(gco,'UserData'),...
                      struct('type',{'()'},'subs',{{i}})));
              istr = 'tex';
      catch, str='(no cached data to surf)'; end
case {'extend','alt'}
      return
case 'open'
      try,    UD = get(gco,'UserData');
              assignin('base','ans',...
                      get(gco,'UserData'));
              evalin('base','ans')
      catch,  fprintf('%s GUI: can''t find hyperparameter estimates\n',mfilename)
      end
      return
end

set(h,'String',str,'Interpreter',istr);


%==========================================================================
otherwise                                           %-Unknown action string
%==========================================================================
error(['Unknown action string: ',varargin{1}])


%==========================================================================
end


%==========================================================================
% function cb_menu(obj,evt,action,SPM,varargin)
%==========================================================================
function cb_menu(obj,evt,action,SPM,varargin)

switch action
    case {'DesMtx','Files&Factors'}
        try
            filenames = reshape(cellstr(SPM.xY.P),size(SPM.xY.VY));
        catch
            filenames = {};
        end
end

switch action
    case 'DesMtx'
        spm_DesRep('DesMtx',SPM.xX,...
            filenames,...
            SPM.xsDes);
            
    case 'DesOrth'
        spm_DesRep('DesOrth',SPM.xX);

    case 'Files&Factors'
        spm_DesRep('Files&Factors',...
            filenames,...
            SPM.xX.I,SPM.xC,SPM.xX.sF,SPM.xsDes);

    case 'Covs'
        spm_DesRep('Covs',SPM.xX,SPM.xC);

    case 'fMRIDesMtx'
        spm_DesRep('fMRIDesMtx',SPM,varargin{1},varargin{2});

    case 'xVi'
        spm_DesRep('xVi', SPM.xVi);
end
