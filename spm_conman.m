function varargout=spm_conman(varargin)
% Contrast manager: GUI for defining valid contrasts
% FORMAT varargout=spm_conman(varargin)
%       - An embedded callback, multi-function function
%       - For detailed programmers comments,
%         see format specifications in main body of program (below user help)
%__________________________________________________________________________
%
% The contrast manager is a user interface for the selection and definition
% of contrasts for a given SPM design. At present, the contrast manager
% provides only interactive GUI facilities.
%
% See also: spm_getSPM.m, spm_contrasts.m
%
%==========================================================================
% U s i n g   t h e   S P M   C o n t r a s t   M a n a g e r   G U I
%==========================================================================
%
% The contrast manager graphicsl user interface (GUI) is a dialog box for
% the specification and selection of contrasts. The contrast selection 
% interface is presented initially, pressing the "Define new contrast..." 
% button pops up the contrast definition interface:
%
%__________________________________________________________________________
%                                      ConMan: Contrast selection interface
%
% The contrast selection interface consists of:
%
% * "Show" statistic type radio-buttons: "t-contrasts", "F-contrasts", "All":
%   Selects the types of contrast that are presented for selection in the
%   contrast selection list box. Depending on the use of the contrast
%   manager, some options may be unavailable.
%
% * List-box of currently defined contrasts.
%   Each contrast is listed by number, type ("T" or "F"), and name. Only
%   contrasts of types specified by the settings of the "show"
%   radiobuttons are shown.
%   Select contrasts by clicking on their entry in the list-box. To
%   select multiple contrasts (for conjunctions or masking, if
%   appropriate), the standard techniques of drag selection and
%   shift-clicking can be used for selecting a set of adjacent contrasts,
%   control-click to select individual contrasts separated in the list.
%   Selected contrasts are highlit in black.
%
% * Image of the design matrix:
%   A grey-scale representation of the design matrix, to aid interpretation
%   of the contrasts.
%
%   The design matrix is "surfable": Clicking (and holding or dragging)
%   around the design matrix image reports the corresponding value of the
%   design matrix ('normal' click - "left" mouse button usually), the
%   image filename ('extend' mouse click - "middle" mouse), or parameter
%   name ('alt' click - "right" mouse). Double clicking the design matrix
%   image extracts the design matrix into the base MATLAB workspace.
%   (Surfing is handled by spm_DesRep.m)
%
% * Parameter estimability bar
%   Under the design matrix the parameter estimability is displayed as
%   a 1xp matrix of grey and white squares. Parameters that are not
%   uniquely specified by the model are shown with a grey patch.
%
%   Recall that only uniquely estimable parameters can be examined
%   individually with [0,...,0,1,0,...,0] type contrats.
%
%   Surfing the estimability image reports the parameter names and
%   their estimability. Double clicking extracts the estimability
%   vector into the base MatLab workspace.
%
% * Contrast weights plot/image:
%   The weights of the selected contrast(s) are imaged above the design
%   matrix, labelled by their contrast number. t-contrasts are displayed
%   as bar-charts, F-contrasts have their weights matrix c' depicted as a
%   grey-scale image.
%
%   Again, the contrast representation is "surfable": Clicking and
%   dragging over the contrast image reports the contrast name, type and
%   number, and the value of the contrast at the mouse location.
%
% * "Define new contrast..." button:
%   Pops up the contrast definition interface (described below)
%
% * "Reset" button:
%   Clears the current contrast selection.
%
% * "Done" button:
%   Press when the required contrasts have been selected.
%
% * Status line:
%   This indicates how many contrasts have been selected, how
%   multi-contrast selections will be handled, and whether you can press
%   "Done" yet!
%
% In addition, the dialog has a figure ContextMenu, accessed by
% right-clicking in the figure background: In addition to providing
% repeats of the "Define new contrast...", "Reset" & "Done" buttons
% described above, there are two additional options:
%
%   - crash out: this bails out of the contrast manager in a nice way!
%   - rename:    This permits a single contrast to be re-named. You
%            must select the contrast to be renamed before pulling
%                up the context menu for this option to be available.
%
%__________________________________________________________________________
%                                     ConMan: Contrast definition interface
%
% To define a contrast, you must specify:
%  1) a name
%  2) the statistic type: "t-contrast" for SPM{T} or "F-contrast" for SPM{F}
%  3) a valid set of contrast weights
%     (for F-contrasts, this can also be generated given a reduced
%     (design as a partition of the existing design
%
% The contrast definition interface consists of:
%
% * A lilac "name" edit widget for naming the new contrast
%   Type the name of the contrast in this widget.
%   Press return or move the focus to enter the name.
%
% * Radio-buttons for selecting statistic type: "t-contrast" or "F-contrast"
%
% * A large lilac edit widget for entering "contrast weights matrix"
%   - Note that contrast weights should be entered transposed, with
%     contrast weights in rows.
%   - Zero's will be automatically added to the right hand side of
%     contrast weights as needed to give contrast weights matching the
%     number of parameters. For example, if you are interested in
%     contrasting the first two conditions of a design four parameters
%     (design matrix with 4 columns), you need only enter "+1 -1". The
%     contrast manager will parse this to [+1 -1 0 0].
%   - For t-contrasts, this will only accept a single line of input,
%     since contrast weights c' for an SPM{T} must be a row-vector.
%     Pressing <return> or moving the focus (by clicking on another GUI
%     element, such as the "Submit" button) will enter the contrast
%     weights for parsing.
%   - For F-contrasts, the widget accepts multi-line input.
%     Pressing <return> will move to a new line. Press <ctrl>-<return> or
%     move the focus (by clicking on another GUI element, such as the
%     "Submit" button) to enter the contrast weights for parsing.
%   - Text entered in the "contrast weights matrix" is evaluated in the
%     base workspace. Thus, matlab functions can be used in the widget,
%     and base workspace variables can be referred to. (See the help for
%     spm_input.m for more tips on using evaluated input widgets.)
%
% * For F-contrasts, a "columns for reduced design" edit widget appears:
%   - Here you can specify SPM{F}s by specifying the reduced design as
%     a sub-partition of the current design matrix.
%   - Enter the indicies of the design matrix columns you wish to retain
%     in the reduced design.
%   - Pressing <return> or moving the focus (by clicking on another GUI
%     element, such as the "Submit" button) will enter the column indicies
%     for parsing.
%   - An F-contrast weights matrix is constructed from the specified
%     partitioning. (The corresponding F-contrast weights are imaged
%     above the design matrix picture. Double click (or "surf") the
%     contrast image to see the actual values of the F-contrast weights
%     matrix.)
%   - Again, text entered in the "columns for reduced design" is
%     evaluated in the base workspace, permitting the use of functions
%     and variables available in the base workspace.
%   - (Note that the F-contrast weights matrix produced may not be the
%     simplest F-contrast possible for this test, and that the F-contrast
%     weights matrix may not be full rank (e.g. may have two rows where
%     one would do). Nontheless, the F-test is equivalent for the
%     specified partitioning.
%
% * "Submit" button:
%   This button can be used to force parsing of the contrast weights (or
%   columns for reduced design).
%
% * contrast parsing errors pane & contrast parsing info pane:
%   - Once the contrast weights matrix has been entered in the GUI, the
%     inout is parsed.
%   - First, each line is evaluated.
%   - Then, the results for each line are checked to ensure thay define
%     valid contrasts, with trailing zeros added as necessary so the
%     contrast length matches the number of parameters.
%   - Errors encountered during parsing, and invalid contrast weights,
%     are reported in the "contrast parsing errors" pane in red text.
%     Usually the contrast cannot be defined if there are any errors!
%   - Information messages regarding contrast parsing appear in the lower
%     "contrast parsing info" pane in green text.
%   - When defining an F-contrast via "columns for reduced design", the
%     string defining the indicies is similarly parsed, with errors and
%     information messages appearing in the two panes.
%
% * Contrast weights plot/image:
%   Depicts the contrast once valid contrast weights have been specified.
%
% * Image of the design matrix:
%   (As described above for the contrast selection interface)
%
% * Parameter estimability bar
%   (As described above for the contrast selection interface)
%
% * "Reset" button:
%    Resets the contrast definition interface, clearing any contrast
%    currently being defined.
%
% * "Cancel" button:
%   Returns to the contrast selection interface without defining a new
%   contrast.
%
% * "OK" button:
%   Once a valid set of contrast weights has been defined, and the
%   contrast named, pressing "OK" defines the contrast and returns to the
%   contrast selection interface, with the newly defined contrast
%   selected.
%
% * Status line:
%   This indicates progress in contrast definition.
%   Once a valid set of contrast weights have been specified, and a
%   the contrast named, then the status line turns green, indicating
%   that the current contrast can be defined by pressing "OK".
%
%
%==========================================================================
% S P M   C o n t r a s t   m a n a g e m e n t
%==========================================================================
%
% Contrasts are stored by SPM in a single structure (See spm_FcUtil.m
% for the definition and handling of the contrast structure.)
%
% Note that the xCon structure for each contrast contains data specific
% to the current experimental design. Therefore, contrast structures
% can only be copied between analyses (to save re-entering contrasts)
% if the designs are *identical*.
%
% Although the contrasts are named by the user, they are referred to
% internally and on the corresponding contrast, ESS and SPM images (see
% spm_getSPM.m) by their contrast number, which indexes them in the
% order in which they were created. Because of this, it can be rather
% involved to delete any but the most recently defined contrast: All
% file references and indices would have to be canonicalised! Thus, no
% "delete" function is provided (as yet).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_conman.m 7738 2019-12-02 12:45:37Z guillaume $


%==========================================================================
% - FORMAT specifications
%==========================================================================
%( This is a multi function function: If the first argument is a string,)
%( then this is the action string, specifying the particular action     )
%( function to take.                                                    )
%
% FORMAT [I,xCon] = spm_conman(SPM,STATmode,n,Prompt,Mcstr,OK2chg)
%
% SPM.xX        - Design Matrix structure
%               - (see spm_spm.m for structure)
%               - fields used directly are:
%    .xKXs      - space structure of smoothed design matrix
%    .name      - cellstr of parameter names
%
% SPM.xCon (in) - Contrast definitions structure array
%                 (see spm_FcUtil.m for structure, rules & handling)
%               - defaults to empty contrast structure
%               - fields used directly are:
%    .name      - contrast name string
%    .STAT      - character describing statistic required: 'T' or 'F'
%    .c         - contrast weights (column) vector / matrix
%
% STATmode      - string indicating STAT modes to allow contrast
%                 selection/definition for:
%               - 'T' to limit to contrasts defined for SPM{t}
%               - 'F' to limit to contrasts defined for SPM{F}
%               - 'T|F' to allow either contrasts for SPM{t} or SPM{F}
%                 (both may be defined, but only one type may be selected)
%               - 'T&F' to allow both contrasts for SPM{t} and SPM{F}
%               - defaults to 'T|F'
%
% n             - Number of contrasts to select, Inf for unlimited
%
% Prompt        - Prompt string
%
% Mcstr         - string to describe multiple contrast selection
%                 E.g. ' for conjunction' will result in the status message
%                      reading 'Selected 2 contrasts for conjunction' when
%                      two contrasts are selected.
%
% OK2chg        - logical, specifying whether the contrast structure can be
%                 changed. If false, then new contrasts cannot be defined, and
%                 existing contrasts cannot be renamed.
%
% I             - Index (or indices) of contrasts selected
%
% xCon (out)    - Contrast definitions structure array (updated)
%
%                           ----------------
%
% [F,cF] = spm_conman('Initialise',...
%                              Vis,SPM,STATmode,n,Prompt,Mcstr,OK2chg)
% Initialise ConMan GUI for contrast selection/definition
% Vis        - Initialisation action:
%              'close' - closes ConMan window
%              'off'   - hides ConMan window
%              'reset' - hides and resets ConMan window
%              'on'    - initialises ConMan window using arguments given
% SPM.xX     - design matrix structure
% SPM.xCon   - contrast definitions structure array
% STATmode   - string indicating STAT modes to allow contrast
% n          - number of contrasts to select, Inf for unlimited
% Prompt     - Prompt string
% Mcstr      - string to describe multiple contrast selection
% OK2chg     - logical, specifying whether contrast structure can be changed
% F          - figure used for contrast manager window
% cF         - Figure which was current before function call
%
% FORMAT spm_conman('ImDesMtx',xX,h)
% Utility function to display design matrix image & setup "surfing"
% xX         - Design Matrix structure
%            - the first of {xX.nX, xX.xKXs.X} is used for display
% .nX        - Desgin matrix already scaled for display
% .xKXs.X    - temporally filtered design matrix (within space structure)
% .name    - px1 CellStr of parameter names
%
% FORMAT spm_conman('ImParEst',xX,h)
% Utility function to display parameter estimability & setup "surfing"
% xX         - Design Matrix structure
% xX.xKXs    - space structure for K*X
% xX.name    - px1 CellStr of parameter names
% h          - handle of axes to use
%
% FORMAT spm_conman('ListCon',hConList,xCon,STAT,I)
% Utility function to list contrasts in ListBox
% hConList   - handle of GUI ListBox object
% xCon       - current contrast definitions structure array
% STAT       - STAT character: 'T' or 'F' of empty (show all)
% I          - indicies of currently selected contrasts
%
% FORMAT spm_conman('GraphCons',xCon,I,F)
% Utility function to display contrast image/bar-graph & setup "surfing"
% xCon       - contrast definitions structure array
% I          - indicies of contrasts to display
% F          - handle of 'ConMan' figure
%
% FORMAT spm_conman('StatusLine',F,str,col)
% Utility function to update ConMan window statusline
% F          - handle of 'ConMan' figure
% str        - string to display
%              (defaults to status of contrast selection)
% col        - colour to use [defaults to 'w' - white
%
% FORMAT spm_conman('Done_CB')
% CallBack for "Done" button on ConMan contrast selection interface
%
% FORMAT spm_conman('ConList_CB')
% CallBack for contrast selection ListBox
%
% FORMAT STAT = spm_conman('TFA')
% CallBack for 'T','F' or 'any' STAT selection RadioButtons
% FORMAT STAT = spm_conman('TFA',F,STAT)
% Initialisation of 'T','F' or 'any' STAT selection RadioButtons
% FORMAT STAT = spm_conman('TFA',F,STAT,STATmode)
% Function to set STAT & STATmode
% Initialisation of 'T','F' or 'any' STAT selection RadioButtons & STATmode
% F          - handle of 'ConMan' figure
% STAT       - STAT character: 'T' or 'F' of empty (all)
% STATmode   - string indicating STAT modes to allow contrast
%
% FORMAT spm_conman('FConMenu_CB')
% CallBack to set state of ConMan contrast selection figure context menu
%
% FORMAT spm_conman('Rename_CB')
% CallBack to handle contrast renaming
%
% FORMAT [c,I,emsg,imsg,msg] = spm_conman('ParseCon',cstr,X,STAT)
% Contrast weights parser: Catch evaluation errors and invalid contrasts
% cstr       - string (or cellstr) to evaluate to get contrast(s)
% X          - design matrix
% STAT       - 'T' or 'F' (for contrast checking)
% c          - contrast weights matrix
% I          - logical validity indicator: indicates which rows of
%              cellstr(cstr) generated valid contrasts which were
%              included in c
% emsg       - cellstr of error messages produced during parsing
% imsg       - cellstr of information messages for valid contrasts
% msg        - cellstr of all messages produced during parsing,
%              one cell per string in cstr
%
% FORMAT [iX0,I,emsg,imsg] = spm_conman('ParseIStr',str,max)
% DesMtx column index parser: Catch eval errors and invalid indices
% str        - string for evaluation to get column indices
% max        - number of columns in design matrix
% iX0        - vector of column indices: '!' if evaluation error
% I          - 0 if evaluation error, 1 if evaluated OK
% emsg       - error message
%              (evaluation errors, non-integer indices, out of range indices)
% imsg       - information message (valid indices)
%
% FORMAT spm_conman('Reset_CB')
% CallBack handler for "reset" button on contrast selection interface
%
% FORMAT spm_conman('D_Setup_CB')
% CallBack handler for "Define new" button:
% Initialises Contrast Definition interface for use
%
% FORMAT spm_conman('D_ConMtx_CB')
% CallBack handler for contrast weights definition widget
% FORMAT spm_conman('D_X1cols_CB')
% Callback handler for F-contrast "columns for reduced design" widget
%
% FORMAT spm_conman('D_Reset_CB')
% CallBack handler for "reset" button on contrast definition interface
%
% FORMAT spm_conman('D_Cancel_CB')
% CallBack handler for "cancel" button on contrast definition interface
%
% FORMAT spm_conman('D_OK_CB')
% CallBack handler for "OK" button on contrast definiton interface
%
% FORMAT spm_conman('D_Status',F)
% Set status line of contrast definition interface
% F          - Handle of ConMan figure [defaults gcbf]
%
% FORMAT [F,H] = spm_conman('CreateFig')
% Creates ConMan dialog figure (definition interface initially hidden)
% F          - Handle of ConMan figure created
% H          - Stricture of Handles:
% .hConList  - handle of contrast selection ListBox
% .hDesMtxAx - handle of axes for design matrix imaging
% .hPrompt   - handle of prompt text object
% .hSTATmode - handle of frame containing "T/F/All" radio buttons
% .hStatLin  - handle of status line text object (in selection interface)
% .hNew      - handle of "Define new contrast" pushbutton
%__________________________________________________________________________


%-Parameters
%==========================================================================
COLOUR   = [.8,.8,1];   %-Question background colour
PJump    = 1;           %-Jumping of pointer to ConMan window


%==========================================================================
% [I,xCon] = spm_conman(SPM,STATmode,n,Prompt,Mcstr,OK2chg)
%==========================================================================
if (nargin==0) || ~ischar(varargin{1})

    %-Condition arguments
    %----------------------------------------------------------------------
    if nargin<6, OK2chg = 0; else OK2chg=varargin{6}; end
    if nargin<5, Mcstr = ''; else Mcstr=varargin{5}; end
    if nargin<4, Prompt='Select contrast(s)...'; else Prompt=varargin{4}; end
    if nargin<3, n=1; else n=varargin{3}; end
    if nargin<2, STATmode='T|F'; else STATmode=varargin{2}; end
    if nargin<1, error('no SPM struct specified'); else SPM = varargin{1};end

    %----------------------------------------------------------------------
    try
        SPM.xCon;
    catch
        SPM.xCon = {};
    end

    %-Setup ConMan window & jump cursor
    [F,cF] = spm_conman('Initialise','on',SPM,STATmode,n,Prompt,Mcstr,OK2chg);

    if PJump
        PLoc = get(0,'PointerLocation');
        FRec = get(F,'Position');
        set(0,'PointerLocation',[FRec(1)+FRec(3)/2, FRec(2)+FRec(2)/2])
    end

    %-Copy SPM structure for later change detecton, eg renamed contrasts
    tmpSPM = SPM;
    
    %-Wait until filenames have been selected
    hDone = findobj(F,'Tag','Done');
    waitfor(hDone,'UserData')

    %-Exit with error if ConManWin deleted
    if ~ishandle(hDone), error('Contrast Manager was quit!'), end

    %-Get xCon, I & exit status
    status   = get(hDone,'UserData');
    hConList = findobj(F,'Tag','ConList');
    Q        = get(hConList,'UserData');
    I        = Q(get(hConList,'Value'));
    SPM      = get(F,'UserData');
    xCon     = SPM.xCon;
    
    %-Save SPM.mat only if SPM structure as changed
    if spm_check_version('matlab','8.0') >= 0, my_isequaln = @isequaln;
    else my_isequaln = @isequalwithequalnans; end
    if ~my_isequaln(tmpSPM,SPM)
        fprintf('\t%-32s: %30s','Saving SPM.mat','...writing');         %-#
        fmt = spm_get_defaults('mat.format');
        s = whos('SPM');
        if s.bytes > 2147483647, fmt = '-v7.3'; end
        save('SPM.mat', 'SPM', fmt);
        fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...SPM.mat saved')%-#
    end

    %-Reset and hide SelFileWin
    spm_conman('Initialise','off');

    %-Return focus to previous figure (if any)
    set(0,'CurrentFigure',cF)

    %-Jump cursor back to previous location
    if PJump, set(0,'PointerLocation',PLoc), end

    %-Exit with error if reset (status=-1)
    if status == -1, error(['reset: ',mfilename,' bailing out!']), end

    %-Output arguments
    varargout={I,xCon};

    return
end

%==========================================================================
% - Callbacks & Utility embedded functions
%==========================================================================

switch lower(varargin{1})
    %======================================================================
    case 'initialise'
    %======================================================================
    % [F,cF] = spm_conman('Initialise',Vis,SPM,STATmode,n,Prompt,Mcstr,OK2chg)

    if nargin<2, Vis='on'; else Vis=varargin{2}; end

    %-Recover ConMan figure handle (if it exists)
    F  = findobj(get(0,'Children'),'Flat','Tag','ConMan');

    cF = get(0,'CurrentFigure');                  %-Save current figure

    switch lower(Vis)
        %------------------------------------------------------------------
        case 'close'
        %------------------------------------------------------------------
            close(F)
            varargout = {[],cF};
            return
        
        %------------------------------------------------------------------
        case {'off','reset'}
        %------------------------------------------------------------------
            varargout = {F,cF};                   %-Return figure handles
            if isempty(F), return, end            %-Return if no ConMan win
            set(F,'Visible','off')                %-Make window Invisible
            if strcmpi(Vis,'reset')
                set(findobj(F,'Tag','Done'),'UserData',-1)
            end
            return
            
        %------------------------------------------------------------------
        case 'on'
        %------------------------------------------------------------------
            %-Sort out arguments
            %--------------------------------------------------------------
            if nargin<8, OK2chg = 0; else OK2chg=varargin{8}; end
            if nargin<7, Mcstr = ''; else Mcstr=varargin{7}; end
            Mcstr = cellstr(Mcstr); if length(Mcstr)<2, Mcstr{2}=''; end
            if nargin<6, Prompt='Select contrast(s)'; else Prompt=varargin{6}; end
            if nargin<5, n=Inf; else n=varargin{5}; end
            if nargin<4, STATmode='T&F'; else STATmode=varargin{4}; end
            if nargin<3, error('insufficient arguments'), end
            SPM  = varargin{3};
            xCon = SPM.xCon;
            xX   = SPM.xX;

            %-Create/find ConMan window
            %--------------------------------------------------------------
            if isempty(F)                         %-Create ConMan win
                [F,H] = spm_conman('CreateFig');
            else                        %-Get handles
                H.hDesMtxAx = findobj(F,'Tag','DesMtxAx');
                H.hParEstAx = findobj(F,'Tag','ParEstAx');
                H.hConList  = findobj(F,'Tag','ConList');
                H.hPrompt   = findobj(F,'Tag','Prompt');
                H.hTFAf     = findobj(F,'Tag','TFAf');
                H.hSTATmode = findobj(F,'Tag','STATmode');
                H.hStatLin  = findobj(F,'Tag','StatusLine');
                H.hNew      = findobj(F,'Tag','New');
            end

            varargout = {F,cF};                   %-Return figure handles

            %-Store required parameters in UserData of various objects
            %--------------------------------------------------------------
            set(F,          'UserData',SPM)
            set(H.hStatLin, 'UserData',Mcstr)

            %-Initialise interface
            %--------------------------------------------------------------
            set(findobj(F,'Tag','Done'),'UserData',0)   %-Init. Done UserData
            STAT = spm_conman('TFA',F,'',STATmode);     %-Init. TFA buttons
            set(H.hPrompt,'String',Prompt,'UserData',n) %-Set prompt
            spm_conman('ImDesMtx',xX,H.hDesMtxAx)       %-Depict DesMtx
            spm_conman('ImParEst',xX,H.hParEstAx)       %-Parameter estimability
            spm_conman('ListCon',H.hConList,xCon,STAT,[]) %-List contrasts
            if OK2chg, tmp='on'; else tmp='off'; end    %-OK to change xCon?
            set(H.hNew,'Enable',tmp)                    %-En/dis-able change UI
            %-****  if isempty(xCon), spm_conman(); end %-Go straight to DNewUI

            %-Popup figure, retaining CurrentFigure
            %--------------------------------------------------------------
            set(get(findobj(F,'Tag','DefineNew'),'UserData'),'Visible','off')
            %-Hide define UI
            figure(F)                             %-PopUp figure
            set(0,'CurrentFigure',cF)             %-Return to prev. figure
            return

        otherwise
            error('Unrecognised ''Vis'' option')
    end


    %======================================================================
    case 'imdesmtx'
    %======================================================================
        % spm_conman('ImDesMtx',xX,h)

        h  = varargin{3};

        %-Picture design matrix
        %------------------------------------------------------------------
        axes(h)
        if isfield(varargin{2},'nKX') && ~isempty(varargin{2}.nKX)
            hDesMtxIm = image((varargin{2}.nKX+1)*32);
        else
            hDesMtxIm = image(...
                (spm_DesMtx('sca',varargin{2}.xKXs.X,varargin{2}.name)+1)*32);
        end
        set(h,'YTick',[],'XTick',[])                    %-No Tick marks
        set(h,'Tag','DesMtxAx','UserData',varargin{2})  %-Reset axis UserData after image
        xlabel('Design matrix')
        set(hDesMtxIm,'UserData',...
            struct('X',varargin{2}.xKXs.X,'Xnames',{varargin{2}.name}))
        set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')


    %======================================================================
    case 'imparest'
    %======================================================================
        % spm_conman('ImParEst',xX,h)

        xX = varargin{2};
        h  = varargin{3};

        %-Picture design parameter estimability
        %------------------------------------------------------------------
        axes(h)
        est  = spm_SpUtil('IsCon',xX.xKXs);
        nPar = length(est);

        hParEstIm = image((est+1)*32);
        set(h,  'XLim',[0,nPar]+.5,'XTick',[1:nPar-1]+.5,'XTickLabel','',...
            'YLim',[0,1]+.5,'YDir','reverse','YTick',[],...
            'Box','on','TickDir','in','XGrid','on','GridLineStyle','-');
        xlabel('parameter estimability')
        set(h,'Tag','ParEstAx')         %-Reset 'Tag' after image cleared it
        set(hParEstIm,'UserData',struct('est',est,'Xnames',{xX.name}))
        set(hParEstIm,'ButtonDownFcn','spm_DesRep(''SurfEstIm_CB'')')


    %======================================================================
    case 'listcon'
    %======================================================================
        % spm_conman('ListCon',hConList,xCon,STAT,I)

        hConList = varargin{2};
        xCon     = varargin{3};
        STAT     = varargin{4};
        if nargin<5
            Q    = get(hConList,'UserData');
            I    = Q(get(hConList,'Value'));
        else
            I    = varargin{5};
        end

        %-Sort out list, filtering by STAT, and display
        %------------------------------------------------------------------
        if isempty(xCon)
            Q = [];
        elseif isempty(STAT)
            Q = 1:length(xCon);
        else
            Q = find(strcmp({xCon(:).STAT},STAT));
        end

        q = find(ismember(Q,I));

        if ~isempty(Q)
            str        = cell(0);
            for i=1:length(Q)
                str{i} = sprintf('%03d {%c} : %s',...
                    Q(i),xCon(Q(i)).STAT,xCon(Q(i)).name);
            end
            FontAngle  = 'Normal';
            FontWeight = 'Normal';
            Enable     = 'on';
        else
            str        = ['no',deblank([' ',STAT]),' contrasts defined'];
            FontAngle  = 'Italic';
            FontWeight = 'Bold';
            Enable     = 'off';
        end

        set(hConList,'String',str,...
            'UserData',Q,...
            'Value',q,...
            'FontAngle',FontAngle,'FontWeight',FontWeight,...
            'Enable',Enable)

        spm_conman('GraphCons',xCon,Q(q),get(hConList,'Parent'))
        spm_conman('StatusLine',get(hConList,'Parent'))


    %======================================================================
    case 'graphcons'
    %======================================================================
        % spm_conman('GraphCons',xCon,I,F)

        if nargin<2, error('insufficient arguments'), end
        xCon = varargin{2};
        if nargin<3, I=[1:length(xCon)]; else I=varargin{3}; end
        if nargin<4, F=[]; else F=varargin{4}; end
        if isempty(F), F=spm_figure('FindWin','ConMan'); end
        if isempty(F), error('can''t find ConMan win'), end

        cF = get(0,'CurrentFigure');        %-Save current figure
        set(0,'CurrentFigure',F);           %-Make F current

        delete(findobj(F,'Tag','ConGrphAx'));


        %-Display contrasts
        %------------------------------------------------------------------
        if isempty(I)
            axes('Position',[0.65 (0.7 + .1*(1-0.9)) 0.3 .1*.9],...
                'Tag','ConGrphAx','Visible','off')
            text(0.5,0.5,'no contrast(s)',...
                'FontSize',spm('FontSize',8),...
                'FontAngle','Italic',...
                'HorizontalAlignment','Center',...
                'VerticalAlignment','Middle')
        else
            nPar   = size(xCon(1).c,1);
            xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
            nCon   = length(I);
            dy     = 0.2/max(nCon,2);
            hConAx = axes('Position',[0.65 (0.70 + dy*.1) 0.30 dy*(nCon-.1)],...
                'Tag','ConGrphAx','Visible','off');
            title('contrast(s)')
            htxt   = get(hConAx,'title'); set(htxt,'Visible','on')

            for ii = nCon:-1:1
                i  = abs(I(ii));
                axes('Position',[0.65 (0.7 + dy*(nCon-ii+.1)) 0.3 dy*.9])
                if xCon(i).STAT == 'T' && size(xCon(i).c,2)==1
                    %-Single vector contrast for SPM{t} - bar
                    yy = [zeros(1,nPar);repmat(xCon(i).c',2,1);zeros(1,nPar)];
                    h = patch(xx,yy,[1,1,1]*.5);
                    set(gca,'Tag','ConGrphAx',...
                        'Box','off','TickDir','out',...
                        'XTick',[],...
                        'XLim', [0,nPar],...
                        'YTick',[-1,0,+1],'YTickLabel','',...
                        'YLim', [min(xCon(i).c),max(xCon(i).c)] + ...
                        [-1 +1] * max(abs(xCon(i).c))/10    )
                else
                    %-F-contrast - image
                    h = image((xCon(i).c'/max(abs(xCon(i).c(:)))+1)*32);
                    set(gca,'Tag','ConGrphAx',...
                        'Box','on','TickDir','out',...
                        'XTick',[],...
                        'XLim', [0,nPar]+0.5,...
                        'YTick',[0:size(xCon(i).c,2)]+0.5,'YTickLabel','',...
                        'YLim', [0,size(xCon(i).c,2)]+0.5   )
                end
                if I(ii)>0, ylabel(num2str(i)), end
                set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
                    'UserData', struct( 'i',        I(ii),...
                    'h',        htxt,...
                    'xCon',     xCon(i)))
            end
        end

        set(0,'CurrentFigure',cF)  %-Reset CurrentFigure to previous figure


    %======================================================================
    case 'statusline'
    %======================================================================
        % spm_conman('StatusLine',F,str,col)

        if nargin<2, F = findobj(get(0,'Children'),'Flat','Tag','ConMan');
        else         F = varargin{2}; end

        if nargin<3
            n = get(findobj(F,'Tag','Prompt'),'UserData');
            m = length(get(findobj(F,'Tag','ConList'),'Value'));

            str = sprintf('Selected %d contrast',m);
            if m~=1, str=[str,'s']; end

            Mcstr = get(findobj(F,'Tag','StatusLine'),'UserData');
            if m>1, str=[str,Mcstr{1}]; else str=[str,Mcstr{2}]; end

            if m==0
                if n<0
                    str = [str,', press "Done" when finished.'];
                else
                    str = [str,'.'];
                end
            else
                if n==1
                    str = [str,', press "Done".'];
                else
                    str = [str,', press "Done" when finished.'];
                end
            end
        else
            str = varargin{3};
        end

        if nargin<4, col='w'; else col=varargin{4}; end

        set(findobj(F,'Tag','StatusLine'),'String',str,'ForegroundColor',col)


    %======================================================================
    case 'done_cb'
    %======================================================================
        % spm_conman('Done_CB')

        F = gcbf;

        n = get(findobj(F,'Tag','Prompt'),'UserData');
        q = get(findobj(F,'Tag','ConList'),'Value');

        if n>0 && isempty(q)    %-Contrast(s) required, but none selected
            if n==1,            str = 'Select a contrast!';
            elseif isfinite(n), str = sprintf('Select %d contrasts!',n);
            else                str = 'Select at least one contrast!';
            end
        elseif length(q)>abs(n) %-Too many contrasts selected
            if n<0, tstr='at most'; else tstr='only'; end
            if abs(n)==1,       str = ['Select',tstr,' one contrast!'];
            else                str = sprintf('Select %s %d contrasts!',tstr,abs(n));
            end
        elseif n>0 && isfinite(n) && length(q)<n
            if n==1,            str = 'Select a contrast!';
            else                str = sprintf('Select %d contrasts!',n);
            end
        else
            str = '';
        end

        if ~isempty(str)        %-error: display error dialog box
            spm('Beep')
            msgbox(str,sprintf('%s: %s...',spm('Version'),mfilename),...
                'error','modal')
        else                    %-OK, set Done UserData tag for handling
            set(findobj(F,'Tag','Done'),'UserData',1)
        end


    %======================================================================
    case 'conlist_cb'
    %======================================================================
        % spm_conman('ConList_CB')

        F        = gcbf;
        hConList = gcbo;

        Q        = get(hConList,'UserData');
        I        = Q(get(hConList,'Value'));
        SPM      = get(F,'UserData');
        xCon     = SPM.xCon;

        spm_conman('GraphCons',xCon,I,F)
        spm_conman('StatusLine',get(hConList,'Parent'))

        if strcmp(get(F,'SelectionType'),'open'), spm_conman('Done_CB'), end


    %======================================================================
    case {'tfa','d_tf'}
    %======================================================================
        % STAT = spm_conman('TFA')
        % STAT = spm_conman('TFA',F,STAT)
        % STAT = spm_conman('TFA',F,STAT,STATmode)

        D = strcmpi(varargin{1},'d_tf');           %-Handling DefineNew UI?

        if nargin<2, F = gcbf; else F=varargin{2}; end

        if nargin<3              %-Called as CallBack of T or F RadioButton
        %------------------------------------------------------------------
            h = gcbo;
            if get(h,'Value')==0
                %-Was selected - don't allow unselection
                set(h,'Value',1)
                varargout={get(h,'UserData')};
                return
            else
                %-Get new STAT
                STAT = get(h,'UserData');
            end
        else
            STAT = varargin{3};
        end

        if ~(nargin<4)                      %-Called to set STAT & STATmode
        %------------------------------------------------------------------
            STATmode = varargin{4};
            b_set = 1;
        else
            STATmode = get(findobj(F,'Tag','STATmode'),'UserData');
            b_set = 0;
        end

        %-Check STATmode & STAT, set temporary STATmode & STAT indicies
        %------------------------------------------------------------------
        STATinfo = struct(...                         %-STAT info structure
            'mode',     { 'T',  'F',  'T|F',    'T&F'},...
            'vSTAT',    {{'T'},{'F'},{'T','F'},{'T','F',''}},...
            'defSTAT',  { 'T',  'F',  'T',      ''},...
            'Enable',   {   {'on', 'off','off'},...
            {'off','on' ,'off'},...
            {'on', 'on', 'off'},...
            {'on', 'on', 'on' }}    );
        i        = find(strcmp(STATmode,{STATinfo.mode})); %-STATmode index
        if isempty(i), error('unknown STATmode'), end      %-Check STATmode valid
        if D && i==4, i=3; end                          %-Treat 'T&F' as 'T|F'?
        if isempty(STAT), STAT=STATinfo(i).defSTAT; end %-Set STAT as default(?)
        j        = find(strcmp(STAT,{'T','F',''}));     %-STAT index
        if isempty(j), error('unknown STAT'); end       %-Check known STAT
        if ~any(strcmp(STAT,STATinfo(i).vSTAT))         %-Check STAT is
            error('Invalid STAT for this STATmode')     % valid for
        end                                             % this STATmode

        %-Set STAT buttons (& Dis/Enable according to STATmode if b_setEnable)
        %------------------------------------------------------------------
        if ~D, Tag='TFA'; else Tag='D_TF'; end
        H = flipud(findobj(F,'Tag',Tag));
        set(H(j),           'Value',1)
        set(H(setdiff([1:length(H)],j)),'Value',0)
        if b_set
            %-Set RadioButton 'Enable' & store STATmode
            set(findobj(F,'Tag','STATmode'),'UserData',STATmode)
            for k=1:length(H), set(H(k),'Enable',STATinfo(i).Enable{k}), end
        end

        if ~D               %-Additional UI setup for main contrast manager
        %------------------------------------------------------------------

            %-Reset ConList for new STAT if callback for TFA button
            %--------------------------------------------------------------
            if nargin<3
                SPM      = get(F,'UserData');
                xCon     = SPM.xCon;
                spm_conman('ListCon',findobj(F,'Tag','ConList'),xCon,STAT)
            end

        else  %-Additional UI setup for  DNew contrast definition interface
        %------------------------------------------------------------------

            %-Get handles of control objects
            %--------------------------------------------------------------
            hD_ConMtx  = findobj(F,'Tag','D_ConMtx');
            hD_X1cols  = findobj(F,'Tag','D_X1cols');
            hD_ConErrs = findobj(F,'Tag','D_ConErrs');
            hD_ConInfo = findobj(F,'Tag','D_ConInfo');
            HD_Ttxt    = findobj(F,'Tag','D_Ttxt');
            HD_Ftxt    = findobj(F,'Tag','D_Ftxt');

            %-Set interface for new STAT
            %--------------------------------------------------------------
            set(hD_ConMtx,'String',{},'UserData',[])            %-Clear ConMtx box
            set(hD_X1cols,'String','')                          %-Clear X1cols box
            set([hD_ConErrs,hD_ConInfo],'String',{},'Value',[]) %-Clear parsing boxes
            spm_conman('GraphCons',[],[],F)                     %-Clear contrast plot
            spm_conman('D_Status',F)                            %-Set StatusLine

            switch STAT
                case 'T'
                    set(hD_ConMtx,'Max',1)
                    set(HD_Ttxt,'Visible','on')
                    set([hD_X1cols;HD_Ftxt],'Visible','off')
                case 'F'
                    set(hD_ConMtx,'Max',2)
                    set(HD_Ttxt,'Visible','off')
                    set([hD_X1cols;HD_Ftxt],'Visible','on')
                otherwise
                    error('illegal case')
            end

        end

        %-Return STAT
        %------------------------------------------------------------------
        varargout = {STAT};


    %======================================================================
    case 'fconmenu_cb'
    %======================================================================
        % spm_conman('FConMenu_CB')

        if strcmp(get(findobj(gcbf,'Tag','New'),'Enable'),'off')
            set(findobj(gcbo,'Tag','CM_New'),'Enable','off')
            set(findobj(gcbo,'Tag','CM_Ren'),'Enable','off')
        else
            set(findobj(gcbo,'Tag','CM_New'),'Enable','on')
            if length(get(findobj(gcbf,'Tag','ConList'),'Value'))==1
                set(findobj(gcbo,'Tag','CM_Ren'),'Enable','on')
            else
                set(findobj(gcbo,'Tag','CM_Ren'),'Enable','off')
            end
        end


    %======================================================================
    case 'rename_cb'
    %======================================================================
        % spm_conman('Rename_CB')

        F = gcbf;

        hConList = findobj(F,'Tag','ConList');
        Q        = get(hConList,'UserData');
        i        = get(hConList,'Value');

        %-Make sure there's only one selected contrast
        %------------------------------------------------------------------
        if length(i)~=1
            msgbox('Can only rename a single contrast!',...
                sprintf('%s: %s...',spm('Version'),mfilename),...
                'error','modal')
            return
        end

        %-Get contrast structure array and indicies of current contrast
        %------------------------------------------------------------------
        SPM      = get(F,'UserData');
        xCon     = SPM.xCon;
        I      = Q(i);

        %-Get new name
        %------------------------------------------------------------------
        str = sprintf('Enter new name for contrast %d (currently "%s"):',I,xCon(I).name);
        nname  = inputdlg(str,'SPM: Rename contrast',1,{xCon(I).name},'off');
        if isempty(nname), return, end

        %-Change name in ConList
        %------------------------------------------------------------------
        tmp    = get(hConList,'String');
        tmp{i} = strrep(tmp{i},xCon(I).name,nname{1});
        set(hConList,'String',tmp)

        %-Change name in contrast structure
        %------------------------------------------------------------------
        xCon(I).name = nname{1};
        SPM.xCon = xCon;
        set(F,'UserData',SPM)


    %======================================================================
    case 'parsecon'                      %-Parse (cell)string into contrast
    %======================================================================
        % [c,I,emsg,imsg,msg] = spm_conman('ParseCon',cstr,X,STAT)
        % X is raw DesMtx or space structure

        %-Sort out parameters
        %------------------------------------------------------------------
        if nargin<4, STAT='F'; else STAT=varargin{4}; end
        cstr = varargin{2};
        X    = varargin{3};
        p    = spm_SpUtil('size',X,2);

        %-Preliminary parsing of contrast string (cstr) into cellstr
        %------------------------------------------------------------------
        if isempty(cstr)
            varargout = {[],0,{'    <- !empty input'},{},{}};
            return
        elseif iscell(cstr)
            if numel(cstr)~=1, cstr=cstr(:); end
            c    = cstr;
        elseif ischar(cstr)
            cstr = cellstr(cstr);
            c    = cstr;
        elseif isnumeric(cstr)
            if ndims(cstr)>2, error('matrices only!'), end
            c    = num2cell(cstr,2);
            cstr = cell(size(c));
            for i=1:numel(c), cstr{i}=num2str(c{i}); end
        else
            error('contrast input must be string or number')
        end

        %-Evaluate individual lines of contrast matrix input
        %------------------------------------------------------------------
        I   = zeros(size(c,1),1);
        msg = cell(size(c)); [msg{:}] = deal(' (OK)');
        for i=1:size(c,1)
            if ischar(c{i})
                c{i} = evalin('base',['[',c{i},']'],'''!''');
            end
            if ischar(c{i})
                msg{i} = '!evaluation error';
            else
                if isempty(c{i})
                    msg{i}=' (empty line  - ignored)';
                elseif all(c{i}(:)==0)
                    if size(c{i},1)==1, str='vector'; else str='matrix'; end
                    c{i}=[]; msg{i}=[' (zero ',str,' - ignored)'];
                elseif STAT=='T' && size(c{i},1)>1
                    c{i}='!'; msg{i}='!vector required';
                elseif size(c{i},2)>p
                    c{i}='!'; msg{i}=sprintf('!too long - only %d prams',p);
                else
                    if size(c{i},2)<p
                        tmp = p-size(c{i},2);
                        c{i}=[c{i}, zeros(size(c{i},1),tmp)];
                        if size(c{i},1)>1, str=' column'; else str=''; end
                        if tmp>1,          str=[str,'s']; end
                        msg{i} = sprintf(' (right padded with %d zero%s)',tmp,str);
                    end
                    if ~spm_SpUtil('allCon',X,c{i}')
                        c{i}='!'; msg{i}='!invalid contrast';
                    end
                end
            end
            I(i)=~ischar(c{i});
        end

        %-Construct contrast matrix, validity indicator, and collate parsing messages
        %------------------------------------------------------------------
        c    = cat(1,c{find(I)});
        msg  = [char(cstr), repmat('  <-  ',size(msg,1),1), char(msg)];
        emsg = msg(find(~I),:);
        imsg = msg(find( I),:);

        if all(I) && STAT=='T' && size(c,1)>1 %-Check for vector t-contrasts!
            I=zeros(size(I)); emsg={'!vector required'}; imsg={};
        end

        %-Return arguments
        %------------------------------------------------------------------
        varargout = {c',I,emsg,imsg,msg};



    %======================================================================
    case 'parseistr'                                   %-Parse index string
    %======================================================================
        % [iX0,I,emsg,imsg] = spm_conman('ParseIStr',str,max)
        % str should be a string (row)vector

        %-Sort out parameters
        %------------------------------------------------------------------
        str = varargin{2};
        mx  = varargin{3};

        %-Process input string
        %------------------------------------------------------------------
        I = evalin('base',['[',str,']'],'''!''');

        if ischar(I)
            varargout = {'!',0,[str,'  <- !evaluation error'],''};
            return
        end

        %-Construct list of valid indicies
        %------------------------------------------------------------------
        ok  = ismember(I(:)', 1:mx);
        iX0 = I(ok);

        %-Construct diagnostic info messages
        %------------------------------------------------------------------
        str = {}; msg={};
        if any(ok)
            str = [str; num2str(I(ok))];
            msg = [msg; '  <-  (OK)'];
        end
        tmp = ( I<1 | I>mx );           %-Out of range
        if any(tmp)
            str = [str; num2str(I(tmp))];
            msg = [msg; sprintf('  <-  (ignored - not in [1:%d]',mx)];
        end
        tmp = ( ~tmp & ~ok );           %-Non integer in range
        if any(tmp)
            str = [str; num2str(I(tmp))];
            msg = [msg; '  <-  (ignored - non-integer)'];
        end

        %-Return arguments
        %------------------------------------------------------------------
        varargout = {iX0,1,'',cellstr([char(str),char(msg)])};


    %======================================================================
    case 'reset_cb'
    %======================================================================
        % spm_conman('Reset_CB')

        hConList = findobj(gcbf,'Tag','ConList');
        SPM     = get(gcf,'UserData');
        xCon = SPM.xCon;
        STAT     = get(findobj(gcbf,'Tag','TFA','Value',1),'UserData');

        spm_conman('ListCon',hConList,xCon,STAT,[])


    %======================================================================
    case 'd_setup_cb'
    %======================================================================
        % spm_conman('D_Setup_CB')

        F        = gcbf;
        STAT     = get(findobj(F,'Tag','TFA','Value',1),'UserData');
        STATmode = get(findobj(F,'Tag','STATmode'),'UserData');

        set(F,'UIContextMenu',[])                 %-Disable Fig ContextMenu
        H = get(findobj(F,'Tag','DefineNew'),'UserData');   %-Get define UI handles
        set(findobj(H,'flat','Tag','D_name'),'String','')   %-Clear name
        %set(findobj(H,'flat','Tag','D_ConMtx'),'UserData',[])  %-Clear con definition
        set(H,'Visible','on')                     %-Show UI
        SPM = get(F, 'UserData');

        spm_conman('D_TF',F,STAT,STATmode);       %-Setup rest of define UI


    %======================================================================
    case {'d_conmtx_cb','d_x1cols_cb'}
    %======================================================================
        % spm_conman('D_ConMtx_CB')
        % spm_conman('D_X1cols_CB')

        fcn = find(strcmpi(varargin{1},{'d_conmtx_cb','d_x1cols_cb'}));

        F   = findobj('Tag', 'ConMan');
        hD_ConMtx  = findobj(F,'Tag','D_ConMtx');
        hD_X1cols  = findobj(F,'Tag','D_X1cols');

        if strcmpi(get(hD_ConMtx, 'Enable'), 'on')
            if fcn == 1
                str = get(hD_ConMtx,'String');
            elseif fcn == 2
                str = get(hD_X1cols,'String');
            else
                % This shouldn't happen
            end
        else
            % This shouldn't happen
        end

        %-Extract info from holding objects
        %------------------------------------------------------------------
        xX   = get(findobj(F,'Tag','DesMtxAx'),'UserData');
        STAT = get(findobj(F,'Tag','D_TF','Value',1),'UserData');


        if fcn==1                         %-Parse string from ConMtx widget
        %------------------------------------------------------------------
            set(hD_X1cols,'String','')
            
            [c,I,emsg,imsg] = spm_conman('ParseCon',str,xX.xKXs,STAT);
            if all(I)
                DxCon = spm_FcUtil('Set','',STAT,'c',c,xX.xKXs);
            else
                DxCon = [];
            end

        elseif fcn==2          %-Process column indicies from X1cols widget
        %------------------------------------------------------------------
            set(hD_ConMtx,'String','')

            nPar              = spm_SpUtil('size',xX.xKXs,2);
            [iX0,I,emsg,imsg] = spm_conman('ParseIStr',str,nPar);

            if I
                try %-try-catch block for any errors in spm_FcUtil!
                    DxCon = spm_FcUtil('Set','',STAT,'iX0',iX0,xX.xKXs);
                    if STAT=='T' && size(DxCon.c,2)>1
                        I = 0; emsg = {'! t-contrasts must be vectors'};
                    end
                catch
                    I    = 0;
                    emsg = lasterr;
                end
            end
        end

        %-Define the contrast or report errors...
        %------------------------------------------------------------------
        set(findobj(F,'Tag','D_ConErrs'),'String',emsg,'Value',[])
        set(findobj(F,'Tag','D_ConInfo'),'String',imsg,'Value',[])
        if all(I)
            set(hD_ConMtx,'UserData',DxCon);        %-Store contrast
            spm_conman('GraphCons',DxCon,-1,F)      %-Depict contrast
        else
            set(hD_ConMtx,'UserData',[]);           %-Clear contrast store
            spm_conman('GraphCons',[],[],F)         %-Clear contrast plot
        end
        spm_conman('D_Status',F)                    %-Set StatusLine


    %======================================================================
    case 'd_reset_cb'
    %======================================================================
        % spm_conman('D_Reset_CB')

        STAT = get(findobj(gcbf,'Tag','TFA','Value',1),'UserData');
        set(findobj(gcbf,'Tag','D_name'),'String','')     %-Clear name
        set(findobj(gcbf,'Tag','D_ConMtx'),'UserData',[]) %-Contrast definition
        spm_conman('D_TF',gcbf,STAT);               %-Setup rest of define UI


    %======================================================================
    case 'd_cancel_cb'
    %======================================================================
        % spm_conman('D_Cancel_CB')

        set(get(findobj(gcbf,'Tag','DefineNew'),'UserData'),'Visible','off')
        set(gcbf,'UIContextMenu',findobj(gcbf,'Tag','ConMan_ConMen'))
        spm_conman('StatusLine')


    %======================================================================
    case 'd_ok_cb'
    %======================================================================
        % spm_conman('D_OK_CB')

        F = gcbf;

        name  = get(findobj(F,'Tag','D_name'),'String');
        DxCon = get(findobj(F,'Tag','D_ConMtx'),'UserData');
        STAT  = get(findobj(F,'Tag','D_TF','Value',1),'UserData');

        dNam = ~isempty(name);
        dCon = ~isempty(DxCon);

        if ~(dNam && dCon)
            spm('Beep')
            str = { 'contrast name not defined!','',...
                'no valid contrast defined!',''};
            msgbox(str([dNam+1,dCon+3]),...
                sprintf('%s: %s...',spm('Version'),mfilename),...
                'error','modal')
            return
        end

        %-Append new contrast to xCon structure of ConMan figure 'UserData'
        %------------------------------------------------------------------
        DxCon.name = name;
        if ~strcmp(DxCon.STAT,STAT), error('STAT & DxCon.STAT mismatch!'), end
        SPM      = get(F,'UserData');
        xCon     = SPM.xCon;
        
        % Due to an earlier bug in SPM, a spurious field 
        % PSTAT may have been created. Remove this if it exists.
        if isfield(xCon,'PSTAT')
            xCon=rmfield(xCon,'PSTAT');
        end
            
        if isempty(xCon)
            xCon = DxCon;
        else
            xCon = [xCon, DxCon];
        end
        SPM.xCon = xCon;
        set(F,'UserData',SPM);
        
        %-Redisplay the new list of contrasts, with the new one selected
        %------------------------------------------------------------------
        hConList = findobj(F,'Tag','ConList');
        I = length(xCon);

        %-Use this code to add the new contrast to a multiple selection, if allowed
        % Q        = get(hConList,'UserData');
        % I        = Q(get(hConList,'Value'));
        % n        = get(findobj(F,'Tag','Prompt'),'UserData');
        % if abs(n)>1, I=[I,length(xCon)]; else, I=length(xCon); end

        spm_conman('TFA',F,xCon(end).STAT);                     %-Set STAT
        spm_conman('ListCon',hConList,xCon,xCon(end).STAT,I)    %-ListCon

        %-Hide the DefineNew UI
        %------------------------------------------------------------------
        set(get(findobj(gcbf,'Tag','DefineNew'),'UserData'),'Visible','off')
        set(gcbf,'UIContextMenu',findobj(gcbf,'Tag','ConMan_ConMen'))

        
    %======================================================================
    case 'd_status'
    %======================================================================
        % spm_conman('D_Status',F)

        if nargin<2, F=gcbf; else F=varargin{2}; end
        str  = {' not',''};
        dNam = ~isempty(get(findobj(F,'Tag','D_name'),'String'));
        dCon = ~isempty(get(findobj(F,'Tag','D_ConMtx'),'UserData'));
        if dNam && dCon, ok='on'; col='g'; else ok='off'; col='w'; end
        spm_conman('StatusLine',F,...
            sprintf('name%s defined, contrast%s defined',str{dNam+1},str{dCon+1}),...
            col)
        %set(findobj(F,'Tag','D_OK'),'Enable',ok)     %-Enable "OK" button?


    %======================================================================
    case 'createfig'
    %======================================================================
        % [F,H] = spm_conman('CreateFig')
        % Handle Structure - H.{hConList,hDesMtxAx,hPrompt,hSTATmode,hStatLin,hNew}

        cF = get(0,'CurrentFigure');        %-Save current figure

        %-Create window, compute scaling for screen size
        %------------------------------------------------------------------
        WS = spm('WinScale');               %-Window scaling factors
        FS = spm('FontSizes');              %-Scaled font sizes
        PF = spm_platform('fonts');         %-Font names (for this platform)
        S0 = spm('WinSize','0',1);          %-Screen size (of the current monitor)

        F = figure('IntegerHandle','off',...
            'Tag','ConMan',...
            'Name','SPM contrast manager','NumberTitle','off',...
            'Position',[S0(1)+S0(3)/2,S0(2)+S0(4)/2,0,0] + [-250,-200,500,400].*WS,...
            'Resize','off',...
            'Color',[1 1 1]*.7,...
            'MenuBar','none',...
            'DefaultTextColor','k',...
            'DefaultTextFontName',PF.helvetica,...
            'DefaultTextFontSize',FS(10),...
            'DefaultAxesFontName',PF.helvetica,...
            'DefaultUicontrolBackgroundColor',[1 1 1]*.7,...
            'DefaultUicontrolFontName',PF.helvetica,...
            'DefaultUicontrolFontSize',FS(10),...
            'DefaultUicontrolInterruptible','on',...
            'Colormap',gray(64),...
            'Renderer',spm_get_defaults('renderer'),...
            'Visible','off');


        %-Draw GUI objects
        %------------------------------------------------------------------
        hPrompt = uicontrol(F,'Style','Text','Tag','Prompt','UserData',[],...
            'String','<Prompt not set yet>',...
            'FontName',PF.times,...
            'FontWeight','Bold',...
            'FontAngle','Italic',...
            'FontSize',FS(16),...
            'ForegroundColor','k',...
            'BackgroundColor',[1,1,1]*.7,...
            'HorizontalAlignment','Center',...
            'Position',[020 370 280 025].*WS);

        %                           ----------------
        %-T/F/all buttons...

        hSTATmode = uicontrol(F,'Style','Frame','Tag','STATmode',...
            'Position',[040 340 260 028].*WS);
        uicontrol(F,'Style','Text','String','show',...
            'FontName',PF.times,'FontAngle','Italic','FontSize',FS(8),...
            'ForegroundColor','w',...
            'Position',[045 365 030 010].*WS);
        hT = uicontrol(F,'Style','RadioButton','String','t-contrasts','Tag','TFA',...
            'ToolTipString','...to show only contrasts for SPM{t}',...
            'FontSize',FS(9),...
            'ForegroundColor','k',...
            'UserData','T',...
            'Position',[044 343 105 020].*WS);
        hF = uicontrol(F,'Style','RadioButton','String','F-contrasts','Tag','TFA',...
            'ToolTipString','...to show only contrasts for SPM{F}',...
            'FontSize',FS(9),...
            'ForegroundColor','k',...
            'UserData','F',...
            'Position',[149 343 105 020].*WS);
        hA = uicontrol(F,'Style','RadioButton','String','all','Tag','TFA',...
            'ToolTipString','...to show all defined contrasts',...
            'FontSize',FS(9),...
            'ForegroundColor','k',...
            'UserData','',...
            'Position',[254 343 041 020].*WS);
        set([hT,hF,hA],'CallBack','spm_conman(''TFA'');',...
            'Interruptible','off','BusyAction','Queue')


        %                           ----------------
        %-Contrast list...

        uicontrol(F,'Style','Text','Tag','ConListHdr',...
            'String','### {type} : name',...
            'FontSize',FS(8),'FontAngle','Italic',...
            'HorizontalAlignment','Left',...
            'BackgroundColor','w',...
            'Position',[042 320 256 016].*WS);
        
        hh = uicontextmenu('Tag','ConMan_ConMenu');
        uimenu(hh,'Label','Rename selected contrast...',...
            'Tag','CM_Ren',...
            'CallBack', 'spm_conman(''Rename_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        set(hh,'CallBack','spm_conman(''FConMenu_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        
        hConList = uicontrol(F,'Style','ListBox','Tag','ConList',...
            'ToolTipString',['Select contrast(s) - drag/shift-click',...
            '/ctrl-click to select multiple contrasts'],...
            'String',{'list','not','set','yet'},...
            'Max',2,...
            'CallBack','spm_conman(''ConList_CB'')',... 
            'UIContextMenu',hh,...
            'Interruptible','off','BusyAction','Queue',...
            'BackgroundColor','w',...
            'Position',[040 080 260 240].*WS);

        %                           ----------------
        %-Control buttons & status area...

        hNew = uicontrol(F,'Style','Pushbutton','String','Define new contrast...',...
            'Tag','New',...
            'ToolTipString','define new contrast',...
            'ForegroundColor','b',...
            'Callback','spm_conman(''D_Setup_CB'')',...
            'Enable','on',...
            'Position',[040 050 150 022].*WS);

        uicontrol(F,'Style','Pushbutton','String','Reset',...
            'ToolTipString','reset selection',...
            'ForegroundColor','r',...
            'Callback','spm_conman(''Reset_CB'')',...
            'Position',[195 050 050 022].*WS);

        uicontrol(F,'Style','Pushbutton','String','Done',...
            'ToolTipString','done - press when selected contrast(s)',...
            'ForegroundColor','m',...
            'Tag','Done','UserData',1,...
            'Callback','spm_conman(''Done_CB'')',...
            'Interruptible','off','BusyAction','Cancel',...
            'Position',[250 050 050 022].*WS);

        uicontrol(F,'Style','Frame','Tag','StatusArea',...
            'Position',[010 010 480 030].*WS);

        hStatLin = uicontrol(F,'Style','Text','Tag','StatusLine',...
            'String','<Not set yet>',...
            'FontAngle','Italic',...
            'HorizontalAlignment','Center',...
            'ForegroundColor','w',...
            'Position',[020 015 430 020].*WS);

        %                           ----------------
        %-Axes for design matrix and parameter estimability...

        hDesMtxAx = axes('Parent',F,'Tag','DesMtxAx',...
            'Position',[0.65 0.30 0.30 0.40],...
            'Color','w',...
            'Box','on','XTick',[],'YTick',[]);

        hParEstAx = axes('Parent',F,'Tag','ParEstAx',...
            'Position',[0.65 0.18 0.30 0.05],...
            'Color','w',...
            'Box','on','XTick',[],'YTick',[]);

        %                           ----------------
        %-Figure UICOntextMenu

        h = uicontextmenu('Tag','ConMan_ConMen');
        set(F,'UIContextMenu',h)
        uimenu(h,'Label','Define new contrast...',...
            'Tag','CM_New',...
            'CallBack','spm_conman(''D_Setup_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        uimenu(h,'Label','Rename selected contrast...',...
            'Tag','CM_Ren',...
            'CallBack','spm_conman(''Rename_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        uimenu(h,'Label','Reset','Separator','on',...
            'CallBack','spm_conman(''Reset_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        uimenu(h,'Label','Done',...
            'CallBack','spm_conman(''Done_CB'')',...
            'Interruptible','off','BusyAction','Cancel');
        uimenu(h,'Label','crash out','Separator','on',...
            'CallBack','spm_conman(''Initialise'',''reset'');',...
            'Interruptible','off','BusyAction','Cancel');
        set(h,'CallBack','spm_conman(''FConMenu_CB'')',...
            'Interruptible','off','BusyAction','Cancel');


        %-Draw contrast definition GUI
        %------------------------------------------------------------------
        H = [];                     %-Save handles for visibility switching

        h = uicontrol(F,'Style','Frame','Tag','DefineNew',...
            'Position',[010 045 300 350].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','Tag','D_Prompt','UserData',[],...
            'String','define contrast...',...
            'FontName',PF.times,...
            'FontWeight','Bold',...
            'FontAngle','Italic',...
            'FontSize',FS(14),...
            'ForegroundColor','b',...
            'HorizontalAlignment','Center',...
            'Position',[020 360 280 030].*WS);
        H = [H,h];

        %                           ----------------
        %-name
        h = uicontrol(F,'Style','Frame','Position',[020 335 280 033].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','name',...
            'FontSize',FS(10),...
            'FontAngle','Italic',...
            'ForegroundColor','w',...
            'HorizontalAlignment','Center',...
            'Position',[025 355 045 020].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Edit','Tag','D_name',...
            'ToolTipString','enter name for contrast',...
            'HorizontalAlignment','Left',...
            'BackgroundColor',COLOUR,...
            'Callback','spm_conman(''D_Status'')','Interruptible','off',...
            'Position',[080 340 215 022].*WS);
        H = [H,h];

        %                           ----------------
        %-type
        h = uicontrol(F,'Style','Frame','Position',[020 295 280 033].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','type',...
            'FontSize',FS(10),...
            'FontAngle','Italic',...
            'ForegroundColor','w',...
            'HorizontalAlignment','Center',...
            'Position',[025 315 040 020].*WS);
        H = [H,h];

        hDT = uicontrol(F,'Style','RadioButton','String','t-contrast','Tag','D_TF',...
            'ToolTipString','...to define contrast for SPM{t}',...
            'FontSize',FS(9),...
            'ForegroundColor','k',...
            'UserData','T',...
            'Position',[080 300 105 022].*WS);
        hDF = uicontrol(F,'Style','RadioButton','String','F-contrast','Tag','D_TF',...
            'ToolTipString','...to define contrast for SPM{F}',...
            'FontSize',FS(9),...
            'ForegroundColor','k',...
            'UserData','F',...
            'Position',[190 300 105 022].*WS);
        set([hDT,hDF],'CallBack','spm_conman(''D_TF'');',...
            'Interruptible','off','BusyAction','Queue')
        H = [H,hDT,hDF];

        %                           ----------------
        %-contrast
        h = uicontrol(F,'Style','Frame','Position',[020 080 280 208].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','contrast',...
            'FontSize',FS(10),...
            'FontAngle','Italic',...
            'ForegroundColor','w',...
            'HorizontalAlignment','Center',...
            'Position',[025 275 055 020].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','contrast',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[030 255 045 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','weights',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[030 245 045 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','vector','Tag','D_Ttxt',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[030 235 045 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','matrix','Tag','D_Ftxt',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[030 235 045 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Edit','Tag','D_ConMtx',...
            'ToolTipString','enter contrast',...
            'HorizontalAlignment','Left',...
            'BackgroundColor',COLOUR,...
            'Max',2,...
            'CallBack','spm_conman(''D_ConMtx_CB'')',...
            'Interruptible','off','BusyAction','Queue',...
            'UserData',[],...
            'Position',[080 200 215 082].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','or','Tag','D_Ftxt',...
            'FontAngle','Italic',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Left',...
            'Position',[025 205 030 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','columns for','Tag','D_Ftxt',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[022 190 070 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Text','String','reduced design','Tag','D_Ftxt',...
            'FontSize',FS(6),...
            'HorizontalAlignment','Right',...
            'Position',[022 180 070 008].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Edit','Tag','D_X1cols',...
            'ToolTipString',...
            'enter column indicies of reduced design matrix X0',...
            'HorizontalAlignment','Left',...
            'BackgroundColor',COLOUR,...
            'CallBack','spm_conman(''D_X1cols_CB'')',...
            'Interruptible','off','BusyAction','Queue',...
            'Position',[090 180 155 020].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Pushbutton','String','...submit',...
            'Tag', 'submit',...
            'FontSize',FS(8),...
            'ForegroundColor','c',...
            'Position',[245 180 050 020].*WS);
        H = [H,h];

        %                           ----------------
        %-Errors & info boxes...
        h = uicontrol(F,'Style','ListBox','Tag','D_ConErrs',...
            'ToolTipString','contrast parsing errors',...
            'FontName',PF.courier,'FontSize',FS(7),...
            'ForegroundColor','r',...
            'BackgroundColor',[1 1 1]*.7,...
            'Enable','on','Max',2,'Value',[],...
            'Position',[027 127 268 042].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','ListBox','Tag','D_ConInfo',...
            'ToolTipString','contrast parsing info',...
            'FontName',PF.courier,'FontSize',FS(7),...
            'ForegroundColor','g',...
            'BackgroundColor',[1 1 1]*.7,...
            'Enable','on','Max',2,'Value',[],...
            'Position',[027 085 268 042].*WS);
        H = [H,h];

        %                           ----------------
        %-Control buttons & status area...
        h = uicontrol(F,'Style','Pushbutton','String','Reset',...
            'Tag','D_Reset',...
            'ToolTipString','reset definition interface',...
            'ForegroundColor','b',...
            'Callback','spm_conman(''D_Reset_CB'')',...
            'Interruptible','off','BusyAction','Cancel',...
            'Position',[140 053 050 022].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Pushbutton','String','Cancel',...
            'Tag','D_Cancel',...
            'ToolTipString','cancel contrast definition',...
            'ForegroundColor','r',...
            'Callback','spm_conman(''D_Cancel_CB'')',...
            'Interruptible','off','BusyAction','Cancel',...
            'Position',[195 053 050 022].*WS);
        H = [H,h];
        h = uicontrol(F,'Style','Pushbutton','String','OK',...
            'Tag','D_OK',...
            'ToolTipString','OK - press to accept newly defined contrast',...
            'ForegroundColor','m',...
            'Callback','spm_conman(''D_OK_CB'')',...
            'Interruptible','off','BusyAction','Cancel',...
            'Position',[250 053 050 022].*WS);
        H = [H,h];

        set(findobj(H,'flat','Tag','DefineNew'),'UserData',H)

        %-Finish up
        %------------------------------------------------------------------
        set(0,'CurrentFigure',cF)
        varargout = {F,struct(  'hConList', hConList,...
            'hDesMtxAx',    hDesMtxAx,...
            'hParEstAx',    hParEstAx,...
            'hPrompt',  hPrompt,...
            'hSTATmode',    hSTATmode,...
            'hStatLin', hStatLin,...
            'hNew',     hNew    )};


    %======================================================================
    otherwise                                              %-unknown action
    %======================================================================
        error(['Illegal Action string: ',varargin{1}])


%==========================================================================
end                                                               % - E N D
%==========================================================================
