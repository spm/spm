function varargout = spm_input(varargin)
% Comprehensive graphical and command line input function
% FORMAT [p,YPos] = spm_input(Prompt,YPos,Type,Labels,Values,Default)
%_______________________________________________________________________
%
% spm_input handles most forms of user input for SPM. (File selection
% is handled by spm_get.m) This can be either via GUI with keyboard
% shortcuts (in the SPM "Interactive" window), or on the command line
% (in the MatLab command window) if requested.
%
% There are five types of input: String, Evaluated, Conditions, Buttons
% and Menus:  These prompt for string input; string input which is
% evaluated to give a numerical result; selection of one item from a
% set of buttons; selection of an item from a menu.
%
% - STRING, EVALUATED & CONDITION input -
% For STRING, EVALUATED and CONDITION input types, a prompt is
% displayed adjacent to an editable text entry widget (with a lilac
% background!).  This entry widget will initially contain any default
% reply. Pressing <RETURN> or <ENTER> will accept the default. Clicking
% in the entry widget allows editing, pressing <RETURN> or <ENTER>
% enters the result. You must enter something, empty answers are not
% accepted.
%
% Basic editing of the entry widget is supported *without* clicking in
% the widget. This enables you to type responses to a sequence of
% questions without having to repeatedly click the mouse. Supported are
% BackSpace and Delete, line kill (^U). Other standard ASCII characters
% are appended to the text in the entry widget. Press <RETURN> or <ENTER>
% to submit your response.
%
% A ContextMenu is provided (in the figure background) giving access to
% relevant utilities including the facility to load input from a file
% (see spm_load.m and examples given below): Click the right button on
% the figure background.
%
% For EVALUATED input, the string submitted is evaluated in the base
% MatLab workspace (see MatLab's `eval` command) to give a numerical
% value. This permits the entry of numerics, matrices, expressions,
% functions or workspace variables. I.e.:
%   i)  - a number, vector or matrix        e.g. "[1 2 3 4]"
%                                                "[1:4]"
%                                                "1:4"
%  ii)  - an expression                     e.g. "pi^2"
%                                                "exp(-[1:36]/5.321)"
% iii)  - a function (that will be invoked) e.g. "spm_load('tmp.dat')"
%         (function must be on MATLABPATH)       "input_cov(36,5.321)"
%  iv)  - a variable from the base workspace
%                                           e.g. "tmp"
%
% The last three options provide a great deal of power: spm_load will
% load a matrix from an ASCII data file and return the results, and is
% easily used from the ContextMenu. The second example assummes a
% custom funcion called input_cov has been written which expects two
% arguments, for example the following file saved as input_cov.m
% somewhere on the MATLABPATH (~/matlab, the matlab subdirectory of
% your home area, and the current directory, are on the MATLABPATH by
% default):
%
%       function [x] = input_cov(n,decay)
%       % data input routine - mono-exponential covariate
%       % FORMAT [x] = input_cov(n,decay)
%       % n     -  number of time points
%       % decay - decay constant
%       x = exp(-[1:n]/decay);
%
% Although this example is trivial, specifying large vectors of empirical 
% data (e.g. reaction times for 72 scans) is efficient and reliable using 
% this device.
%
% Alternatively, variables containing the required responses may be
% setup in the MatLab base workspace before starting an SPM procedure
% E.g.
% >> tmp=exp(-[1:36]/5.321)
%
% CONDITIONS type input is for getting indicator vectors. The features
% of evaluated input described above are complimented as follows:
%   v)  - a compressed list of digits 0-9   e.g. "12121212"
%  ii)  - a list of indicator characters    e.g. "abababab"
%         (chars mapped to whole numbers in alphabetical order)
%         (case insensitive, [A:Z,a:z] only)
% ...in addition the response is checked to ensure integer condition indices.
%
% Errors in string evaluation are handled gracefully, and the user
% prompted to re-enter.
%
% - BUTTON input -
% For Button input, the prompt is displayed adjacent to a small row of
% buttons. Press the approprate button. The default button (if
% available) has a dark outline. Keyboard accelerators are available:
% <RETURN> or <ENTER> selects the default button (if available). Typing
% the first character of the button label (case insensitive) "presses"
% that button. (If these Keys are not unique, then the integer keys
% 1,2,... "press" the appropriate button.)
%
% The CommandLine variant presents a simple menu of buttons and prompts
% for a selection. Any default response is indicated, and accepted if
% an empty line is input.
%
%
% - MENU input -
% For Menu input, the prompt is displayed in a pull down menu widget.
% Using the mouse, a selection is made by pulling down the widget and
% releasing the mouse on the appropriate response. The default response
% (if set) is marked with an asterisk. Keyboard accelerators are
% available as follows: 'f', 'n' or 'd' move forward to next response
% down; 'b', 'p' or 'u' move backwards to the previous response up the
% list; the number keys jump to the appropriate response number;
% <RETURN> or <ENTER> slelects the currently displayed response. If a
% default is available, then pressing <RETURN> or <ENTER> when the
% prompt is displayed jumps to the default response.
%
% The CommandLine variant presents a simple menu and prompts for a selection.
% Any default response is indicated, and accepted if an empty line is
% input.
%
% - Comand line -
% If YPos is 0 or global CMDLINE is true, then the command line is used.
% Negative YPos overrides CMDLINE, ensuring the GUI is used, at
% YPos=abs(YPos). Similarly relative YPos beginning with '!'
% (E.g.YPos='!+1') ensures the GUI is used.
%
% spm_input uses the SPM 'Interactive' window, which is 'Tag'ged
% 'Interactive'. If there is no such window, then the current figure is
% used, or an 'Interactive' window created if no windows are open.
%
%-----------------------------------------------------------------------
% Programers help is contained in the main body of spm_input.m
%-----------------------------------------------------------------------
% See      : input.m   (MatLab Reference Guide)
% See also : spm_get.m (SPM file selector dialog)
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for programers
%=======================================================================
% - generic   - FORMAT [p,YPos] = spm_input(Prompt,YPos,Type,...)
% - string    - FORMAT [p,YPos] = spm_input(Prompt,YPos,'s',DefStr)
% - evaluated - FORMAT [p,YPos] = spm_input(Prompt,YPos,'e',DefStr,n)
% - condition - FORMAT [p,YPos] = spm_input(Prompt,YPos,'c',DefStr,n)
% - button    - FORMAT [p,YPos]= spm_input(Prompt,YPos,'b',Labels,Values,DefItem)
% - menu      - FORMAT [p,YPos]= spm_input(Prompt,YPos,'m',Labels,Values,DefItem)
% - display   - FORMAT     spm_input(Message,YPos,'d',Label)
% - alert     - FORMAT     spm_input(Alert,YPos,'d!',Label)
%
% - yes/no    - FORMAT p = spm_input(Prompt,YPos,'y/n',Values,DefItem)
% - buttons (shortcut) where Labels is a bar delimited string
%             - FORMAT p = spm_input(Prompt,YPos,Labels,Values,DefItem)
%
% -- Parameters (input) --
%
% Prompt   - prompt string
%          - Defaults (missing or empty) to 'Enter an expression'
%
% YPos     - (numeric) vertical position {1 - 12}
%                          - overriden by global CMDLINE
%                          - 0 for command line
%                          - negative to force GUI
%          - (string) relative vertical position E.g. '+1'
%                          - relative to last position used
%                          - overriden by global CMDLINE
%                          - YPos(1)=='!' forces GUI E.g. '!+1'
%                          - '_' is a shortcut for the lowest GUI position
%          - Defaults (missing or empty) to '+1'
%
% Type     - type of interrogation
%                          - 's'tring
%                          - 'e'valuated string
%                          - 'b'uttons
%                          - 'm'enu pulldown
%                          - 'y/n' : Yes or No buttons
%                                    (See shortcuts below)
%                          - bar delimited string : buttons with these labels
%                                    (See shortcuts below)
%          - Defaults (missing or empty) to 'e'
%
% DefStr   - Default string to be placed in entry widget for string and
%            evaluated types
%          - Defaults to ''
%
% n        - Length of vectors required by 'e' & 'c' types
%            Actually, if *any* dimension is n then input is accepted.
%            This enables input of matrices.
%          - Inf implies no restriction
%          - Default (missing or empty) to Inf - no restriction
%
% Labels   - Labels for button and menu types.
%                          - string matrix, one label per row
%                          - bar delimited string
%                            E.g. 'AnCova|Scaling|None'
%
% Values   - Return values corresponding to Labels for button and menu types
%          - j-th row is returned if button / menu item j is selected
%            (row vectors are transposed)
%          - Defaults (missing or empty) to - (button) Labels
%                                           - ( menu ) menu item numbers
%
% DefItem  - Default item number, for button and menu types.
%
% -- Parameters (output) --
% p        - results
% YPos     - Optional second output argument returns GUI position just used
%
%-----------------------------------------------------------------------
% WINDOWS:
%
% spm_input uses the SPM 'Interactive' 'Tag'ged window. If this isn't
% available and no figures are open, an 'Interactive' SPM window is
% created (`spm('CreateIntWin')`). If figures are available, then the
% current figure is used *unless* it is 'Tag'ged.
%
%-----------------------------------------------------------------------
% SHORTCUTS:
%
% Buttons SHORTCUT - If the Type parameter is a bar delimited string, then
% the Type is taken as 'b' with the specified labels, and the next parameter
% (if specified) is taken for the Values.
%
% Yes/No question shortcut - p = spm_input(Prompt,YPos,'y/n') expands
% to p = spm_input(Prompt,YPos,'b','yes|no',...), enabling easy use of
% spm_input for yes/no dialogue. Values defaults to 'yn', so 'y' or 'n'
% is returned as appropriate.
%
%-----------------------------------------------------------------------
% EXAMPLES:
%            ( Specified YPos is overriden if global CMDLINE is )
%            ( true, when the command line versions are used.   )
%
%       p = spm_input
%               Command line input of an evaluated string, default prompt.
%       p = spm_input('Enter a value',1)
%               Evaluated string input, prompted by 'Enter a value', in
%               position 1 of the dialog figure.
%       p = spm_input(str,'+1','e',0.001)
%               Evaluated string input, prompted by contents of string str,
%               in next position of the dialog figure.
%               Default value of 0.001 offered.
%       p = spm_input(str,2,'e',[],5)
%               Evaluated string input, prompted by contents of string str,
%               in second position of the dialog figure.
%               Vector of length 5 required.
%       p = spm_input(str,0,'c','ababab')
%               Condition string input, prompted by contents of string str
%               Uses command line interface.
%               Default string of 'ababab' offered.
%       p = spm_input(str,0,'c','010101')
%               As above, but default string of '010101' offered.
%       [p,YPos] = spm_input(str,'0','s','Image')
%               String input, same position as last used, prompted by str,
%               default of 'Image' offered. YPos returns GUI position used.
%       p = spm_input(str,'-1','y/n')
%               Yes/No buttons for question with prompt str, in position one 
%               before the last used Returns 'y' or 'n'.
%       p = spm_input(str,'-1','y/n',[1,0],2)
%               As above, but returns yes response or 0 for no, with 'no' as
%               the default response
%       p = spm_input(str,4,'AnCova|Scaling')
%               Presents two buttons labelled 'AnCova' & 'Scaling', with 
%               prompt str, in position 4 of the dialog figure. Returns the 
%               string on the depresed button, where buttons can be pressed 
%               with the mouse or by the respective keyboard accelerators
%               'a' & 's' (or 'A' & 'S').
%       p = spm_input(str,-4,'b','AnCova|Scaling',[],2)
%               As above, but makes "Scaling" the default response, and
%               overrides global CMDLINE
%       p = spm_input(str,0,'b','AnCova|Scaling|None',[1,2,3])
%               Prompts for [A]ncova / [S]caling / [N]one in MatLab command
%               window, returns 1, 2, or 3 according to the first character
%               of the entered string as one of 'a', 's', or 'n' (case 
%               insensitive).
%       p = spm_input(str,'!0','m','Single Subject|Multi Subject|Multi Study')
%               Prints the prompt str in a pull down menu containing items
%               'Single Subject', 'Multi Subject' & 'Multi Study'. When OK is
%               clicked p is returned as the index of the  choice, 1,2, or 3 
%               respectively. Uses last used position in GUI, irrespective of 
%               global CMDLINE
%       p = spm_input(str,5,'m',...
%               'Single Subject|Multi Subject|Multi Study',...
%               ['SS';'MS';'SP'],2)
%               As above, but returns strings 'SS', 'MS', or 'SP' according to
%               the respective choice, with 'MS; as the default response.
%       p = spm_input(str,0,'m',...
%               'Single Subject|Multi Subject|Multi Study',...
%               ['SS';'MS';'SP'],2)
%               As above, but the menu is presented in the command window
%               as a numbered list.
%       spm_input('AnCova, GrandMean scaling',0,'d')
%               Displays message in a box in the MatLab command window
%       [null,YPos]=spm_input('Contrast must sum to zero','+1','d!','Warning')
%		Displays message in next GUI position in red and sounds
%               the keyboard bell. Message is labelled as a 'Warning'
%               Position used is returned in YPos.
%-----------------------------------------------------------------------
% UTILITY FUNCTIONS:
%
% [iCond,msg] = spm_input('!iCond',str,n)
% Parser for special 'c'ondition type: Handles digit strings and
% strings of indicator chars.
% str     - input string
% n       - length of condition vector required [defaut Inf - no restriction]
% iCond   - Integer condition indicator vector
% msg     - status message
%
% [CmdLine,YPos] = spm_input('!CmdLine',YPos)
% Sorts out whether to use CmdLine or not & canonicalises YPos
% CmdLine - Binary flag
% YPos    - Position index
%
% FORMAT Finter = spm_input('!GetWin',F)
% Locates (or creates) figure to work in
% F       - Interactive Figure, defaults to 'Interactive'
% Finter  - Handle of figure to use
%
% FORMAT [PLoc,cF] = spm_input('!PointerJump',RRec,F,XDisp)
% Raise window & jump pointer over question
% RRec  - Response rectangle of current question
% F     - Interactive Figure, Defaults to 'Interactive'
% XDisp - X-displacement of cursor relative to RRec
% PLoc  - Pointer location before jumping
% cF    - Current figure before making F current.
%
% FORMAT [PLoc,cF] = spm_input('!PointerJumpBack',PLoc,cF)
% Replace pointer and reset CurrentFigure back
% PLoc  - Pointer location before jumping
% cF    - Previous current figure
%
% FORMAT spm_input('!PrntPrmpt',Prompt)
% Print prompt for CmdLine questioning
% Prompt - Prompt string
%
% FORMAT [Frec,QRec,PRec,RRec] = spm_input('!InputRects',YPos,rec,F)
% Returns rectangles (pixels) used in GUI
% YPos  - Position index
% rec   - Rectangle specifier: String, one of 'Frec','QRec','PRec','RRec'
%         Defaults to '', which returns them all.
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% FRec  - Position of interactive window
% QRec  - Position of entire question
% PRec  - Position of prompt
% RRec  - Position of response
%
% FORMAT spm_input('!DeleteInputObj',F)
% Deltes input objects (only) from figure F
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
%
% FORMAT [CPos,hCPos] = spm_input('!CurrentPos',F)
% Returns currently used GUI question positions & their handles
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% CPos  - Vector of position indices
% hCPos - (n x CPos) matrix of object handles
%
% FORMAT h = spm_input('!FindInputObj',F)
% Returns handles of input GUI objects in figure F
% F - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
% h - vector of object handles
%
% FORMAT [NPos,CPos,hCPos] = spm_input('!NextPos',YPos,F,CmdLine)
% Returns next position index, specified by YPos
% YPos    - Absolute (integer) or relative (string) position index
%           Defaults to '+1'
% F       - Interactive Figure, defaults to spm_figure('FindWin','Interactive')
% CmdLine - Command line? Defaults to spm_input('!CmdLine',YPos)
% NPos    - Next position index
% CPos & hCPos - as for !CurrentPos
%
% FORMAT NPos = spm_input('!SetNextPos',YPos,F,CmdLine)
% Sets up for input at next position index, specified by YPos. This utility
% function can be used stand-alone to implicitly set the next position
% by clearing positions NPos and greater.
% YPos    - Absolute (integer) or relative (string) position index
%           Defaults to '+1'
% F       - Interactive Figure, defaults to spm_figure('FindWin','Interactive')
% CmdLine - Command line? Defaults to spm_input('!CmdLine',YPos)
% NPos    - Next position index
%
% FORMAT MPos = spm_input('!MaxPos',F,FRec3)
% Returns maximum position index for figure F
% F     - Interactive Figure, Defaults to spm_figure('FindWin','Interactive')
%	  Not required if FRec3 is specified
% FRec3 - Length of interactive figure in pixels
% 
% FORMAT spm_input('!EditableKeyPressFcn',h,ch)
% KeyPress callback for GUI string / eval input
%
% FORMAT spm_input('!ButtonKeyPressFcn',h,Keys,DefItem,ch)
% KeyPress callback for GUI buttons
%
% FORMAT spm_input('!PullDownKeyPressFcn',h,ch,DefItem)
% KeyPress callback for GUI pulldown menus
%
% FORMAT spm_input('!dScroll',h,str)
% Scroll text string in object h
% h      - handle of text object
% Prompt - Text to scroll (Defaults to 'UserData' of h)
%
%-----------------------------------------------------------------------
% SUBFUNCTIONS:
%
% FORMAT [p,msg] = sf_ecEval(str,Type,n)
% Common code for evaluating 'e' & 'c' types.
%_______________________________________________________________________
% Andrew Holmes

%-Parameters
%=======================================================================
COLOR = [.8,.8,1];	%-Question background colour
PJump = 1;		%-Jumping of pointer to question?


%-Condition arguments
%=======================================================================
if nargin<1 | isempty(varargin{1}), Prompt='Enter an expression';
else, Prompt=varargin{1}; end

if Prompt(1)=='!'
	%-Utility functions have Prompt starting with '!'
	Type = Prompt;
else
	%-Condition arguments for question types
	if nargin<6, DefItem=[];  else, DefItem=varargin{6}; end
	if nargin<5, Values=[];   else, Values=varargin{5};  end
	if nargin<4, Labels='';   else, Labels=varargin{4};  end
	if nargin<3 | isempty(varargin{3}), Type='e'; else, Type=varargin{3}; end
	if strcmp(Type,'y/n')
		Type    = 'b';
		DefItem = Values;
		Values  = Labels; if isempty(Values), Values='yn'; end
		Labels  = 'yes|no';
	end
	if any(Type=='|')
		DefItem = Values;
		Values  = Labels;
		Labels  = Type;
		Type    = 'b';
	end
	if nargin<2|isempty(varargin{2}), YPos='+1'; else, YPos=varargin{2}; end
	if iscell(Labels), Labels=char(Labels); end
end


if any(Type(1)=='descbm')
%=======================================================================
%-Setup for GUI use
[CmdLine,YPos] = spm_input('!CmdLine',YPos);
if ~CmdLine
	%-Locate (or create) figure to work in
	Finter = spm_input('!GetWin');

	%-Find out which Y-position to use, setup for use
	YPos = spm_input('!SetNextPos',YPos,Finter,CmdLine);

	%-Determine position of objects
	[FRec,QRec,PRec,RRec] = spm_input('!InputRects',YPos,'',Finter);
end
R2 = YPos;
end % (any(Type(1)=='descbm'))

if Type(1)=='d'
%-Display message
%=======================================================================
if strcmp(Type,'d!'), spm('Beep'), dCol='r'; else, dCol='k'; end
if CmdLine
	fprintf('\n     +-%s%s+',Labels,repmat('-',1,57-length(Labels)))
	Prompt = [Prompt,' '];
	while length(Prompt)>0
		tmp = length(Prompt);
		if tmp>56, tmp=min([max(find(Prompt(1:56)==' ')),56]); end
		fprintf('\n     | %s%s |',Prompt(1:tmp),repmat(' ',1,56-tmp))
		Prompt(1:tmp)=[];
	end
	fprintf('\n     +-%s+\n',repmat('-',1,57))
else
	if ~isempty(Labels), Prompt = [Labels,': ',Prompt]; end
	figure(Finter)
	%-Create text axes and edit control objects
	%---------------------------------------------------------------
	h = uicontrol(Finter,'Style','Text',...
		'String',Prompt(1:min(length(Prompt),56)),...
		'Tag',['GUIinput_',int2str(YPos)],...
		'HorizontalAlignment','Center',...
		'ForegroundColor',dCol,...
		'UserData',Prompt,...
		'Position',QRec);
	if length(Prompt)>56
		pause(1)
		set(h,'ToolTipString',Prompt)
		spm_input('!dScroll',h)
		uicontrol(Finter,'Style','PushButton','String','>',...
			'ToolTipString','press to scroll message',...
			'Tag',['GUIinput_',int2str(YPos)],...
			'UserData',h,...
			'CallBack',[...
			 'set(gcbo,''Visible'',''off''),',...
			 'spm_input(''!dScroll'',get(gcbo,''UserData'')),',...
			 'set(gcbo,''Visible'',''on'')'],...
			'Position',[QRec(1)+QRec(3)-10,QRec(2),15,QRec(4)]);
	end
end
return

elseif any(Type(1)=='esc')
%-String and evaluated string input.
%=======================================================================
DefString = Labels; if ~ischar(DefString), DefString=num2str(DefString); end
n = Values; if isempty(Values), n=Inf; end

if CmdLine
	spm_input('!PrntPrmpt',Prompt)
	if ~isempty(DefString)
		Prompt=[Prompt,' (Default: ',DefString,' )']; end
	str = input([Prompt,' : '],'s'); if isempty(str), str=DefString; end

	%-Eval Types 'e' & 'c' in Base workspace, catch errors
	if any(Type(1)=='ec')
		[p,msg] = sf_ecEval(str,Type,n);
		while isempty(p)
			spm('Beep'), fprintf('! %s : %s\n',mfilename,msg)
			str = input([Prompt,' : '],'s');
			if isempty(str), str=DefString; end
			[p,msg] = sf_ecEval(str,Type,n);
		end
	else
		p = str; msg='';
	end
	if ~isempty(msg), fprintf('\t%s\n',msg), end

else

	%-Create text and edit control objects
	%---------------------------------------------------------------
	hPrmpt = uicontrol(Finter,'Style','Text',...
		'String',Prompt,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'HorizontalAlignment','Right',...
		'Position',PRec);

	%-Edit widget: Callback sets UserData to true when edited
	h = uicontrol(Finter,'Style','Edit',...
		'String',DefString,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'UserData',[],...
		'CallBack','set(gcbo,''UserData'',1)',...
		'Horizontalalignment','Left',...
		'BackgroundColor',COLOR,...
		'Position',RRec);

	%-Figure ContextMenu for shortcuts
	hM = uicontextmenu('Parent',Finter);
	uimenu(hM,'Label','load from text file',...
		'CallBack',[...
		'set(get(gcbo,''UserData''),''UserData'',1,''String'',',...
			'[''spm_load('''''',spm_get(1),'''''')''])'],...
		'UserData',h)
	uimenu(hM,'Label','help on spm_input',...
		'CallBack','spm_help(''spm_input.m'')')
	uimenu(hM,'Label','crash out','Separator','on',...
		'CallBack','delete(get(gcbo,''UserData''))',...
		'UserData',[hPrmpt,h,hM])
	set(Finter,'UIContextMenu',hM)
	
	%-Setup FigureKeyPressFcn to enable editing of entry widget
	% without clicking in it.
	%-If edit widget is clicked in, then it grabs all keyboard input.
	%-Store handle of entry widget in figure UserData
	set(Finter,'KeyPressFcn',[...
	    'spm_input(''!EditableKeyPressFcn'',',...
	    'findobj(gcf,''Tag'',''GUIinput_',int2str(YPos),''',',...
	    	'''Style'',''edit''),',...
	    'get(gcbf,''CurrentCharacter''))'])

	%-Bring window to fore & jump pointer to edit widget
	[PLoc,cF] = spm_input('!PointerJump',RRec,Finter);


	%-Wait for edit, evaluate 'e' & 'c' Types in Base workspace, catch errors
	%---------------------------------------------------------------
	waitfor(h,'UserData')
	if ~ishandle(h), error(['Input window cleared whilst waiting ',...
		'for response: Bailing out!']), end
	str = get(h,'String');
	if any(Type(1)=='ec')
		[p,msg] = sf_ecEval(str,Type,n);
		while isempty(p)
			spm('Beep')
			set(h,'Style','Text',...
				'String',msg,...
				'ForegroundColor','r')
			pause(2)
			set(h,'Style','Edit',...
				'String',str,...
				'Horizontalalignment','Left',...
				'ForegroundColor','k',...
				'UserData',0)
			waitfor(h,'UserData')
			if ~ishandle(h), error(['Input window cleared ',...
				'whilst waiting for response: Bailing out!']),end
			str = get(h,'String');
			[p,msg] = sf_ecEval(str,Type,n);
		end
	else
		p = str; msg='';
	end

	%-Fix edit window, clean up, reposition pointer, set CurrentFig back
	set(h,'Style','Text','HorizontalAlignment','Center',...
		'ToolTipString',msg,...
		'BackgroundColor',[.7,.7,.7]), drawnow
	delete(hM)
	set(Finter,'UserData',[],'KeyPressFcn','')
	spm_input('!PointerJumpBack',PLoc,cF)

end % (if CmdLine)

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	if Type(1)=='e'
		spm_log(['spm_input : ',Prompt,': (',str,')'],p);
	else
		spm_log(['spm_input : ',Prompt,':',p]);
	end
end
varargout = {p,YPos};
return


elseif any(Type(1)=='bm')
%-More complicated input - Types 'b'utton and 'm'enu
%=======================================================================
%-Preliminaries for both types

%-Condition arguments
%-----------------------------------------------------------------------
%-Convert Option string to string matrix if required
if any(Labels=='|')
	OptStr=Labels;
	BarPos=find([OptStr=='|',1]);
	Labels=OptStr(1:BarPos(1)-1);
	for Bar = 2:sum(OptStr=='|')+1
		Labels=str2mat(Labels,OptStr(BarPos(Bar-1)+1:BarPos(Bar)-1));
	end
end

%-Set default Values for the Labels, numeric for 'm'enu Type input
if isempty(Values)
	if Type=='m'
		Values=[1:size(Labels,1)]';
	else
		Values=Labels;
	end
end

%-Make sure Values are in rows
if size(Values,1)==1, Values = Values'; end

%-Check some labels were specified
if isempty(Labels), error('No Labels specified'), end

%-Check numbers of Labels and Values match
if (size(Labels,1)~=size(Values,1))
	error('Labels & Values incompatible sizes'), end

if size(Labels,1)==1
%=======================================================================
%-Only one choice : just return the appropriate Value
	p = Values;


elseif Type(1)=='b'
%=======================================================================
	%-Process button type
	NoLabels = size(Labels,1);
	
	%-Make unique character keys for the Labels, ignoring case
	%---------------------------------------------------------------
	Keys=Labels(:,1)';
	if any(~diff(abs(sort(lower(Keys)))))
		Keys=int2str(1:NoLabels);
		Labs=Labels;
	else
		Labs=Labels;
		Labs(:,1)=[];
	end

	if ~isempty(DefItem)
		DefString = Keys(DefItem);
	else
		DefItem=0;
		DefString = ' ';
	end

	if CmdLine
		%-Display question prompt
		spm_input('!PrntPrmpt',Prompt)
		%-Build prompt
		%-------------------------------------------------------
		if ~isempty(Labs) 
			Prmpt = ['[',Keys(1),']',deblank(Labs(1,:)),' '];
			for lab = 2:NoLabels
			    Prmpt=[Prmpt,...
			    	'/ [',Keys(lab),']',deblank(Labs(lab,:)),' '];
			end
		else
			Prmpt = ['[',Keys(1),'] '];
			for lab = 2:NoLabels
			    Prmpt=[Prmpt,'/ [',Keys(lab),'] ']; end
		end
		if DefItem
		    Prmpt = [Prmpt,...
		    	' (Default: ',deblank(Labels(DefItem,:)),' )'];
		end

		%-Ask for user response
		%-------------------------------------------------------
		str = input([Prmpt,'? '],'s');
		if isempty(str), str=DefString; end
	
		while ~any(lower(Keys)==lower(str(1)))
			fprintf('%c',7)
			str = input([Prmpt,'? '],'s');
			if isempty(str), str=DefString; end
		end
		fprintf('\n')
	
		k = find(lower(Keys)==lower(str(1)));
		p = Values(k,:); if ischar(p), p=deblank(p); end
	
	else

		Tag = ['GUIinput_',int2str(YPos)];

		%-Create text and edit control objects
		%-'UserData' of prompt contains answer
		%-------------------------------------------------------
		hPrmpt = uicontrol(Finter,'Style','Text',...
			'String',Prompt,...
			'Tag',Tag,...
			'UserData',[],...
			'HorizontalAlignment','Right',...
			'Position',PRec);
	
		dX = RRec(3)/NoLabels;
	
		%-Store button # in 'Max' property
		%-Store handle of prompt string in buttons UserData
		%-Button callback sets UserData of prompt string to
		% number of pressed button
		cb = ['set(get(gcbo,''UserData''),''UserData'',',...
			'get(gcbo,''Max''))'];
		H = [];
		XDisp = [];
		for lab=1:NoLabels
			if lab==DefItem
				%-Default button, outline it
				h = uicontrol(Finter,'Style','Frame',...
					'BackGroundColor','k',...
					'ForeGroundColor','k',...
					'Position',...
				[RRec(1)+(lab-1)*dX RRec(2)-1 dX RRec(4)+2]);
				XDisp = (lab-1/3)*dX;
				H = [H,h];
			end
			h = uicontrol(Finter,'Style','Pushbutton',...
				'String',deblank(Labels(lab,:)),...
				'Tag',Tag,...
				'Max',lab,...
				'UserData',hPrmpt,...
				'BackgroundColor',COLOR,...
				'Callback',cb,...
				'Position',...
				[RRec(1)+(lab-1)*dX+1 RRec(2) dX-2 RRec(4)]);
			H = [H,h];
		end
	
		%-Callback for KeyPress, to store valid button # in
		% UserData of Prompt, DefItem if (DefItem~=0)
		% & return (ASCII-13) is pressed
		set(Finter,'KeyPressFcn',...
			['spm_input(''!ButtonKeyPressFcn'',',...
			'findobj(gcf,''Tag'',''',Tag,''',',...
				'''Style'',''text''),',...
			'''',lower(Keys),''',',num2str(DefItem),',',...
			'get(gcbf,''CurrentCharacter''))'])
	
		%-Bring window to fore & jump pointer to default button
		[PLoc,cF] = spm_input('!PointerJump',RRec,Finter,XDisp);

		%-Wait for button press, process results
		%-------------------------------------------------------
		waitfor(hPrmpt,'UserData')
		if ~ishandle(hPrmpt), error(['Input window cleared whilst ',...
			'waiting for response: Bailing out!']), end
		k = get(hPrmpt,'UserData');
		p = Values(k,:); if ischar(p), p=deblank(p); end
		
		%-Display answer and clean up window
		delete(H)
		set(Finter,'KeyPressFcn','')
		uicontrol(Finter,'Style','Text',...
			'String',deblank(Labels(k,:)),...
			'Tag',Tag,...
			'Horizontalalignment','Center',...
			'BackgroundColor',[.7,.7,.7],...
			'Position',RRec);
		spm_input('!PointerJumpBack',PLoc,cF)

	end


elseif Type(1)=='m'
%=======================================================================
	NoLabels = size(Labels,1);
	%-Process pull down menu type
	if CmdLine
		spm_input('!PrntPrmpt',Prompt)
		NoLabels = size(Labels,1);
		for lab = 1:NoLabels
			fprintf('\t%2d : %s\n',lab,Labels(lab,:)), end
		Prmpt = ['Menu choice (1-',int2str(NoLabels),')'];
		if DefItem
		    Prmpt=[Prmpt,' (Default: ',num2str(DefItem),')'];
		end

		k = input([Prmpt,' ? ']);
		if DefItem & isempty(k), k=DefItem; end
		while ~any([1:NoLabels]==k)
			fprintf('%c\t!Out of range',7)
			k = input([Prmpt,' ? ']);
			if DefItem & isempty(k), k=DefItem; end
		end
		fprintf('\n')
	
	else

		Tag = ['GUIinput_',int2str(YPos)];

		Labs=[repmat(' ',NoLabels,2),Labels];
		if DefItem
			Labs(DefItem,1)='*';
			H = uicontrol(Finter,'Style','Frame',...
				'BackGroundColor','k',...
				'ForeGroundColor','k',...
				'Position',QRec+[-1,-1,+2,+2]);
		else
			H = [];
		end
		cb = ['if (get(gcbo,''Value'')>1),',...
			'set(gcbo,''UserData'',''Selected''), end'];
		hPopUp = uicontrol(Finter,'Style','PopUp',...
			'HorizontalAlignment','Left',...
			'ForegroundColor','k',...
			'BackgroundColor',COLOR,...
			'String',str2mat([Prompt,'...'],Labs),...
			'Tag',Tag,...
			'UserData',DefItem,...
			'CallBack',cb,...
			'Position',QRec);

		%-Callback for KeyPresses
		cb=['spm_input(''!PullDownKeyPressFcn'',',...
			'findobj(gcf,''Tag'',''',Tag,'''),',...
			'get(gcf,''CurrentCharacter''))'];
		set(Finter,'KeyPressFcn',cb)

		%-Bring window to fore & jump pointer to menu widget
		[PLoc,cF] = spm_input('!PointerJump',RRec,Finter);

		%-Wait for menu selection
		%-------------------------------------------------------
		waitfor(hPopUp,'UserData')
		if ~ishandle(hPopUp), error(['Input window cleared whilst ',...
			'waiting for response: Bailing out!']), end
		k = get(hPopUp,'Value')-1;
	
		%-Display answer, cleanup window
		delete(H)
		set(Finter,'KeyPressFcn','')
		set(hPopUp,'Style','Text',...
			'Horizontalalignment','Center',...
			'String',deblank(Labels(k,:)),...
			'BackgroundColor',[.7,.7,.7])
		spm_input('!PointerJumpBack',PLoc,cF)
	end

	p = Values(k,:); if ischar(p), p=deblank(p); end

end % (if Type(1)==...) switch for 'b' or 'm' within any(Type(1)=='bm')

%-Set return value
varargout = {p,YPos};

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_input : ',Prompt,': (',deblank(Labels(k,:)),')'],p); end

%-Return
%-----------------------------------------------------------------------
return
end


%=======================================================================
% U T I L I T Y   F U N C T I O N S 
%=======================================================================

switch lower(Type), case '!icond'
%=======================================================================
% [iCond,msg] = spm_input('!iCond',str,n)
% Parse condition spec strings:
%	'2 3 2 3', '0 1 0 1', '2323', '0101', 'abab', 'R A R A'
if nargin<3, n=Inf; else, n=varargin{3}; end
if nargin<2, i=''; else, i=varargin{2}; end
if isempty(i), varargout={[],'empty input'}; return, end
msg = ''; i=i(:)';

if ischar(i)
	if i(1)=='0' & all(ismember(unique(i(:)),setstr(abs('0'):abs('9'))))
		%-Leading zeros in a digit list
		msg = sprintf('%s expanded',i);
		z = min(find(diff(i=='0')));
		i = [zeros(1,z), spm_input('!iCond',i(z+1:end))];
	else
		%-Try an eval, for functions & string #s
		i = evalin('base',['[',i,']'],'i');
	end
end

if ischar(i)
	%-Evaluation error from above: see if it's an 'abab' or 'a b a b' type:
	[c,null,i] = unique(lower(i(~isspace(i))));
	if ~all(ismember(c,setstr(abs('a'):abs('z'))))
		i = [];	msg = 'evaluation error';
	else
		msg = sprintf('[%s] mapped to [1:%d]',c,length(c));
	end
elseif ~all(floor(i(:))==i(:))
	i = []; msg = 'must be integers';
elseif length(i)==1 & n>1
	msg = sprintf('%d expanded',i);
	i = floor(i./10.^[floor(log10(i)+eps):-1:0]);
	i = i-[0,10*i(1:end-1)];
end
if ~isempty(i) & isfinite(n) & length(i)~=n
	i = []; msg = sprintf('%d-vector required',n);
end

varargout = {i,msg};


case '!cmdline'
%=======================================================================
% [CmdLine,YPos] = spm_input('!CmdLine',YPos)
%-Sorts out whether to use CmdLine or not & canonicalises YPos
if nargin<2, YPos=''; else, YPos=varargin{2}; end
if isempty(YPos), YPos='+1'; end

%-Global CMDLINE
CmdLine = spm('isGCmdLine');

%-Special YPos specifications
if ischar(YPos)
	if(YPos(1)=='!'), CmdLine=0; YPos(1)=[]; end
elseif YPos==0
	CmdLine=1;
elseif YPos<0
	CmdLine=0;
	YPos=-YPos;
end
varargout = {CmdLine,YPos};


case '!getwin'
%=======================================================================
% Finter = spm_input('!GetWin',F)
%-Locate (or create) figure to work in (Don't use 'Tag'ged figs)
if nargin<2, F='Interactive'; else, F=varargin{2}; end
Finter = spm_figure('FindWin',F);
if isempty(Finter)
	if any(get(0,'Children'))
		if isempty(get(gcf,'Tag')), Finter = gcf;
		else, Finter = spm('CreateIntWin'); end
	else, Finter = spm('CreateIntWin'); end
end
varargout = {Finter};


case '!pointerjump'
%=======================================================================
% [PLoc,cF] = spm_input('!PointerJump',RRec,F,XDisp)
%-Raise window & jump pointer over question
if nargin<4, XDisp=[]; else, XDisp=varargin{4}; end
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, error('Insufficient arguments'), else, RRec=varargin{2}; end
F = spm_figure('FindWin',F);
PLoc = get(0,'PointerLocation');
cF   = get(0,'CurrentFigure');
if ~isempty(F)
	figure(F)
	FRec = get(F,'Position');
	if isempty(XDisp), XDisp=RRec(3)*4/5; end
	if PJump, set(0,'PointerLocation',...
		floor([(FRec(1)+RRec(1)+XDisp), (FRec(2)+RRec(2)+RRec(4)/3)]));
	end
end
varargout = {PLoc,cF};


case '!pointerjumpback'
%=======================================================================
% spm_input('!PointerJumpBack',PLoc,cF)
%-Replace pointer and reset CurrentFigure back
if nargin<4, cF=[]; else, F=varargin{3}; end
if nargin<2, error('Insufficient arguments'), else, PLoc=varargin{2}; end
if PJump, set(0,'PointerLocation',PLoc), end
cF = spm_figure('FindWin',cF);
if ~isempty(cF), set(0,'CurrentFigure',cF); end
return


case '!prntprmpt'
%=======================================================================
% spm_input('!PrntPrmpt',Prompt)
%-Print prompt for CmdLine questioning
if nargin<2, Prompt=''; else, Prompt=varargin{2}; end
if isempty(Prompt), Prompt='Enter an expression'; end
fprintf('\n%s\n\t%s\n%s\n',repmat('=',1,70),Prompt,repmat('=',1,70))
return


case '!inputrects'
%=======================================================================
% [Frec,QRec,PRec,RRec,Sz,Se] = spm_input('!InputRects',YPos,rec,F)
if nargin<4, F='Interactive'; else, F=varargin{4}; end
if nargin<3, rec=''; else, rec=varargin{3}; end
if nargin<2, YPos=1; else, YPos=varargin{2}; end
F = spm_figure('FindWin',F);
if isempty(F), error('Figure not found'), end

Units = get(F,'Units');
set(F,'Units','pixels')
FRec = get(F,'Position');
set(F,'Units',Units);
Xdim = FRec(3); Ydim = FRec(4);

WS   = spm('GetWinScale');
Sz   = 2*round(22*min(WS)/2);	%-Height
Pd   = Sz/2;			%-Pad
Se   = 2*round(25*min(WS)/2);	%-Seperation
Yo   = round(2*min(WS));	%-Y offset for responses

a = 5.5/10;
y = Ydim - Se*YPos;
QRec   = [Pd            y         Xdim-2*Pd        Sz];	%-Question
PRec   = [Pd            y     floor(a*Xdim)-2*Pd   Sz];	%-Prompt
RRec   = [ceil(a*Xdim)  y+Yo  floor((1-a)*Xdim)-Pd Sz];	%-Response
% MRec = [010           y         Xdim-50          Sz];	%-Menu PullDown
% BRec = MRec + [Xdim-50+1, 0+1, 50-Xdim+30, 0];	%-Menu PullDown OK butt.

if ~isempty(rec)
	varargout = {eval(rec)};
else
	varargout = {FRec,QRec,PRec,RRec,Sz,Se};
end


case '!deleteinputobj'
%=======================================================================
% spm_input('!DeleteInputObj',F)
if nargin<2, F='Interactive'; else, F=varargin{2}; end
h = spm_input('!FindInputObj',F);
delete(h(h>0))


case {'!currentpos','!findinputobj'}
%=======================================================================
% [CPos,hCPos] = spm_input('!CurrentPos',F)
% h            = spm_input('!FindInputObj',F)
% hPos contains handles: Columns contain handles corresponding to Pos
if nargin<2, F='Interactive'; else, F=varargin{2}; end
F = spm_figure('FindWin',F);

%-Find tags and YPos positions of 'GUIinput_' 'Tag'ged objects
H    = [];
YPos = [];
for h = get(F,'Children')'
	tmp = get(h,'Tag');
	if ~isempty(tmp)
		if strcmp(tmp(1:min(length(tmp),9)),'GUIinput_')
			H    = [H, h];
			YPos = [YPos, eval(tmp(10:end))];
		end
	end
end

switch lower(Type), case '!findinputobj'
	varargout = {H};
case '!currentpos'
	if nargout<2
		varargout = {max(YPos),[]};
	elseif isempty(H)
		varargout = {[],[]};
	else
		%-Sort out 
		tmp     = sort(YPos);
		CPos    = tmp(find([1,diff(tmp)]));
		nPos    = length(CPos);
		nPerPos = diff(find([1,diff(tmp),1]));
		hCPos   = zeros(max(nPerPos),nPos);
		for i = 1:nPos
			hCPos(1:nPerPos(i),i) = H(YPos==CPos(i))';
		end
		varargout = {CPos,hCPos};
	end
end

case '!nextpos'
%=======================================================================
% [NPos,CPos,hCPos] = spm_input('!NextPos',YPos,F,CmdLine)
%-Return next position to use
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, YPos='+1'; else, YPos=varargin{2}; end
if nargin<4, [CmdLine,YPos]=spm_input('!CmdLine',YPos);
	else, CmdLine=varargin{4}; end

F = spm_figure('FindWin',F);

%-Get current positions
if nargout<3
	CPos = spm_input('!CurrentPos',F);
	hCPos = [];
else
	[CPos,hCPos] = spm_input('!CurrentPos',F);
end

if CmdLine
	NPos = 0;
else
	MPos = spm_input('!MaxPos',F);
	if ischar(YPos)
		%-Relative YPos
		%-Strip any '!' prefix from YPos
		if(YPos(1)=='!'), YPos(1)=[]; end
		if strcmp(YPos,'_')
			%-YPos='_' means bottom
			NPos=MPos;
		else
			YPos = max([0,CPos])+eval(YPos);
		end
	else
		%-Absolute YPos
		YPos=abs(YPos);
	end
	NPos = min(max(1,YPos),MPos);
end
varargout = {NPos,CPos,hCPos};

case '!setnextpos'
%=======================================================================
% NPos = spm_input('!SetNextPos',YPos,F,CmdLine)
%-Set next position to use
if nargin<3, F='Interactive'; else, F=varargin{3}; end
if nargin<2, YPos='+1'; else, YPos=varargin{2}; end
if nargin<4, [CmdLine,YPos]=spm_input('!CmdLine',YPos);
	else, CmdLine=varargin{4}; end

%-Find out which Y-position to use
[NPos,CPos,hCPos] = spm_input('!NextPos',YPos,F,CmdLine);

%-Delete any previous inputs using positions NPos and after
if any(CPos>=NPos), h=hCPos(:,CPos>=NPos); delete(h(h>0)), end

varargout = {NPos};

case '!maxpos'
%=======================================================================
% MPos = spm_input('!MaxPos',F,FRec3)
%
if nargin<3
	if nargin<2, F='Interactive'; else, F=varargin{2}; end
	F = spm_figure('FindWin',F);
	if isempty(F)
		FRec3=spm('WinSize','Interactive')*[0;0;0;1];
	else
		%-Get figure size
		Units = get(F,'Units');
		set(F,'Units','pixels')
		FRec3 = get(F,'Position')*[0;0;1;0];
		set(F,'Units',Units);
	end
end

Se   = 2*round(25*min(spm('GetWinScale'))/2);
MPos = floor((FRec3-5)/Se);

varargout = {MPos};


case '!editablekeypressfcn'
%=======================================================================
% spm_input('!EditableKeyPressFcn',h,ch)
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcbf,'KeyPressFcn','','UserData',[]), return, end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=varargin{3};end

tmp = get(h,'String');

if isempty(ch)
	%- shift / control / &c. pressed
	return
elseif any(abs(ch)==[32:126])
	tmp = [tmp, ch];
elseif abs(ch)==21
	%- ^U - kill
	tmp = '';
elseif any(abs(ch)==[8,127])
	%-BackSpace or Delete
	if length(tmp), tmp(length(tmp))=''; end
elseif abs(ch)==13
	if ~isempty(tmp), set(h,'UserData',1), end
else
	%-Illegal character
	return
end
set(h,'String',tmp)


case '!buttonkeypressfcn'
%=======================================================================
% spm_input('!ButtonKeyPressFcn',h,Keys,DefItem,ch)
%-Callback for KeyPress, to store valid button # in UserData of Prompt,
% DefItem if (DefItem~=0) & return (ASCII-13) is pressed

%-Condition arguments
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcf,'KeyPressFcn','','UserData',[]), return, end
if nargin<3, error('Insufficient arguments'); else, Keys=varargin{3}; end
if nargin<4, DefItem=0; else, DefItem=varargin{4}; end
if nargin<5, ch=get(gcf,'CurrentCharacter'); else, ch=varargin{5}; end

if isempty(ch)
	%- shift / control / &c. pressed
	return
elseif (DefItem & ch==13)
	But = DefItem;
else
	But = find(ch==Keys);
end
if ~isempty(But), set(h,'UserData',But), end


case '!pulldownkeypressfcn'
%=======================================================================
% spm_input('!PullDownKeyPressFcn',h,ch,DefItem)
if nargin<2, error('Insufficient arguments'), else, h=varargin{2}; end
if isempty(h), set(gcf,'KeyPressFcn',''), return, end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=varargin{3};end
if nargin<4, DefItem=get(h,'UserData'); else, ch=varargin{4}; end

Pmax = get(h,'Max');
Pval = get(h,'Value');

if Pmax==1, return, end

if isempty(ch)
	%- shift / control / &c. pressed
	return
elseif abs(ch)==13
	if Pval==1
		if DefItem, set(h,'Value',max(2,min(DefItem+1,Pmax))), end
	else
		set(h,'UserData','Selected')
	end
elseif any(ch=='bpu')
	%-Move "b"ack "u"p to "p"revious entry
	set(h,'Value',max(2,Pval-1))
elseif any(ch=='fnd')
	%-Move "f"orward "d"own to "n"ext entry
	set(h,'Value',min(Pval+1,Pmax))
elseif any(ch=='123456789')
	%-Move to entry n
	set(h,'Value',max(2,min(eval(ch)+1,Pmax)))
else
	%-Illegal character
end


case '!dscroll'
%=======================================================================
% spm_input('!dScroll',h,Prompt)
%-Scroll text in object h
if nargin<2, return, else, h=varargin{2}; end
if nargin<3, Prompt = get(h,'UserData'); else, Prompt=varargin{3}; end
tmp = Prompt;
if length(Prompt)>56
	while length(tmp)>56
		tic, while(toc<0.1), pause(0.05), end
		tmp(1)=[];
		set(h,'String',tmp(1:min(length(tmp),56)))
	end
	pause(1)
	set(h,'String',Prompt(1:min(length(Prompt),56)))
end


otherwise
%=======================================================================
warning('Invalid type / action')

%=======================================================================
end


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function [p,msg] = sf_ecEval(str,Type,n)
%-----------------------------------------------------------------------
if nargin<3, n=Inf; end
if nargin<2, Type='e'; end
if nargin<1, str=''; end
if isempty(str), p=[]; msg='empty input'; return, end
switch Type(1)
case 'e'
	p = evalin('base',['[',str,']'],'[]');
	if isempty(p), msg = 'evaluation error'; else, msg=''; end
case 'c'
	[p,msg] = spm_input('!iCond',str);
end
if ~isempty(p) & isfinite(n) & ~any(size(p)==n)
	p = []; msg = sprintf('%d-vector required',n);
end
