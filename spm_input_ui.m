function [R1,R2,R3,R4] = spm_input(P1,P2,P3,P4,P5,P6)
% Comprehensive graphical and command line input function
% FORMAT p = spm_input(Prompt,YPos,Type,Labels,Values,Defaults)
%_______________________________________________________________________
%
% spm_input handles most forms of user input for SPM. (File selection
% is handled by spm_get.m) This can be either via GUI with keyboard
% shortcuts (in the SPM "Interactive" window), or on the command line
% (in the MatLab command window) if requested.
%
% There are four types of input: String, Evaluated, Buttons and Menus:
% These prompt for string input; string input which is evaluated to
% give a numerical result; selection of one item from a set of buttons;
% selection of an item from a menu.
%
% - STRING & EVALUATED input -
% For String and Evaluated input, a prompt is displayed adjacent to an
% editable text entry widget (with a lilac background!). This entry
% widget will initially contain the defualt reply, if specified.
% Pressing <RETURN> or <ENTER> will accept the default. Clicking in the
% entry widget allows editing, pressing <RETURN> or <ENTER> enters the
% result. You must enter something, empty answers are not accepted.
% (Note that once you start editing you must make a change before
% <RETURN> or <ENTER> will submit your answer - a quirk of MatLab I'm
% afraid.)
%
% The entry widget can be pasted into in the usual manner for your
% platform. On UNIX, this is a Motif widget, so select the text you
% wish to enter and use the middle mouse button to paste the it into
% the entry widget.
%
% Basic editing of the entry widget is supported *without* clicking in
% the widget. This enables you to type responses to a sequence of
% questions without having to repeatedly click the mouse. Supported are
% BackSpace and Delete, line kill (^U). Other standard ASCII characters
% are appended to the text in the entry widget. Press <RETURN> or <ENTER>
% to submit your response.
%
% For evaluated input, the string submitted is evaluated in the base
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
% load a matrix from an ASCII data file and return the results. The
% second example assummes a custom funcion called input_cov has been
% written which expects two arguments, for example the following file
% saved as input_cov.m somewhere on the MATLABPATH (~/matlab, the
% matlab subdirectory of your home area, and the current directory, are
% on the MATLABPATH by default):
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
% setup in the MatLab base workspace before invoking the SPM routine.
% E.g.
% >> tmp=exp(-[1:36]/5.321)
%
% Errors in string evaluation are handled gracefully, and the user
% prompted to re-enter.
%
% In the CommandLine variant evaluated input is evaluated in the
% spm_input's workspace, so base workspace variables are not available.
% Any default response is indicated, and accepted if an empty line is
% input.
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
% ( - Technical note - : This form of interactive GUI forces the code to )
% ( wait for a response in a computer intensive loop. To avoid           )
% ( hammering the CPU when queries are left unattended, spm_input has    )
% ( "sleep" modes. Wake spm_input by clicking in any MatLab window. For  )
% ( GUI menu selections using spm_input, you have 120s to make a         )
% ( selection before sleep mode sets in. For GUI button selections, you  )
% ( have 10s to release the chosen button, after which a further GUI     )
% ( event will accept the selection.                                     )
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
% - generic    - FORMAT p = spm_input(Prompt,YPos,Type,...)
% - string     - FORMAT p = spm_input(Prompt,YPos,'s',DefStr)
% - evaluated  - FORMAT p = spm_input(Prompt,YPos,'e',DefStr)
% - button     - FORMAT p = spm_input(Prompt,YPos,'b',Labels,Values,DefItem)
% - menu       - FORMAT p = spm_input(Prompt,YPos,'m',Labels,Values,DefItem)
%
% - yes/no     - FORMAT p = spm_input(Prompt,YPos,'y/n',Values,DefItem)
% - buttons (shortcut) where Labels is a bar delimited string
%              - FORMAT p = spm_input(Prompt,YPos,Labels,Values,DefItem)
%
% - Parameters -
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
% p        - results
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
% EXAMPLES :
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
%       p = spm_input(str,'0','s')
%               String input, same position as last used, prompted by str,
%               default of 'Image' offered
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
%_______________________________________________________________________
% Andrew Holmes

%-Parameters
%=======================================================================
COLOR = [.8,.8,1];	%-Question background colour


%-Condition arguments
%=======================================================================
if nargin<1, Prompt=''; else, Prompt=P1; end
if isempty(Prompt), Prompt='Enter an expression'; end

if Prompt(1)=='!'
	%-Utility functions have Prompt starting with '!'
	Type = Prompt;
else
	%-Condition arguments for question types
	if nargin<6, DefItem=[];  else, DefItem=P6; end
	if nargin<5, Values=[]; else, Values=P5; end
	if nargin<4, Labels=[]; else, Labels=P4; end
	if nargin<3, Type='';   else, Type=P3; end
	if isempty(Type), Type='e'; end
	if strcmp(Type,'y/n')
		Type = 'b';
		DefItem=Values;
		Values=Labels; if isempty(Values), Values='yn'; end
		Labels='yes|no';
	end
	if any(Type=='|')
		DefItem=Values;
		Values=Labels;
		Labels=Type;
		Type='b';
	end
	if nargin<2, YPos=[];   else, YPos=P2; end
	if isempty(YPos), YPos='+1'; end
end


if any(Type(1)=='esbm')
%=======================================================================
% Preliminaries for all question types

%-Command line? -> If global CMDLINE is true & YPos > 0, or if YPos is empty
%-----------------------------------------------------------------------
CmdLine = spm('isGCmdLine');
if isstr(YPos), if YPos(1)=='!', CmdLine=0; YPos(1)=[]; end
elseif YPos==0, CmdLine = 1;
elseif YPos <0, CmdLine = 0; YPos=abs(YPos);
end

%-Setup
%-----------------------------------------------------------------------
if CmdLine
	%-Print out prompt to command window
	%---------------------------------------------------------------
	fprintf('\n'), fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
	fprintf('\t%s\n',Prompt)
	fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
else
	%-Locate (or open) figure to work in
	%---------------------------------------------------------------
	Finter = spm_figure('FindWin','Interactive');
	if isempty(Finter), if any(get(0,'Children')), Finter = gcf;
				else, Finter = spm('CreateIntWin'); end
	end

	%-Find out which Y-position to use - CmdLine over-rides
	%---------------------------------------------------------------
	[CPos,hCPos] = spm_input('!CurrentPos',Finter);
	if isstr(YPos)
		YPos = max(1,(max([0,max(CPos)])+eval(YPos)));
	end

	%-Delete any previous inputs using positions YPos and after
	%---------------------------------------------------------------
	% delete(findobj(Finter,'Tag',['GUIinput_',int2str(YPos)]))
	if any(CPos>=YPos)
		h = hCPos(:,CPos>=YPos);
		delete(h(h>0))
	end

	%-Determine position of objects
	%---------------------------------------------------------------
	[FRec,QRec,PRec,RRec] = spm_input('!InputRects',YPos,'',Finter);
	
	%-Bring window to fore & place pointer over control object
	%---------------------------------------------------------------
	figure(Finter)
	CPos = get(0,'PointerLocation');
	set(0,'PointerLocation',...
		[(FRec(1) +RRec(1) +RRec(3)/5),...
			(FRec(2) +RRec(2) +RRec(4)/2)]);
end % (if CmdLine)
end % (any(Type(1)=='esbm'))



if (Type(1)=='e')|(Type(1)=='s')
%-String and evaluated string input.
%=======================================================================
DefString = Labels; if ~isstr(DefString), DefString=num2str(DefString); end

if CmdLine
	if ~isempty(DefString)
		Prompt=[Prompt,' (Default: ',DefString,' )']; end
	str = input([Prompt,' : '],'s'); if isempty(str), str=DefString; end
	if Type(1)=='e', p=eval(['[',str,']'],'sprintf(''<ERROR>'')');
		else, p=str; end

	%-Catch eval errors for Type 'e', make sure 
	%-(Input can do this by itself, but doesn't return the typed string)
	while (strcmp(p,'<ERROR>') & Type(1)=='e') | isempty(p)
		fprintf('%c! spm_input : ',7)
		if isempty(p), fprintf('enter something!\n')
			else,  fprintf('evaluation error\n'), end
		str = input([Prompt,' : '],'s');
		if isempty(str), str=DefString; end
		if Type(1)=='e', p=eval(['[',str,']'],'sprintf(''<ERROR>'')');
			else, p=str; end
	end % (while)

else

	%-Create text and edit control objects
	%---------------------------------------------------------------
	uicontrol(Finter,'Style','Text',...
		'String',Prompt,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'HorizontalAlignment','Right',...
		'Position',PRec);
	
	%-Construct callback for edit widget - if evaluating input,
	% do it in the callback, so that base workspace variables can be used
	if Type(1)=='e'
		cb = ['set(gco,''UserData'',',...
			'eval([''['',get(gco,''String''),'']''],',...
				' ''sprintf(''''<ERROR>'''')''))'];
	else
		cb = 'set(gco,''UserData'',get(gco,''String'') )';
	end

	%-Edit widget: Callback sets UserData to result when edited
	h = uicontrol(Finter,'Style','Edit',...
		'String',DefString,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'UserData',[],...
		'CallBack',cb,...
		'Horizontalalignment','Left',...
		'BackgroundColor',COLOR,...
		'Position',RRec);

	%-Setup FigureKeyPressFcn to enable editing of entry widget
	% without clicking in it. Return for entry handled in actual
	% callback string (rather than callback function) so that evals
	% are made in the base workspace where users may have setup
	% variables.
	%-If edit widget is clicked in, then it grabs all keyboard input
	%-Store handle of entry widget in figure UserData
	if Type(1)=='e'
	    cb = ['set(get(gcf,''UserData''),''UserData'',',...
	    'eval([''['',get(get(gcf,''UserData''),''String''),'']''],',...
			' ''sprintf(''''<ERROR>'''')''))'];
	else
	    cb = ['set(get(gcf,''UserData''),''UserData'',',...
			'get(get(gcf,''UserData''),''String'') )'];
	end
	set(Finter,'UserData',h)
	set(Finter,'KeyPressFcn',[...
		'if abs(get(gcf,''CurrentCharacter''))==13,',...
			cb,',',...
		'else,',...
			'spm_input(''!EditableKeyPressFcn'',',...
				'get(gcf,''UserData''),',...
				'get(gcf,''CurrentCharacter'')),',...
		'end'])

	%-Wait for edit, evaluate string if input not a string variable
	%---------------------------------------------------------------
	while isempty(get(h,'UserData'))
		waitforbuttonpress; pause(0.1), pause(0.1), end

	p = get(h,'UserData');
	if Type(1)=='e'
		while strcmp(p,'<ERROR>')	%-Catch eval errors
			fprintf('%c',7)
			set(h,'UserData',[])
			tmp = uicontrol(Finter,'Style','Text',...
				'String','<ERROR>',...
				'Tag',['GUIinput_',int2str(YPos)],...
				'Horizontalalignment','Center',...
				'BackgroundColor',COLOR,...
				'ForegroundColor','r',...
				'Position',RRec);
			drawnow, pause(2)
			delete(tmp)
			while isempty(get(h,'UserData'));
			    waitforbuttonpress; pause(0.1), pause(0.1), end
			p = get(h,'UserData');
		end % (while)
	end

	%-Fix edit window, clean up, reposition pointer
	set(h,'Style','Text','HorizontalAlignment','Center',...
		'BackgroundColor',[.7,.7,.7]), drawnow
	set(Finter,'UserData',[],'KeyPressFcn','')
	set(0,'PointerLocation',CPos)

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
R1 = p;
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
		end % (while)
		fprintf('\n')
	
		k = find(lower(Keys)==lower(str(1)));
		p = Values(k,:); if isstr(p), p=deblank(p); end
	
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
	
		%-Tag buttons with their button number
		%-Store handle of prompt string in buttons UserData
		%-Button callback sets UserData of prompt string to
		% number of pressed button
		cb = ['set(get(gco,''UserData''),''UserData'',',...
			'eval(get(gco,''Tag'')))'];
		H = [];
		for lab=1:NoLabels
			if lab==DefItem
				%-Default button, outline it
				h = uicontrol(Finter,'Style','Frame',...
					'BackGroundColor','k',...
					'ForeGroundColor','k',...
					'Position',...
				[RRec(1)+(lab-1)*dX RRec(2)-1 dX 022]);
				H = [H,h];
			end
			h = uicontrol(Finter,'Style','Pushbutton',...
				'String',deblank(Labels(lab,:)),...
				'Tag',num2str(lab),...
				'UserData',hPrmpt,...
				'BackgroundColor',COLOR,...
				'Callback',cb,...
				'Position',...
				[RRec(1)+(lab-1)*dX+1 RRec(2) dX-2 020]);
			H = [H,h];
		end
	
		%-Callback for KeyPress, to store valid button # in
		% UserData of Prompt, 0 if dmDel & return (ASCII-13) 
		% is pressed
		if DefItem
			cb=['if abs(get(gcf,''CurrentCharacter''))==13,',...
				'set(findobj(gcf,''Tag'',''',...
				Tag,'''),''UserData'',0),',...
			'elseif any(lower(get(gcf,''CurrentCharacter''))',...
			        '==''',lower(Keys),'''),',...
			        'set(findobj(gcf,''Tag'',''',...
			        Tag,'''),''UserData'',',...
			        'find(lower(get(gcf,',...
			        '''CurrentCharacter''))==''',...
			        lower(Keys),''')),',...
		    	'end'];
		else
			cb=['if any(lower(get(gcf,''CurrentCharacter''))',...
			    '==''',lower(Keys),'''),',...
			    'set(findobj(gcf,''Tag'',''',...
			    Tag,'''),''UserData'',',...
			    'find(lower(get(gcf,',...
			    '''CurrentCharacter''))==''',...
			    lower(Keys),''')),end'];
		end
		set(Finter,'KeyPressFcn',cb)
	
		%-Wait for button press, process results
		%-------------------------------------------------------
		while isempty(get(hPrmpt,'UserData'))
			waitforbuttonpress;
			%-Give 10s grace for release of button
			%-After this time spm_input "sleeps", and the user
			% will have to press another button to have choice
			% registered
			% ( This loop hammers the CPU, but there's no other )
			% ( way to do an interruptible pause without        ) 
			% ( requiring a second buttonpress event! NB: pause )
			% ( (ML4.2c) only works for integer seconds!        )
			tic, while(toc<10)
				if ~isempty(get(hPrmpt,'UserData'))
					break, end
				pause(0.25)
			end
		end
		k = get(hPrmpt,'UserData');
		%-Sort out default (if specified)
		if k==0; k=DefItem; end
	
		p = Values(k,:); if isstr(p), p=deblank(p); end
		
		%-Display answer and clean up window
		delete(H)
		set(Finter,'KeyPressFcn','')
		uicontrol(Finter,'Style','Text',...
			'String',deblank(Labels(k,:)),...
			'Tag',Tag,...
			'Horizontalalignment','Center',...
			'BackgroundColor',[.7,.7,.7],...
			'Position',RRec);
		set(0,'PointerLocation',CPos)

	end % (if CmdLine)


elseif Type(1)=='m'
%=======================================================================
	NoLabels = size(Labels,1);
	%-Process pull down menu type
	if CmdLine
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

		Labs=[setstr(ones(NoLabels,2)*' '),Labels];
		if DefItem
			Labs(DefItem,1)='*';
			cb = ['if (get(gco,''Value'')>1),',...
				'set(gco,''UserData'',''Selected''), end'];
		else
			cb = ['if (get(gco,''Value'')>1),',...
				'set(gco,''UserData'',''Selected''), end'];
		end
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

		%-Wait for menu selection
		%-------------------------------------------------------
		while ~isstr(get(hPopUp,'UserData'))
			%-Give approx 120s grace for menu selection
			%-After this time spm_input "sleeps", waiting for
			% a buttonpress to reactivate.
			% ( This loop hammers the CPU, but there's no other )
			% ( way to do an interruptible pause for PullDowns. )
			% ( NB: pause (ML4.2c) only works for integer       )
			% ( seconds!                                        )
			tic, while(toc<120)
				if isstr(get(hPopUp,'UserData')), break, end
				pause(0.25)
			end
			if ~isstr(get(hPopUp,'UserData'))
				waitforbuttonpress; end
		end
		
		k = get(hPopUp,'Value')-1;
	
		%-Display answer, cleanup window
		set(Finter,'KeyPressFcn','')
		set(hPopUp,'Style','Text',...
			'Horizontalalignment','Center',...
			'String',deblank(Labels(k,:)),...
			'BackgroundColor',[.7,.7,.7])
		set(0,'PointerLocation',CPos)
	end % (if CmdLine)

	p = Values(k,:); if isstr(p), p=deblank(p); end

end % (if Type(1)==...) switch for 'b' or 'm' within any(Type(1)=='bm')

%-Set return value
R1 = p;

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_input : ',Prompt,': (',deblank(Labels(k,:)),')'],p); end

%-Return
%-----------------------------------------------------------------------
return


%=======================================================================
% U T I L I T Y   F U N C T I O N S 
%=======================================================================

elseif strcmp(lower(Type),lower('!InputRects'))
%=======================================================================
% [Frec,QRec,PRec,RRec] = spm_input('!InputRects',YPos,rec,F)
if nargin<4, F='Interactive'; else, F=P4; end
if nargin<3, rec=''; else, rec=P3; end
if nargin<2, YPos=1; else, YPos=P2; end
F = spm_figure('FindWin',F);
if isempty(F), error('Figure not found'), end

Units = get(F,'Units');
set(F,'Units','pixels')
FRec = get(F,'Position');
set(F,'Units',Units);
Xdim = FRec(3); Ydim = FRec(4);

a = 5.5/10;
y = Ydim - 25*YPos;
QRec = [010     y      Xdim-20      020];
PRec = [010     y  (a*Xdim -20)     020];
RRec = [a*Xdim  y  ((1-a)*Xdim -10) 020];
MRec = [010     y      Xdim-50      020];
BRec = MRec + [Xdim-50+1, 0+1, 50-Xdim+30, 0];

if ~isempty(rec)
	R1 = eval(rec);
else
	R1 = FRec;
	R2 = QRec;
	R3 = PRec;
	R4 = RRec;
	R5 = MRec;
	R6 = BRec;
end
return



elseif strcmp(lower(Type),lower('!CurrentPos'))
%=======================================================================
% [CPos,hCPos] = spm_input('!CurrentPos',F)
% hPos contains handles: Columns contain handles corresponding to Pos
if nargin<2, F='Interactive'; else, F=P2; end
F = spm_figure('FindWin',F);
if isempty(F), return, end

%-Find tags and YPos positions of 'GUIinput_' 'Tag'ged objects
H    = [];
YPos = [];
for h = get(F,'Children')'
	tmp = get(h,'Tag');
	if ~isempty(tmp)
		if strcmp(tmp(1:min(length(tmp),9)),'GUIinput_')
			H    = [H, h];
			YPos = [YPos, eval(tmp(10:length(tmp)))];
		end
	end
end

if nargout<2
	CPos  = max([YPos,0]);
	hCPos = [];
elseif isempty(H)
	CPos  = [];
	hCPos = [];
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
end
R1 = CPos;
R2 = hCPos;
return

elseif strcmp(lower(Type),lower('!EditableKeyPressFcn'))
%=======================================================================
% spm_input('!EditableKeyPressFcn',h,ch)
if nargin<2, error('Insufficient arguments'), else, h=P2; end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=P3; end

tmp = get(h,'String');

if any(abs(ch)==[32:126])
	tmp = [tmp, ch];
elseif abs(ch)==21
	%- ^U - kill
	tmp = '';
elseif any(abs(ch)==[8,127])
	%-BackSpace or Delete
	if length(tmp), tmp(length(tmp))=''; end
elseif abs(ch)==13
	fprintf('Return shouldn''t be passed to this routine!')
else
	%-Illegal character
	return
end
set(h,'String',tmp)
return


elseif strcmp(lower(Type),lower('!PullDownKeyPressFcn'))
%=======================================================================
% spm_input('!PullDownKeyPressFcn',h,ch,DefItem)
if nargin<2, error('Insufficient arguments'), else, h=P2; end
if nargin<3, ch=get(get(h,'Parent'),'CurrentCharacter'); else, ch=P3; end
if nargin<4, DefItem=get(h,'UserData'); else, ch=P4; end

Pmax = get(h,'Max');
Pval = get(h,'Value');

if Pmax==1, return, end

if abs(ch)==13
	if Pval==1
		if DefItem, set(h,'Value',max(2,min(DefItem+1,Pmax))), end
	else
		set(h,'UserData','Selected')
	end
elseif ch=='b' | ch=='p' | ch=='u'
	%-Move "b"ack "u"p to "p"revious entry
	set(h,'Value',max(2,Pval-1))
elseif ch=='f' | ch=='n' | ch=='d'
	%-Move "f"orward "d"own to "n"ext entry
	set(h,'Value',min(Pval+1,Pmax))
elseif any(ch=='123456789')
	%-Move to entry n
	set(h,'Value',max(2,min(eval(ch)+1,Pmax)))
else
	%-Illegal character
end
return


else
%=======================================================================
error('Invalid Type')

end
