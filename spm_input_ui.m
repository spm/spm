function p = spm_input(Prompt,YPos,Type,Labels,Values)
% Comprehensive graphical and command line input function
%
% FORMAT p = spm_input(Prompt,YPos,Type,Labels,Values)
%
% Prompt   - prompt string
% YPos     - vertical position {1 - 12}
%                                  - overriden by global CMDLINE
%                                  - 0 for command line
%                                  - negative to force GUI
% Type     - type of interrogation - 's'tring, 'e'valuated string
%                                    'b'uttons, pull down 'm'enu
%                                  - 'e' is the default (missing or empty)
% Labels   - Labels for button and menu types.
%            Default value for string and eval types.
% Values   - [optional] Return values corresponding to Labels
% p        - results
%_______________________________________________________________________
%
% spm_input prompts for input by creating objects in Figure 2 and
% setting the background color to a delicate lilac (bleurgh).  All
% processes stop until something is entered. After the query, the
% (uneditable) response is displayed on a light gray background.
% Previous spm_input transactions at YPos are deleted.
%
% If YPos is 0 or global CMDLINE is true, then the command line is used.
% Negative YPos ensures the GUI is used, at YPos=abs(YPos).
%
% There are four types of query: string, evaluated string, buttons and menu,
% specified by Type strings beginning with the letters 's','e','b' & 'm'
% respectively.
%
% STRING - The user is prompted for string input, using a one line
% editable text window, or on the command line (using input). The string
% is returned when <RETURN> or <ENTER> is pressed.
%
% EVALUATED STRING [DEFAULT] - As for Type STRING, a string is requested.
% This string is evaluated, permitting the entry of numerics, matrices,
% expressions or functions. I.e.:
%     i)    - a number, vector or matrix e.g. "[1 2 3 4]", "[1:4]" or "1:4"
%     ii)   - an expression                     e.g.  "pi^2" 
%     iii)  - a function (that will be invoked) e.g. "data_input1"
% (Evaluation is in the workspace of spm_input.) Errors in string
% evaluation are handled gracefully, and the user prompted to re-enter.
%
% DEFAULT VALUES - For both STRING & EVALUATED input, the Labels
% parameter may be used to specify a default value. In the GUI version,
% this is written into the edit widget on creation, and can be accepted
% by pressing return. In the command line version, the dafault value is
% indicated, and is accepted if an empty line is input.
%
% BUTTONS - (Type(1)=='b') - The user is presented with buttons
% labelled with the strings contained in the rows of the string matrix
% Labels. The label of the depressed button is returned. Keyboard
% shortcuts are avaliable: Each button can be pressed by typing the first
% character of the label (case insensitive). (If these Keys are not
% unique, then the shortcut keys for the buttons are the integer keys
% 1,2,... respectively.) The command line version prompts the user to
% input one of these unique character keys, and returns the corresponding
% label.
%
% POPUP MENU OF OPTIONS - (Type(1)=='m') - The user is presented with a
% popup menu, labelled with the prompt. Pulling down the menu reveals the
% options as specified in Labels. When the OK button is pressed the
% number of the row of the Labels matrix which contains the chosen label
% is returned. (A beep sounds if the user attempts to 'OK' the Prompt!)
% The command line version presents a character menu in the MatLab
% window.
%
% VALUES - A vector (or matrix) of alternative return values
% corresponding to the Labels matrix can be specified. The length (row
% dimension) must match the row dimension of the Labels parameter. In
% this case, the button and menu types return the element (row) of Values
% corresponding to the label chosen.
%
% Labels SHORTCUT - Labels may be specified as a bar delimited string.
% A Labels vector containing the character '|' is expanded into a string
% matrix where the |'s mark the line ends. I.e. 'AnCova|Scaling|None'
% expands to ['AnCova ';'Scaling';'None   '].
%
% Buttons SHORTCUT - If the Type parameter is a bar delimited string, then
% the Type is taken as 'b' with the specified labels, and the next parameter
% (if specified) is taken for the Values.
%
% Yes/No question shortcut - p = spm_input(Prompt,YPos,'y/n') expands to 
% p = spm_input(Prompt,YPos,'b','yes|no','yn'), enabling easy use of
% spm_input for yes/no dialogue. p is returned as 'y' or 'n'.
%
% EXAMPLES - (Specified YPos is overriden if global CMDLINE is true, when
%              the command line versions are used.)
% 	p = spm_input
%		Command line input of an evaluated string, with default prompt.
% 	p = spm_input(str,1)
% 		Evaluated string input, prompted by string str, in position 1
% 		of the dialog figure.
% 	p = spm_input(str,1,'e',0.001)
% 		Evaluated string input, prompted by string str, in position 1
% 		of the dialog figure. Default value of 0.001 offered.
% 	p = spm_input(str,2,'s')
%		String input, position 2 of dialog figure, prompted by str.
% 	p = spm_input(str,3,'y/n')
% 		Yes/No buttons for question in strm in position 3 of dialog.
% 		Returns 'y' or 'n'.
% 	p = spm_input(str,4,'AnCova|Scaling')
% 		Presents two buttons labelled 'AnCova' & 'Scaling', with 
% 		prompt str, in position 4 of the dialog figure. Returns the 
% 		string on the depresed button, where buttons can be pressed 
% 		with the mouse or by the respective keyboard accelerators
%               'a' & 's' (or 'A' & 'S').
% 	p = spm_input(str,4,'b','AnCova|Scaling')
% 		As above.
% 	p = spm_input(str,0,'b','AnCova|Scaling|None',[1,2,3])
%		Prompts for [A]ncova / [S]caling / [N]one in MatLab command
%		window, returns 1, 2, or 3 according to the first character
% 		of the entered string as one of 'a', 's', or 'n' (case 
%               insensitive).
% 	p = spm_input(str,5,'m','Single Subject|Multi Subject|Multi Study')
% 		Prints the prompt str in a pull down menu containing items
% 		'Single Subject', 'Multi Subject' & 'Multi Study'. When OK is
% 		clicked p is returned as the index of the  choice, 1,2, or 3 
% 		respectively.
% 	p = spm_input(str,5,'m',...
% 		'Single Subject|Multi Subject|Multi Study',['SS';'MS';'SP'])
% 		As above, but returns strings 'SS', 'MS', or 'SP' according to
% 		the respective choice.
% 	p = spm_input(str,0,'m',...
% 		'Single Subject|Multi Subject|Multi Study',['SS';'MS';'SP'])
% 		As above, but the menu is presented in the command window as a 
% 		numbered list.
%
%
% See      : input.m (MATLAB Reference Guide)
% See also : spm_get.m and 'Entering variables' in the help facility
%
%__________________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%

% Version History
% - Karl Friston  - V1 -  ??/94
% - Andrew Holmes - V2 -  02/95 - Added button and menu variations, command
%                                 line options, logging & error checking.
% - Andrew Holmes - V2b - 04/95 - Added default value option for string &
%                                 eval types.


%-Condition arguments
%-----------------------------------------------------------------------
if nargin<5 Values=[]; end
if nargin<4, Labels=[]; end
if nargin<3, Type=''; end
if isempty(Type), Type='e'; end
if strcmp(Type,'y/n'), Type = 'b'; Labels='yes|no'; Values='yn'; end
if any(Type=='|'), Values=Labels; Labels=Type; Type='b'; end
if nargin<2, YPos=[]; end
if isempty(YPos), YPos=0; end
if nargin<1, Prompt=''; end
if isempty(Prompt), Prompt='Select '; end

%-Check Type
if ~any(Type(1)=='esbm'), error('Unrecognised question Type'), end

% Background Color
%-----------------------------------------------------------------------
COLOR = [.8,.8,1];


%-Command line? -> If global CMDLINE is true & YPos > 0, or if YPos is empty
%-----------------------------------------------------------------------
CmdLine = 0; global CMDLINE;
if ~isempty(CMDLINE), CmdLine = CMDLINE; end
if YPos==0, CmdLine = 1; end
if YPos<0, CmdLine = 0; YPos=abs(YPos); end

%-Setup
%-----------------------------------------------------------------------
if CmdLine
	%-Print out prompt to command window
	%-------------------------------------------------------------------
	fprintf('\n'), fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
	fprintf('\t%s\n',Prompt)
	fprintf('%c',setstr(ones(1,70)*'=')), fprintf('\n')
else
	%-Determine position of objects
	%-------------------------------------------------------------------
	figure(2)
	set(2,'Units','pixels')
	FigPos = get(2,'Position');
	Xdim = FigPos(3); Ydim = FigPos(4);
	
	a = 5.5/10;
	y = Ydim - 30*YPos;
	PPos = [10,     y, (a*Xdim -20),     20];
	RPos = [a*Xdim, y, ((1-a)*Xdim -10), 20];
	
	%-Delete any previous inputs using position YPos
	%-------------------------------------------------------------------
	delete(findobj(2,'Tag',['GUIinput_',int2str(YPos)]))

	%-Place pointer over control object
	%-------------------------------------------------------------------
	set(0,'PointerLocation',...
		[(FigPos(1) +RPos(1) +RPos(3)/5),  (FigPos(2) +RPos(2) +RPos(4)/2)]);
end % (if CmdLine)



if (Type(1)=='e')|(Type(1)=='s')

DefString = Labels; if ~isstr(DefString), DefString=num2str(DefString); end

%-String and evaluated string input.
%=======================================================================

if CmdLine
	if ~isempty(DefString), Prompt=[Prompt,' (Default: ',DefString,' )']; end
	str = input([Prompt,' : '],'s'); if isempty(str), str=DefString; end
	if Type(1)=='e', p=eval(['[',str,']'],'sprintf(''<ERROR>'')');
		else, p=str; end

	%-Catch eval errors
	%-(Input can do this by itself, but doesn't return the typed string)
	while strcmp(p,'<ERROR>') & Type(1)=='e'
		fprintf('%c! spm_input : evaluation error, please retype\n',7)
		str = input([Prompt,' : '],'s'); if isempty(str), str=DefString; end
		if Type(1)=='e', p=eval(['[',str,']'],'sprintf(''<ERROR>'')');
			else, p=str; end
	end % (while)

else

	%-Create text and edit control objects
	%-------------------------------------------------------------------
	uicontrol(2,'Style','Text',...
		'String',Prompt,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'HorizontalAlignment','Right',...
		'Position',PPos);
	
	%-Construct callback for edit widget - if evaluating input,
	% do it in the callback, so that base workspace variables can be used
	if Type(1)=='e'
		cb = ['set(gco,''UserData'',',...
			'eval([''['',get(gco,''String''),'']''],',...
				' ''sprintf(''''<ERROR>'''')''))'];
	else
		cb = 'set(gco,''UserData'',get(gco,''String'') )';
	end

		%-Edit widget: Callback sets UserData to 1 when edited
	h = uicontrol(2,'Style','Edit',...
		'String',DefString,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'UserData',[],...
		'CallBack',cb,...
		'Horizontalalignment','Left',...
		'BackgroundColor',COLOR,...
		'Position',RPos);

	%-Setup return-keypress CallBack to process contents of edit widget
	%-Necessary for quick acceptance of default values.
	% (Edit widget callback is only executed if an edit is made.)
	set(2,'UserData',h)
	if Type(1)=='e'
		cb = ['set(get(gcf,''UserData''),''UserData'',',...
			'eval([''['',get(get(gcf,''UserData''),''String''),'']''],',...
				' ''sprintf(''''<ERROR>'''')''))'];
	else
		cb = ['set(get(gcf,''UserData''),''UserData'',',...
			'get(get(gcf,''UserData''),''String'') )'];
	end
	set(2,'KeyPressFcn',...
		['if abs(get(gcf,''CurrentCharacter''))==13, ',cb,', end'])

	%-Wait for edit, evaluate string if input not a string variable
	%-------------------------------------------------------------------
	while ~length(get(h,'UserData')), pause(0.5); end

	p = get(h,'UserData');
	if Type(1)=='e'
		while strcmp(p,'<ERROR>')	%-Catch eval errors
			set(h,'String','<ERROR> re-enter...'), drawnow
			fprintf('%c',7), pause(2)
			set(h,'String','','UserData',[])
			while ~length(get(h,'UserData')); pause(0.5); end
			p = get(h,'UserData');
		end % (while)
	end

	%-Fix edit window & clean up
	set(h,'Style','Text','HorizontalAlignment','Center',...
		'BackgroundColor',[.7,.7,.7]), drawnow
	set(2,'UserData',[],'KeyPressFcn','')

end % (if CmdLine)

%-Log the transaction
%---------------------------------------------------------------------------
if exist('spm_log')==2
	if Type(1)=='e'
		spm_log(['spm_input : ',Prompt,': (',str,')'],p);
	else
		spm_log(['spm_input : ',Prompt,':',p]);
	end
end

return

end % (if Type(1))


%-More complicated input - Types 'b'utton and 'm'enu
%=======================================================================

%-Condition arguments
%-----------------------------------------------------------------------
%-Convert Option string to string matrix if required
if any(Labels=='|')
	OptStr=Labels;
	BarPos=find([OptStr=='|',1]);
	Labels=OptStr(1:BarPos(1)-1);
	for Bar = 2:sum(OptStr=='|')+1
		Labels=str2mat(Labels,OptStr(BarPos(Bar-1)+1:BarPos(Bar)-1)); end
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

if isempty(Labels), error('No Labels specified'), end
if (size(Labels,1)~=size(Values,1))
	error('Labels & Values incompatible sizes'), end


if Type(1)=='b'
%=======================================================================
NoLabels = size(Labels,1);

%-Make unique character keys for the Labels, ignoring case
%-----------------------------------------------------------------------
Keys=Labels(:,1)';
if any(~diff(abs(sort(lower(Keys)))))
	Keys=int2str(1:NoLabels);
	Labs=Labels;
else
	Labs=Labels;
	Labs(:,1)=[];
end

if CmdLine

	%-Build prompt
	%-------------------------------------------------------------------
%	%-Original version
%	Prmpt = ['[',Keys(1),']',Labs(1,:),' '];
%	for lab = 2:NoLabels
%		Prmpt=[Prmpt,'/ [',Keys(lab),']',Labs(lab,:),' ']; end
%	Prmpt = [Prmpt,'? '];

	%-TNichols patch
	if ~isempty(Labs) 
	    Prmpt = ['[',Keys(1),']',Labs(1,:),' '];
	    for lab = 2:NoLabels
		    Prmpt=[Prmpt,'/ [',Keys(lab),']',Labs(lab,:),' ']; end
	    Prmpt = [Prmpt,'? '];
	else
	    Prmpt = ['[',Keys(1),'] '];
	    for lab = 2:NoLabels
		    Prmpt=[Prmpt,'/ [',Keys(lab),'] ']; end
	    Prmpt = [Prmpt,'? '];
	end

	%-Ask for user response
	%-------------------------------------------------------------------
	str = input(Prmpt,'s');
	str = [str(str~=' '),' ']; %-Deblank & ensure not empty

	while ~any(lower(Keys)==lower(str(1)))
		fprintf('%c',7)
		str = input(Prmpt,'s'); str = [str(str~=' '),' '];
	end % (while)
	fprintf('\n')

	k = find(lower(Keys)==lower(str(1)));
	p = Values(k,:);

else

	%-Create text and edit control objects
	%-------------------------------------------------------------------
	hPrmpt = uicontrol(2,'Style','Text',...
		'String',Prompt,...
		'Tag',['GUIinput_',int2str(YPos)],...
		'HorizontalAlignment','Right',...
		'Position',PPos);

	dX = RPos(3)/NoLabels;

	H = [];
	for lab=1:NoLabels
		h = uicontrol(2,'Style','Pushbutton',...
			'String',deblank(Labels(lab,:)),...
			'Tag','',...
			'UserData',lab,...
			'BackgroundColor',COLOR,...
			'Callback','set(gco,''Tag'',''Pressed'')',...
			'Position',[RPos(1)+(lab-1)*dX RPos(2) dX 020]);
		H = [H,h];
	end

	%-Enable keyboard shortcuts - use an invisible frame
	hF = uicontrol(2,'Style','Frame',...
		'Visible','off',...
		'Tag','','UserData',[],...
		'Position',PPos);
	H = [H,hF];
	%-Set callback for KeyPress, to store char in UserData of frame
	% Store handle of invisible frame in figures UserData
	cb = ['set(get(gcf,''UserData''),',...
		'''Tag'',''Pressed'',''UserData'',get(gcf,''CurrentCharacter''))'];
	set(2,'KeyPressFcn',cb,'UserData',h)

	%-Wait for button press, process results
	%-------------------------------------------------------------------
	while isempty(findobj(H,'Tag','Pressed')), pause(0.01), end
	h = findobj(H,'Tag','Pressed');
	k = get(h,'UserData');
	if isstr(k), k = find(lower(Keys)==lower(k)); end

	while isempty(k) %-Keyboard input, key pressed wasn't one of the Keys
		set(h,'Tag','')
		while isempty(findobj(H,'Tag','Pressed')), pause(0.01), end
		h = findobj(H,'Tag','Pressed');
		k = get(h,'UserData');
		if isstr(k), k = find(lower(Keys)==lower(k)); end
	end % (while)

	p = Values(k,:); if isstr(p), p=deblank(p); end
	delete(H)
	set(2,'KeyPressFcn','','UserData',[])
	uicontrol(2,'Style','Text','Tag','',...
		'String',deblank(Labels(k,:)),...
		'Tag',['GUIinput_',int2str(YPos)],...
		'Horizontalalignment','Center',...
		'BackgroundColor',[.7,.7,.7],...
		'Position',RPos);

end % (if CmdLine)

elseif Type(1)=='m'
%=======================================================================

if CmdLine
	NoLabels = size(Labels,1);
	for lab = 1:NoLabels
		fprintf('\t%2d : %s\n',lab,Labels(lab,:)), end
	k = input(['Menu choice (1-',int2str(NoLabels),') ? ']);
	while ~any([1:NoLabels]==k)
		fprintf('%c\t!Out of range',7)
		k = input(['Menu choice (1-',int2str(NoLabels),') ? ']);
	end % (while)
	fprintf('\n')

else

	MPos = [PPos(1), PPos(2), Xdim-50, 20];
	BPos = MPos + [Xdim-50+1, 0+1, 50-Xdim+30, 0];
	RPos = MPos + [0,0,30,0];

	hPopUp = uicontrol(2,'Style','PopUp',...
		'HorizontalAlignment','Left',...
		'ForegroundColor','k',...
		'BackgroundColor',COLOR,...
		'String',str2mat(Prompt,Labels),...
		'Position',MPos);
	H = hPopUp;

	cb = ['if (get(get(gco,''UserData''),''Value'')>1),',...
			'set(gco,''Tag'',''Pressed''),',...
			'else, fprintf(''%c'',7), end'];
	hOK = uicontrol(2,'Style','Pushbutton','Tag','',...
		'String','OK',...
		'UserData',hPopUp,...
		'Callback',cb,...
		'Position',BPos);
	H = [H,hOK];

	%-Wait for button press
	%-------------------------------------------------------------------
	while isempty(get(hOK,'Tag')), pause(0.01), end
	k = get(hPopUp,'Value')-1;

	uicontrol(2,'Style','Text',...
		'String',deblank(Labels(k,:)),...
		'Tag',['GUIinput_',int2str(YPos)],...
		'Horizontalalignment','Center',...
		'BackgroundColor',[.7,.7,.7],...
		'Position',RPos);
	delete(H)

end % (if CmdLine)

p = Values(k,:); if isstr(p), p=deblank(p); end

end % (if Type(1)==...)

%-Log the transaction
%-----------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_input : ',Prompt,': (',deblank(Labels(k,:)),')'],p); end
