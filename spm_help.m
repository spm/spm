function R1=spm_help(Action,P2,P3)
% SPM Interactive help and manual facility
% FORMAT spm_help
%_______________________________________________________________________
%
% spm_help sets up a GUI help system for the SPM package. A layout of
% the SPM components and special topics is drawn. The buttons display
% the relevant help in the Graphics window. Referenced routines and
% previous topics are listed in the Interactive window, permitting
% recursive navigation of the package architecture and it's help.
% 
% From the graphics window the help topics can be printed. Multi page
% topics are paged on screen, and printed page by page. (See
% spm_figure)
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%
%
%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take. Recall   )
%( MatLab's command-function duality: `spm_help Menu` is equivalent to )
%( `spm('Menu')`.                                                      )
%
% FORMAT spm_help
% Defaults to spm_help('Menu')
%
% FORMAT spm_help('Menu')
% Finds or creates (spm_help('CreateWin')) the Help window, and then
% makes it visible for use.
%
% FORMAT spm_help('Quit')
% Clears the Graphics and Interactive windows, and hides the Help window.
%
% FORMAT spm_help('CreateWin')
% Creates the help window ('Tag'ed 'Help'), which mirrors the layout of
% the SPM Menu window.
%
% FORMAT spm_('Topic',Fname)
% Fname     - Name of file from which to display help
% Loads file Fname (which must be on the MATLABPATH), parses it for
% references to other spm_* routines, displays these in the Interactive
% window (with callbacks to display their help) (creates one if
% necessary), and sets up a "Previous topics" pulldown in the
% Interactive window along with a window for direct editing of the
% Fname topic to display. Then calls spm_help('Disp') to display the
% help portion of the file in the Graphics window.
%
% FORMAT spm_help('Disp',Fname,S)
% Fname     - Name of file from which to display help
% S         - [Optional] String vector containing a previously read in
%             contents of file Fname
% Displays the help for the given file in the Graphics window (creating
% one if required). Paginates and provides paging buttons if necessary.
%
% FORMAT spm_help('Page')
% Callback handling function for the paging buttons in multi-page help
% displays.
% 
%_______________________________________________________________________


%-Condition arguments for all uses
%-----------------------------------------------------------------------
if nargin==0, Action='Menu'; end

if strcmp(lower(Action),lower('Menu'))
%=======================================================================
% spm_help('Menu')

%-Find Help window, create one if none exists, show it
%-----------------------------------------------------------------------
Fhelp = spm_figure('FindWin','Help');
if isempty(Fhelp), Fhelp=spm_help('CreateWin'); end
set(Fhelp,'Visible','on')
return


elseif strcmp(lower(Action),lower('Quit'))
%=======================================================================
spm_clf('Graphics')
spm_clf('Interactive')
Fhelp = spm_figure('FindWin','Help');
if isempty(Fhelp), return, end
set(Fhelp,'Visible','off')
return


elseif strcmp(lower(Action),lower('CreateWin'))
%=======================================================================
% Fhelp = spm_help('CreateWin');

%-Create help window
%-----------------------------------------------------------------------
Fhelp = spm_figure('FindWin','Help');
if ~isempty(Fhelp), close(Fhelp), end

%-Create window
S     = get(0,'ScreenSize');
A     = [S(3)/1152 S(4)/900 S(3)/1152 S(4)/900];
S4    = [108 429 400 445].*A;
Fhelp = figure('Tag','Help',...
	'Name','SPM Help facility',...
	'Color',spm('Colour'),...
	'Position',S4,...
	'NumberTitle','off',...
	'Resize','off',...
	'Visible','off');
R1    = Fhelp;

%-Create help buttons with callbacks
%----------------------------------------------------------------------------
%-Special overview man pages
uicontrol(Fhelp,'String','About SPM',...
	'Position',[010 410 087 30].*A,...
	'CallBack','spm_help(''Topic'',''spm.man'');',...
	'ForegroundColor','b');
uicontrol(Fhelp,'String','Data format',...
	'Position',[107 410 088 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_format.man'');',...
	'ForegroundColor','b');
uicontrol(Fhelp,'String','Methods',...
	'Position',[205 410 088 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_methods.man'');',...
	'ForegroundColor','b');
uicontrol(Fhelp,'String','Variables',...
	'Position',[303 410 087 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_ui.man'');',...
	'ForegroundColor','b');

uicontrol(Fhelp,'String','PET Overview',...
	'Position',[040 285 140 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_pet.man'');',...
	'ForegroundColor','b');
uicontrol(Fhelp,'String','fMRI Overview',...
	'Position',[220 285 140 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_fmri.man'');',...
	'ForegroundColor','b');

uicontrol(Fhelp,'String','Graphics',...
	'Position',[165 122 070 24].*A,...
	'CallBack','spm_help(''Topic'',''spm_figure.m'');',...
	'ForegroundColor','b');

%-Man pages for specific functions
uicontrol(Fhelp,'String','Realign',...
	'Position',[040 370 080 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_realign.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Normalize',...
	'Position',[150 350 100 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_sn3d.m'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Smooth',...
	'Position',[280 370 080 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_smooth.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Coregister',...
	'Position',[040 330 080 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_coregister.m'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Segment',...
	'Position',[280 330 080 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_segment.m'');',...
	'Interruptible','yes');


uicontrol(Fhelp,'String','PET Statistics',...
	'Position',[040 245 140 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_spm.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','fMRI Statistics',...
	'Position',[040 215 140 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_fmri_spm.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Eigenimages',...
	'Position',[220 245 140 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_svd.man'');',...
	'Interruptible','yes');

uicontrol(Fhelp,'String','SPM{F}',...
	'Position',[045 165 070 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_F.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','SPM{Z}',...
	'Position',[165 165 070 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_projections.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Results',...
	'Position',[285 165 070 30].*A,...
	'CallBack','spm_help(''Topic'',''spm_sections.man'');',...
	'Interruptible','yes');

uicontrol(Fhelp,'String','Analyze',...
	'Position',[020 088 082 024].*A,...
	'CallBack','',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Display',...
	'Position',[112 088 083 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_image.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','<empty>',...
	'Position',[205 088 083 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_help.m'');',...
	'Visible','off',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','<empty>',...
	'Position',[298 088 082 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_help.m'');',...
	'Visible','off',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Mean',...
	'Position',[020 054 082 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_average.m'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','ImCalc',...
	'Position',[112 054 083 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_image_funks.m'');',...
	'Interruptible','yes');
%uicontrol(Fhelp,'String','MRI to PET',...
%	'Position',[205 054 083 024].*A,...
%	'CallBack','spm_help(''Topic'',''spm_mri2pet.m'');',...
%	'Interruptible','yes');
%uicontrol(Fhelp,'String','MRsegment',...
%	'Position',[298 054 082 024].*A,...
%	'CallBack','spm_help(''Topic'',''spm_segment.m'');',...
%	'Interruptible','yes');
uicontrol(Fhelp,'String','Help',...
	'Position',[020 020 082 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_help.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Defaults',...
	'Position',[112 020 083 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_defaults.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String',['<',getenv('USER'),'>'],...
	'Position',[205 020 083 024].*A,...
	'CallBack','spm_help(''Topic'',''spm_button.man'');',...
	'Interruptible','yes');
uicontrol(Fhelp,'String','Quit',...
	'Position',[298 020 082 024].*A,...
	'ForegroundColor','r',...
	'CallBack','spm_help(''Quit'')');
return



elseif strcmp(lower(Action),lower('Topic'))
%=======================================================================
% spm_help('Topic',Fname)
if nargin<2, P2='spm'; else, Fname=P2; end

%-Callback from PrevTopics passes integer Fname
if ~isstr(Fname)
	%-Current object contains new Fname
	Fname = get(gco,'String');
	n     = get(gco,'Value');
	if  n==1, return, end		% Don't process the prompt!
	Fname = deblank(Fname(n,:));
end

%-Find (or create) window to print in
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), Finter = spm('CreateIntWin'); end
set(Finter,'Pointer','Watch')

%-Load text file
%-----------------------------------------------------------------------
fid = fopen(Fname,'r');
if fid<0
	errordlg(['File ',Fname,' not found'],'SPMerror','on');
	set(Finter,'Pointer','Arrow')
	return
end
S   = setstr(fread(fid))';
fclose(fid);

%-Display the current help comments
%-----------------------------------------------------------------------
spm_help('Disp',Fname,S);

if isempty(findobj(Finter,'Tag','SPMhelp'))
	%-Create control objects
	%---------------------------------------------------------------
	spm_figure('Clear',Finter);
	A = spm('GetWinScale');
	set(Finter,'Name','SPM Help')
	uicontrol(Finter,'Style','Frame','Tag','hAxes',...
		'Position',[001 345 400 050].*A);
	uicontrol(Finter,'Style','Text','Tag','SPMhelp',...
		'String','Routines referenced by ',...
		'HorizontalAlignment','Right',...
		'Position',[005 370 165 020].*A)
	uicontrol(Finter,'Style','Edit','Tag','Fname',...
		'String','<filename>',...
		'ForegroundColor','r','BackgroundColor',[.8,.8,1],...
		'HorizontalAlignment','Left',...
		'CallBack','spm_help(''Topic'',get(gco,''String''))',...
		'Position',[170 370 225 020].*A)
	uicontrol(Finter,'Style','PopUp','Tag','PrevTopics',...
		'String',str2mat('Previous Topics...','spm.man'),...
		'HorizontalAlignment','Left',...
		'ForegroundColor','r',...
		'UserData',1,...
		'Callback',...
		'spm_help(''Topic'',get(gco,''Value''))',...
		'Position',[170-1 350-1 225+1 020].*A)
	uicontrol(Finter,'Style','Frame','Tag','StatusArea',...
		'Position',[001 001 400 030].*A);
	uicontrol(Finter,'Style','Text','Tag','StatusLine',...
		'String','Select routine to display help',...
		'HorizontalAlignment','Center',...
		'ForegroundColor','w',...
		'Position',[020 005 360 020].*A)
	
	figure(Finter)
	hAxes = axes('Position',[020 035 320 280]./[400 395 400 395],...
			'Units','Points','Visible','off');
	set(findobj(Finter,'Tag','hAxes'),'UserData',hAxes);
end

%-Sort out control objects
%-----------------------------------------------------------------------
set(findobj(Finter,'Tag','Fname'),'String',Fname);

%-Sort out previous topics pulldown
%-----------------------------------------------------------------------
hPTopics   = findobj(Finter,'Tag','PrevTopics');
PrevTopics = get(hPTopics,'String');
Prompt     = PrevTopics(1,:); PrevTopics(1,:)=[];
PrevTopics = str2mat(Fname,PrevTopics);
%-Delete any replications of Fname within PrevTopics
if size(PrevTopics,1)>1
	IDRows=[0, all( PrevTopics(2:size(PrevTopics,1),:)'==...
		setstr(ones(size(PrevTopics,1)-1,1)*PrevTopics(1,:))' ) ];
	PrevTopics(IDRows,:)=[];
end
%-Truncate to 20 items max.
if size(PrevTopics,1)>20 PrevTopics(21:size(PrevTopics,1),:)=[]; end
%-Update popup
set(hPTopics,'String',str2mat(Prompt,PrevTopics),'Value',1)


%-Sort out axes for printing area
%-----------------------------------------------------------------------
hAxes = get(findobj(Finter,'Tag','hAxes'),'UserData');
delete(get(hAxes,'Children'))
axes(hAxes)
y     = floor(get(hAxes,'Position'));
y0    = y(3);
set(hAxes(1),'Ylim',[0,y0])

%-Find subroutines and create text objects and their CallBacks
%-----------------------------------------------------------------------
q     = findstr(S,'spm_');
y     = y0;
x     = 0;
U     = Fname;
for i = 1:length(q)
    d = [0:32] + q(i);
    Q = S(d(d <= length(S)));
    d = find((Q == ';') | (Q == '(') | (Q == 10) | (Q == '.') | (Q == ' '));
	if length(d)
	  Q   = [Q(1:(min(d) - 1)) '.m'];
	  if exist(Q) == 2 & ~length(findstr(U,Q))
	    text(x,y,Q,...
	    	'ButtonDownFcn',['spm_help(''Topic'',''' Q ''')'])
	    y = y - 18;
	    U = [U ';' Q];
	  end
    end
    if y < 0; x = x + 0.5; y = y0; end
end

set(Finter,'Pointer','Arrow')
return


elseif strcmp(lower(Action),lower('Disp'))
%=======================================================================
% spm_help('Disp',Fname,S)
if nargin<3, S=''; else, S=P3; end
if nargin<2, P2='spm'; else, Fname=P2; end

%-Find (or create) window to print in
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph = spm_figure('Create','Graphics'); end
spm_clf('Graphics')

%-Parse text file/string
%-----------------------------------------------------------------------
if isempty(S)
	fid = fopen(Fname,'r');
	if fid<0
		errordlg(['File ',Fname,' not found'],'SPMerror','on');
		return
	end
	S   = setstr(fread(fid))';
	fclose(fid);
end	
q     = min([length(S),findstr(S,setstr([10 10]))]);	% find empty lines
q     = find(S(1:q(1)) == 10);				% find line breaks

figure(Fgraph)
spm_clf(Fgraph)
hAxes = axes('Position',[0.1,0.05,0.8,0.9],...
		'Units','Points','Visible','off');
y     = floor(get(hAxes(1),'Position'));
y0    = y(3);
set(hAxes(1),'Ylim',[0,y0])
text(-0.1,y0,[Fname,' - help'],'FontSize',16,'FontWeight','bold');
y     = y0 - 24;


%-Loop over pages & lines of text
%-----------------------------------------------------------------------
Vis     = 'on';
FmtLine = 1;
for i = 1:(length(q) - 1)
	d = S((q(i) + 1):(q(i + 1) - 1));
	if d(1) == abs('%');
		%-For some reason, '|' characters cause a CR.
		d = strrep(d,'|','I');
		h = text(0,y,d(2:length(d)),...
			'FontName','Courier','FontSize',10,'Visible',Vis);
		if FmtLine
			set(h,'FontWeight','bold',...
				'FontName','Times','FontSize',12);
			y=y-5;
			FmtLine=0;
		end
		y = y - 7;
	end
	if y<0 %-Start new page
		text(0.5,-10,['Page ',num2str(length(hAxes))],...
			'FontSize',8,'FontAngle','Italic',...
			'Visible',Vis)
		hAxes = [hAxes,axes('Position',[0.1,0.05,0.8,0.9],...
			'Units','Points','Visible','off')];
		set(hAxes(length(hAxes)),'Ylim',[0,y0])
		y     = y0;
		Vis   = 'off';
	end
end
if strcmp(Vis,'off')
	%-Label last page
	text(0.5,-10,['Page ',num2str(length(hAxes))],...
		'FontSize',8,'FontAngle','Italic',...
		'Visible',Vis)
end


%-Create control object if paging is required
%----------------------------------------------------------------------------
if length(hAxes)>1
	uicontrol(gcf,'Style','Pushbutton','String','next',...
		'Callback','spm_help(''Page'',''+1'')',...
		'Position',[520 020 60 020].*spm('GetWinScale'),...
		'Tag','NextPage','UserData',hAxes);
	uicontrol(gcf,'Style','Pushbutton','String','previous',...
		'Callback','spm_help(''Page'',''-1'')',...
		'Position',[020 020 60 020].*spm('GetWinScale'),...
		'Visible','off',...
		'Tag','PrevPage','UserData',1);
end
return

elseif strcmp(lower(Action),lower('Page'))
%=======================================================================
% spm_help('Page',move)
if nargin<2, move=1; else, move=P2; end
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), error('No Graphics window'), end

hNextPage = findobj(Fgraph,'Tag','NextPage');
hPrevPage = findobj(Fgraph,'Tag','PrevPage');
hAxes     = get(hNextPage,'UserData');
Cpage     = get(hPrevPage,'UserData');
nPages    = length(hAxes);

if isstr(move)
	Npage = Cpage + eval(move);
else
	Npage = move;
end

if (Npage<1)|(Npage>nPages), return, end

set(get(hAxes(Cpage),'Children'),'Visible','off')
set(get(hAxes(Npage),'Children'),'Visible','on')
set(hPrevPage,'UserData',Npage)

if Npage==1, set(hPrevPage,'Visible','off')
else, set(hPrevPage,'Visible','on'), end

if Npage==nPages, set(hNextPage,'Visible','off')
else, set(hNextPage,'Visible','on'), end


return


else
%=======================================================================
error('Unknown action string')

%=======================================================================
end
