function varargout=spm_conman(varargin)
% 
% FORMAT [I,xCon] = spm_conman(xX,xCon,STATmode,n,Prompt,Mcstr,OK2chg)
% FORMAT varargout=spm_conman(varargin)
%       - An embedded callback, multi-function function
%       - For detailed programmers comments, see the format specifications
%         below
%_______________________________________________________________________
%
%    ^     
%   / \    
%  / ! \   Help still under construction...
% /_____\  
% 
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%=======================================================================
% - FORMAT specifications for embedded CallBack functions
%=======================================================================
%
%_______________________________________________________________________
% Andrew Holmes


%-Parameters
%=======================================================================
COLOUR   = [.8,.8,1];	%-Question background colour
PJump    = 1;		%-Jumping of pointer to ConMan window


if (nargin==0) | ~ischar(varargin{1})
%=======================================================================
% [I,xCon] = spm_conman(xX,xCon,STATmode,n,Prompt,Mcstr,OK2chg)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<7, OK2chg = 0; else, OK2chg=varargin{7}; end
if nargin<6, Mcstr = ''; else, Mcstr=varargin{6}; end
if nargin<5, Prompt='Select contrast(s)...'; else, Prompt=varargin{5}; end
if nargin<4, n=1; else, n=varargin{4}; end
if nargin<3, STATmode='T|F'; else, STATmode=varargin{3}; end
if nargin<2, xCon=[]; else, xCon=varargin{2}; end
if nargin<1, error('no design specified'), else, xX=varargin{1}; end

%-----------------------------------------------------------------------

%-Setup ConMan window & jump cursor
[F,cF] = spm_conman('Initialise','on',xX,xCon,STATmode,n,Prompt,Mcstr,OK2chg);
if PJump
	PLoc = get(0,'PointerLocation');
	FRec = get(F,'Position');
	set(0,'PointerLocation',[FRec(1)+FRec(3)/2, FRec(2)+FRec(2)/2])
end

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
xCon     = get(F,'UserData');

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

%=======================================================================
% - Callbacks & Utility embedded functions
%=======================================================================

%=======================================================================
switch lower(varargin{1}), case 'initialise'
%=======================================================================
% [F,cF] = spm_conman('Initialise',Vis,xX,xCon,STATmode,n,Prompt,Mcstr,OK2chg)

if nargin<2, Vis='on'; else, Vis=varargin{2}; end

%-Recover ConMan figure number (if it exists)
F  = findobj(get(0,'Children'),'Flat','Tag','ConMan');

cF = get(0,'CurrentFigure');				%-Save current figure

switch lower(Vis), case 'close'
	close(F)
	varargout = {[],cF};
	return
case {'off','reset'}
	varargout = {F,cF};				%-Return figure handles
	if isempty(F), return, end			%-Return if no ConMan win
	set(F,'Visible','off')				%-Make window Invisible
	if strcmp(lower(Vis),'reset')
		set(findobj(F,'Tag','Done'),'UserData',-1)
	end
	return
case 'on'
	%-Sort out arguments
	%---------------------------------------------------------------
	if nargin<9, OK2chg = 0; else, OK2chg=varargin{9}; end
	if nargin<8, Mcstr = ''; else, Mcstr=varargin{8}; end
	Mcstr = cellstr(Mcstr); if length(Mcstr)<2, Mcstr{2}=''; end
	if nargin<7, Prompt='Select contrast(s)'; else, Prompt=varargin{7}; end
	if nargin<6, n=Inf; else, n=varargin{6}; end
	if nargin<5, STATmode='T&F'; else, STATmode=varargin{5}; end
	if nargin<4, error('insufficient arguments'), end
	xCon = varargin{4};
	xX   = varargin{3};

	%-Create/find ConMan window
	%---------------------------------------------------------------
	if isempty(F)					%-Create ConMan win
		[F,H] = spm_conman('CreateFig');
	else						%-Get handles
		H.hDesMtxAx = findobj(F,'Tag','DesMtxAx');
		H.hParEstAx = findobj(F,'Tag','ParEstAx');
		H.hConList  = findobj(F,'Tag','ConList');
		H.hPrompt   = findobj(F,'Tag','Prompt');
		H.hTFAf     = findobj(F,'Tag','TFAf');
		H.hSTATmode = findobj(F,'Tag','STATmode');
		H.hStatLin  = findobj(F,'Tag','StatusLine');
		H.hNew      = findobj(F,'Tag','New');
	end
	varargout = {F,cF};				%-Return figure handles
	
	%-Store required parameters in UserData of various objects
	%---------------------------------------------------------------
	set(F,			'UserData',xCon)
	set(H.hStatLin,		'UserData',Mcstr)	%-**

	%-Initialise interface
	%---------------------------------------------------------------
	set(findobj(F,'Tag','Done'),'UserData',0)	%-Init. Done UserData
	STAT = spm_conman('TFA',F,'',STATmode);		%-Init. TFA buttons
	set(H.hPrompt,'String',Prompt,'UserData',n)	%-Set prompt
	spm_conman('ImDesMtx',xX,H.hDesMtxAx)		%-Depict DesMtx
	spm_conman('ImParEst',xX,H.hParEstAx)		%-Parameter estimability
	spm_conman('ListCon',H.hConList,xCon,STAT,[])	%-List contrasts
	if OK2chg, tmp='on'; else, tmp='off'; end	%-OK to change xCon?
	set(H.hNew,'Enable',tmp)			%-En/dis-able change UI
%-****	if isempty(xCon), spm_conman(); end		%-Go straight to DNewUI

	%-Popup figure, retaining CurrentFigure
	%---------------------------------------------------------------
	set(get(findobj(F,'Tag','DefineNew'),'UserData'),'Visible','off')
							%-Hide define UI
	figure(F)					%-PopUp figure
	set(0,'CurrentFigure',cF)			%-Return to prev. figure
	return

otherwise
	error('Unrecognised ''Vis'' option')
end


%=======================================================================
case 'imdesmtx'
%=======================================================================
% spm_conman('ImDesMtx',xX,h)

h  = varargin{3};

%-Picture design matrix
%-----------------------------------------------------------------------
axes(h)
if isfield(varargin{2},'nKX') & ~isempty(varargin{2}.nKX)
	hDesMtxIm = image((varargin{2}.nKX+1)*32);
else
	hDesMtxIm = image(...
		(spm_DesMtx('sca',varargin{2}.xKXs.X,varargin{2}.Xnames)+1)*32);
end
set(h,'YTick',[],'XTick',[])			%-No Tick marks
set(h,'Tag','DesMtxAx','UserData',varargin{2})	%-Reset axis UserData after image
xlabel('Design matrix')
set(hDesMtxIm,'UserData',...
	struct('X',varargin{2}.xKXs.X,'Xnames',{varargin{2}.Xnames}))
set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')')


%=======================================================================
case 'imparest'
%=======================================================================
% spm_conman('ImParEst',xX,h)

xX = varargin{2};
h  = varargin{3};

%-Picture design parameter estimability
%-----------------------------------------------------------------------
axes(h)
est  = spm_SpUtil('IsCon',xX.xKXs);
nPar = length(est);

hParEstIm = image((est+1)*32);
set(h,	'XLim',[0,nPar]+.5,'XTick',[1:nPar-1]+.5,'XTickLabel','',...
	'YLim',[0,1]+.5,'YDir','reverse','YTick',[],...
	'Box','on','TickDir','in','XGrid','on','GridLineStyle','-');
xlabel('parameter estimability')
set(h,'Tag','ParEstAx')			%-Reset 'Tag' after image cleared it
set(hParEstIm,'UserData',struct('est',est,'Xnames',{xX.Xnames}))
set(hParEstIm,'ButtonDownFcn','spm_DesRep(''SurfEstIm_CB'')')


%=======================================================================
case 'listcon'
%=======================================================================
% spm_conman('ListCon',hConList,xCon,STAT,I)

hConList = varargin{2};
xCon     = varargin{3};
STAT     = varargin{4};
if nargin<5
	Q = get(hConList,'UserData');
	I = Q(get(hConList,'Value'));
else
	I = varargin{5};
end

%-Sort out list, filtering by STAT, and display
%-----------------------------------------------------------------------
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


%=======================================================================
case 'graphcons'
%=======================================================================
% spm_conman('GraphCons',xCon,I,F)

if nargin<2, error('insufficient arguments'), end
xCon = varargin{2};
if nargin<3, I=[1:length(xCon)]; else, I=varargin{3}; end
if nargin<4, F=[]; else, F=varargin{4}; end
if isempty(F), F=spm_figure('FindWin','ConMan'); end
if isempty(F), error('can''t find ConMan win'), end

cF = get(0,'CurrentFigure');		%-Save current figure
set(0,'CurrentFigure',F);		%-Make F current

delete(findobj(F,'Tag','ConGrphAx'));


%-Display contrasts
%-----------------------------------------------------------------------
if isempty(I)
	axes('Position',[0.65 (0.7 + .1*(1-0.9)) 0.3 .1*.9],...
		'Tag','ConGrphAx','Visible','off')
	text(0.5,0.5,'no contrast(s)',...
		'FontSize',spm('FontSize',8),...
		'FontAngle','Italic',...
		'HorizontalAlignment','Center',...
		'VerticalAlignment','Middle')

else

	nPar = size(xCon(1).c,1);
	xx   = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
	nCon = length(I);
	
	dy    = 0.2/max(nCon,2);
	for ii = nCon:-1:1
	    i  = abs(I(ii));
	    axes('Position',[0.65 (0.7 + dy*(nCon-ii+.1)) 0.3 dy*.9])
	    if xCon(i).STAT == 'T' & size(xCon(i).c,2)==1
		%-Single vector contrast for SPM{t} - bar
		yy = [zeros(1,nPar);repmat(xCon(i).c',2,1);zeros(1,nPar)];
		patch(xx,yy,[1,1,1]*.5)
		set(gca,'Box','off','XTick',[],'YTick',[])
		set(gca,'Tag','ConGrphAx')
	    else
		%-F-contrast - image
		sca = 1/max(abs(xCon(i).c(:)));
		image((xCon(i).c'*sca+1)*32)
		set(gca,'Box','on','XTick',[],'YTick',[])
		set(gca,'Tag','ConGrphAx')
	    end
	   if I(ii)>0, ylabel(num2str(i)), end
	   axis tight
	end
	title('Contrast(s)')

end

set(0,'CurrentFigure',cF)		%-Reset CurrentFigure to previous figure


%=======================================================================
case 'statusline'
%=======================================================================
% spm_conman('StatusLine',F,str,col)

if nargin<2,	F = findobj(get(0,'Children'),'Flat','Tag','ConMan');
	else,	F = varargin{2}; end

if nargin<3
	n = get(findobj(F,'Tag','Prompt'),'UserData');
	m = length(get(findobj(F,'Tag','ConList'),'Value'));
	
	str = sprintf('Selected %d contrast',m);
	if m~=1, str=[str,'s']; end
	
	Mcstr = get(findobj(F,'Tag','StatusLine'),'UserData');
	if m>1, str=[str,Mcstr{1}]; else, str=[str,Mcstr{2}]; end
	
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

if nargin<4, col='w'; else, col=varargin{4}; end

set(findobj(F,'Tag','StatusLine'),'String',str,'ForegroundColor',col)


%=======================================================================
case 'done_cb'
%=======================================================================
% spm_conman('Done_CB')

F = gcbf;

n = get(findobj(F,'Tag','Prompt'),'UserData');
q = get(findobj(F,'Tag','ConList'),'Value');


if n>0 & isempty(q)         %-Contrast(s) required, but none selected
	if n==1,            str = 'Select a contrast!';
	elseif isfinite(n), str = sprintf('Select %d contrasts!',n);
	else,               str = 'Select at least one contrast!';
	end
elseif length(q)>abs(n)     %-Too many contrasts selected
	if n<0, tstr='at most'; else, tstr='only'; end
	if abs(n)==1,       str = ['Select',tstr,' one contrast!'];
	else,               str = sprintf('Select %s %d contrasts!',tstr,abs(n));
	end
elseif n>0 & isfinite(n) & length(q)<n
	if n==1,            str = 'Select a contrast!';
	else,               str = sprintf('Select %d contrasts!',n);
	end
else
	str = '';
end

if ~isempty(str)             %-error: display error dialog box
	spm('Beep')
	msgbox(str,sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'error','modal')
else                         %-OK, set Done UserData tag for handling
	set(findobj(F,'Tag','Done'),'UserData',1)
end


%=======================================================================
case 'conlist_cb'
%=======================================================================
% spm_conman('ConList_CB')

F        = gcbf;
hConList = gcbo;

Q        = get(hConList,'UserData');
I        = Q(get(hConList,'Value'));
xCon     = get(F,'UserData');

spm_conman('GraphCons',xCon,I,F)
spm_conman('StatusLine',get(hConList,'Parent'))

if strcmp(get(F,'SelectionType'),'open'), spm_conman('Done_CB'), end


%=======================================================================
case {'tfa','d_tf'}
%=======================================================================
% STAT = spm_conman('TFA')
% STAT = spm_conman('TFA',F,STAT)
% STAT = spm_conman('TFA',F,STAT,STATmode)

D = strcmp(lower(varargin{1}),'d_tf');		%-Handling DefineNew UI?

if nargin<2, F = gcbf; else, F=varargin{2}; end

if nargin<3                   %-Called as CallBack of T or F RadioButton
%-----------------------------------------------------------------------
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

if ~(nargin<4)                           %-Called to set STAT & STATmode
%-----------------------------------------------------------------------
	STATmode = varargin{4};
	b_set = 1;
else
	STATmode = get(findobj(F,'Tag','STATmode'),'UserData');
	b_set = 0;
end


%-Check STATmode & STAT, set temporary STATmode & STAT indicies
%-----------------------------------------------------------------------
STATinfo = struct(...					%-STAT info structure
		'mode',		{ 'T',  'F',  'T|F',    'T&F'},...
		'vSTAT',	{{'T'},{'F'},{'T','F'},{'T','F',''}},...
		'defSTAT',	{ 'T',  'F',  'T',      ''},...
		'Enable',	{	{'on', 'off','off'},...
					{'off','on' ,'off'},...
					{'on', 'on', 'off'},...
					{'on', 'on', 'on' }}	);
i        = find(strcmp(STATmode,{STATinfo.mode}));	%-STATmode index
if isempty(i), error('unknown STATmode'), end		%-Check STATmode valid
if D & i==4, i=3; end					%-Treat 'T&F' as 'T|F'?
if isempty(STAT), STAT=STATinfo(i).defSTAT; end		%-Set STAT as default(?)
j        = find(strcmp(STAT,{'T','F',''}));		%-STAT index
if isempty(j), error('unknown STAT'); end		%-Check known STAT
if ~any(strcmp(STAT,STATinfo(i).vSTAT))			%-Check STAT is
	error('Invalid STAT for this STATmode')		% valid for
end							% this STATmode


%-Set STAT buttons (& Dis/Enable according to STATmode if b_setEnable)
%-----------------------------------------------------------------------
if ~D, Tag='TFA'; else, Tag='D_TF'; end
H = flipud(findobj(F,'Tag',Tag));
set(H(j),			'Value',1)
set(H(setdiff([1:length(H)],j)),'Value',0)
if b_set
	%-Set RadioButton 'Enable' & store STATmode
	set(findobj(F,'Tag','STATmode'),'UserData',STATmode)
	for k=1:length(H), set(H(k),'Enable',STATinfo(i).Enable{k}), end
end


if ~D                    %-Additional UI setup for main contrast manager
%-----------------------------------------------------------------------

    %-Reset ConList for new STAT if callback for TFA button
    %-------------------------------------------------------------------
    if nargin<3
        spm_conman('ListCon',findobj(F,'Tag','ConList'),get(F,'UserData'),STAT)
    end

else     %-Additional UI setup for  DNew contrast definition interface
%-----------------------------------------------------------------------

    %-Get handles of control objects
    %-------------------------------------------------------------------
    hD_ConMtx  = findobj(F,'Tag','D_ConMtx');
    hD_X1cols  = findobj(F,'Tag','D_X1cols');
    hD_ConErrs = findobj(F,'Tag','D_ConErrs');
    hD_ConInfo = findobj(F,'Tag','D_ConInfo');
    HD_Ttxt    = findobj(F,'Tag','D_Ttxt');
    HD_Ftxt    = findobj(F,'Tag','D_Ftxt');
    
    %-Set interface for new STAT
    %-------------------------------------------------------------------
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
%-----------------------------------------------------------------------
varargout = {STAT};



%=======================================================================
case 'fconmenu_cb'
%=======================================================================
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



%=======================================================================
case 'rename_cb'
%=======================================================================
% spm_conman('Rename_CB')

F = gcbf;

hConList = findobj(F,'Tag','ConList');
Q        = get(hConList,'UserData');
i        = get(hConList,'Value');

%-Make sure there's only one selected contrast
%-----------------------------------------------------------------------
if length(i)~=1
	msgbox('Can only rename a single contrast!',...
		sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'error','modal')
	return
end

%-Get contrast structure array and indicies of current contrast
%-----------------------------------------------------------------------
xCon   = get(F,'UserData');
I      = Q(i);

%-Get new name
%-----------------------------------------------------------------------
str = sprintf('Enter new name for contrast %d (currently "%s"):',I,xCon(I).name);
nname  = inputdlg(str,'SPM: Rename contrast',1,{''},'off');
if isempty(nname), return, end

%-Change name in ConList
%-----------------------------------------------------------------------
tmp    = get(hConList,'String');
tmp{i} = strrep(tmp{i},xCon(I).name,nname{1});
set(hConList,'String',tmp)

%-Change name in contrast structure
%-----------------------------------------------------------------------
xCon(I).name = nname{1};
set(F,'UserData',xCon)


%=======================================================================
case 'parsecon'                       %-Parse (cell)string into contrast
%=======================================================================
% [c,I,emsg,imsg,msg] = spm_conman('ParseCon',cstr,X,STAT)
% X is raw DesMtx or space structure

%-Sort out parameters
%-----------------------------------------------------------------------
if nargin<4, STAT='F'; else, STAT=varargin{4}; end
cstr = varargin{2};
X    = varargin{3};
p    = spm_SpUtil('size',X,2);

%-Preliminary parsing of contrast string (cstr) into cellstr
%-----------------------------------------------------------------------
if isempty(cstr)
	varargout = {[],0,{'    <- !empty input'},{},{}};
	return
elseif iscell(cstr)
	if size(cstr,1)~=1, cstr=cstr(:); end
	c    = cstr;
elseif isstr(cstr)
	cstr = cellstr(cstr);
	c    = cstr;
elseif isnumeric(cstr)
	if ndims(cstr)>2, error('matrices only!'), end
	c    = num2cell(cstr,2);
	cstr = cell(size(c));
	for i=1:prod(size(c)), cstr{i}=num2str(c{i}); end
else
	error('contrast input must be string or number')
end


%-Evaluate individual lines of contrast matrix input
%-----------------------------------------------------------------------
I   = zeros(size(c,1),1);
msg = cell(size(c)); [msg{:}] = deal(' (OK)');
for i=1:size(c,1)
    if isstr(c{i})
        c{i} = evalin('base',['[',c{i},']'],'''!''');
    end
    if isstr(c{i})
        msg{i} = '!evaluation error';
    else
        if isempty(c{i})
            msg{i}=' (empty line  - ignored)';
        elseif all(c{i}(:)==0)
            if size(c{i},1)==1, str='vector'; else, str='matrix'; end
            c{i}=[]; msg{i}=[' (zero ',str,' - ignored)'];
        elseif STAT=='T' & size(c{i},1)>1
            c{i}='!'; msg{i}='!vector required';
        elseif size(c{i},2)>p
            c{i}='!'; msg{i}=sprintf('!too long - only %d prams',p);
        else
            if size(c{i},2)<p
                tmp = p-size(c{i},2);
                c{i}=[c{i}, zeros(size(c{i},1),tmp)];
                if size(c{i},1)>1, str=' column'; else, str=''; end
                if tmp>1,          str=[str,'s']; end
                msg{i} = sprintf(' (right padded with %d zero%s)',tmp,str);
            end
            if ~spm_SpUtil('allCon',X,c{i}')
                c{i}='!'; msg{i}='!invalid contrast';
            end
        end
    end
    I(i)=~isstr(c{i});
end

%-Construct contrast matrix, validity indicator, and collate parsing messages
%-----------------------------------------------------------------------
c    = cat(1,c{find(I)});
msg  = [char(cstr), repmat('  <-  ',size(msg,1),1), char(msg)];
emsg = msg(find(~I),:);
imsg = msg(find( I),:);

if all(I) & STAT=='T' & size(c,1)>1      %-Check for vector t-contrasts!
	I=zeros(size(I)); emsg={'!vector required'}; imsg={};
end

%-Return arguments
%-----------------------------------------------------------------------
varargout = {c',I,emsg,imsg,msg};



%=======================================================================
case 'parseistr'                     %-Parse index string
%=======================================================================
% [i,msg] = spm_conman('ParseIStr',str,max)
% X is raw DesMtx or space structure
% str should be a string (row)vector


%-Sort out parameters
%-----------------------------------------------------------------------
str = varargin{2};
max = varargin{3};


%-Process input string
%-----------------------------------------------------------------------
I = evalin('base',['[',str,']'],'''!''');

if isstr(I)
	varargout = {[],0,[str,'  <- !evaluation error'],''};
	return
end

%-Construct list of valid indicies
%-----------------------------------------------------------------------
ok = ismember(I(:)',[1:max]);
i  = I(ok);

%-Construct diagnostic info messages (if required)
%-----------------------------------------------------------------------
if nargout>1
	str = ''; msg='';
	if any(ok)
		str = strvcat(str,num2str(I(ok)));
		msg = strvcat(msg,'  <-  (OK)');
	end
	tmp = ( I<1 | I>max );			%-Out of range
	if any(tmp)
		str = strvcat(str,num2str(I(tmp)));
		msg = strvcat(msg,sprintf('  <-  (ignored - not in [1:%d]',max));
	end
	tmp = ( ~tmp & ~ok );			%-Non integer in range
	if any(tmp)
		str = strvcat(str,num2str(I(tmp)));
		msg = strvcat(msg,'  <-  (ignored - non-integer)');
	end
	msg = cellstr([str,msg]);
else
	msg={};
end

%-Return arguments
%-----------------------------------------------------------------------
varargout = {i,msg};


%=======================================================================
case 'reset_cb'
%=======================================================================
% spm_conman('Reset_CB')

hConList = findobj(gcbf,'Tag','ConList');
xCon     = get(gcf,'UserData');
STAT     = get(findobj(gcbf,'Tag','TFA','Value',1),'UserData');

spm_conman('ListCon',hConList,xCon,STAT,[])


%=======================================================================
case 'd_setup_cb'
%=======================================================================
% spm_conman('D_Setup_CB')

F        = gcbf;
STAT     = get(findobj(F,'Tag','TFA','Value',1),'UserData');
STATmode = get(findobj(F,'Tag','STATmode'),'UserData');

set(F,'UIContextMenu',[])				%-Disable Fig ContextMenu
H = get(findobj(F,'Tag','DefineNew'),'UserData');	%-Get define UI handles
set(findobj(H,'flat','Tag','D_name'),'String','')	%-Clear name
%set(findobj(H,'flat','Tag','D_ConMtx'),'UserData',[])	%-Clear con definition
set(H,'Visible','on')					%-Show UI
spm_conman('D_TF',F,STAT,STATmode);			%-Setup rest of define UI


%=======================================================================
case {'d_conmtx_cb','d_x1cols_cb'}
%=======================================================================
% spm_conman('D_ConMtx_CB')
% spm_conman('D_X1cols_CB')

fcn = find(strcmp(lower(varargin{1}),{'d_conmtx_cb','d_x1cols_cb'}));

F   = gcbf;
h   = gcbo;
str = get(h,'String');

hD_ConMtx  = findobj(F,'Tag','D_ConMtx');
hD_X1cols  = findobj(F,'Tag','D_X1cols');


%-Extract info from holding objects
%-----------------------------------------------------------------------
xX   = get(findobj(F,'Tag','DesMtxAx'),'UserData');
STAT = get(findobj(F,'Tag','D_TF','Value',1),'UserData');


if fcn==1                              %-Parse string from ConMtx widget
%-----------------------------------------------------------------------

	set(hD_X1cols,'String','')
	[c,I,emsg,imsg] = spm_conman('ParseCon',str,xX.xKXs,STAT);
	if all(I)
		DxCon = spm_FcUtil('Set','',STAT,'c',c,xX.xKXs);
	else
		DxCon = [];
	end

elseif fcn==2               %-Process column indicies from X1cols widget
%-----------------------------------------------------------------------
	set(hD_ConMtx,'String','')

	nPar     = spm_SpUtil('size',xX.xKXs,2);
	[i,imsg] = spm_conman('ParseIStr',str,nPar);

	try	%-try-catch block for any errors in spm_FcUtil!
		DxCon = spm_FcUtil('Set','',STAT,'iX0',i,xX.xKXs);
		if STAT=='T' & size(DxCon.c,2)>1
			I = 0; emsg = {'! t-contrasts must be vectors'};
		else
			I = 1; emsg = '';
		end			
	catch
		I    = 0;
		emsg = lasterr;
	end
end


%-Define the contrast or report errors...
%-----------------------------------------------------------------------
set(findobj(F,'Tag','D_ConErrs'),'String',emsg,'Value',[])
set(findobj(F,'Tag','D_ConInfo'),'String',imsg,'Value',[])
if all(I)
	set(hD_ConMtx,'UserData',DxCon);		%-Store contrast
	spm_conman('GraphCons',DxCon,-1,F)		%-Depict contrast
else
	set(hD_ConMtx,'UserData',[]);			%-Clear contrast store
	spm_conman('GraphCons',[],[],F)			%-Clear contrast plot
end
spm_conman('D_Status',F)				%-Set StatusLine


%=======================================================================
case 'd_reset_cb'
%=======================================================================
% spm_conman('D_Reset_CB')

STAT = get(findobj(gcbf,'Tag','TFA','Value',1),'UserData');
set(findobj(gcbf,'Tag','D_name'),'String','')		%-Clear name
set(findobj(gcbf,'Tag','D_ConMtx'),'UserData',[])	%-Contrast definition
spm_conman('D_TF',gcbf,STAT);				%-Setup rest of define UI


%=======================================================================
case 'd_cancel_cb'
%=======================================================================
% spm_conman('D_Cancel_CB')

set(get(findobj(gcbf,'Tag','DefineNew'),'UserData'),'Visible','off')
set(gcbf,'UIContextMenu',findobj(gcbf,'Tag','ConMan_ConMen'))
spm_conman('StatusLine')


%=======================================================================
case 'd_ok_cb'
%=======================================================================
% spm_conman('D_OK_CB')

F = gcbf;

name  = get(findobj(F,'Tag','D_name'),'String');
DxCon = get(findobj(F,'Tag','D_ConMtx'),'UserData');
STAT  = get(findobj(F,'Tag','D_TF','Value',1),'UserData');

dNam = ~isempty(name);
dCon = ~isempty(DxCon);

if ~(dNam & dCon)
	spm('Beep')
	str = {	'contrast name not defined!','',...
		'no valid contrast defined!',''};
	msgbox(str([dNam+1,dCon+3]),...
		sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'error','modal')
	return
end


%-Append new contrast to xCon structure of ConMan figure 'UserData'
%-----------------------------------------------------------------------
DxCon.name = name;
if ~strcmp(DxCon.STAT,STAT), error('STAT & DxCon.STAT mismatch!'), end
xCon     = get(F,'UserData');
if isempty(xCon)
	xCon = DxCon;
else
	xCon = [xCon, DxCon];
end
set(F,'UserData',xCon);


%-Redisplay the new list of contrasts, with the new one selected
%-----------------------------------------------------------------------
hConList = findobj(F,'Tag','ConList');
I = length(xCon);

%-Use this code to add the new contrast to a multiple selection, if allowed
% Q        = get(hConList,'UserData');
% I        = Q(get(hConList,'Value'));
% n        = get(findobj(F,'Tag','Prompt'),'UserData');
% if abs(n)>1, I=[I,length(xCon)]; else, I=length(xCon); end

spm_conman('TFA',F,xCon(end).STAT);			%-Set STAT
spm_conman('ListCon',hConList,xCon,xCon(end).STAT,I)	%-ListCon

%-Hide the DefineNew UI
%-----------------------------------------------------------------------
set(get(findobj(gcbf,'Tag','DefineNew'),'UserData'),'Visible','off')
set(gcbf,'UIContextMenu',findobj(gcbf,'Tag','ConMan_ConMen'))


%=======================================================================
case 'd_status'
%=======================================================================
% spm_conman('D_Status',F)

if nargin<2, F=gcbf; else, F=varargin{2}; end
str  = {' not',''};
dNam = ~isempty(get(findobj(F,'Tag','D_name'),'String'));
dCon = ~isempty(get(findobj(F,'Tag','D_ConMtx'),'UserData'));
if dNam & dCon, ok='on'; col='g'; else, ok='off'; col='w'; end
spm_conman('StatusLine',F,...
	sprintf('name%s defined, contrast%s defined',str{dNam+1},str{dCon+1}),...
	col)
%set(findobj(F,'Tag','D_OK'),'Enable',ok)		%-Enable "OK" button?


%=======================================================================
case 'createfig'
%=======================================================================
% [F,H] = spm_conman('CreateFig')
% Handle Structure - H.{hConList,hDesMtxAx,hPrompt,hSTATmode,hStatLin,hNew}

cF = get(0,'CurrentFigure');		%-Save current figure

%-Create window, compute scaling for screen size
%-----------------------------------------------------------------------
WS = spm('WinScale');				%-Window scaling factors
FS = spm('FontSizes');				%-Scaled font sizes
PF = spm_platform('fonts');			%-Font names (for this platform)
S0 = get(0,'ScreenSize');			%-Screen size

F = figure('IntegerHandle','off',...
		'Tag','ConMan',...
		'Name','SPM contrast manager','NumberTitle','off',...
		'Position',[S0(3)/2,S0(4)/2,0,0] + [-250,-200,500,400].*WS,...
		'Resize','off',...
		'Color',[1 1 1]*.7,...
		'MenuBar','none',...
		'DefaultTextColor','k',...
		'DefaultTextFontName',PF.helvetica,...
		'DefaultTextFontSize',FS(10),...
		'DefaultUicontrolBackgroundColor',[1 1 1]*.7,...
		'DefaultUicontrolFontName',PF.helvetica,...
		'DefaultUicontrolFontSize',FS(10),...
		'DefaultUicontrolInterruptible','on',...
		'Colormap',gray(64),...
		'Renderer','zbuffer',...
		'Visible','off');

%-Draw GUI objects
%-----------------------------------------------------------------------
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

hConList = uicontrol(F,'Style','ListBox','Tag','ConList',...
		'ToolTipString',['Select contrast(s) - drag/shift-click',...
			'/ctrl-click to select multiple contrasts'],...
		'String',{'list','not','set','yet'},...
		'Max',2,...
		'CallBack','spm_conman(''ConList_CB'')',...
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

if exist('spm_help.m')==2
	uicontrol(F,'Style','Pushbutton','String','?',...
		'ToolTipString','help on contrasts and the contrast manager',...
		'ForegroundColor','g',...
		'Callback','spm_help(''spm_conman.m'')',...
		'Position',[460 015 020 020].*WS);
end

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
uimenu(h,'Label','help','Separator','on',...
	'CallBack','spm_help(''spm_conman'')',...
	'Interruptible','off','BusyAction','Cancel');
uimenu(h,'Label','crash out','Separator','on',...
	'CallBack','spm_conman(''Initialise'',''reset'');',...
	'Interruptible','off','BusyAction','Cancel');
set(h,'CallBack','spm_conman(''FConMenu_CB'')',...
	'Interruptible','off','BusyAction','Cancel');
	

%-Draw contrast definition GUI
%-----------------------------------------------------------------------
H = [];				%-Save handles for visibility switching

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
		'FontSize',FS(8),...
		'ForegroundColor','c',...
		'Position',[245 180 050 020].*WS);
H = [H,h];

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
%-----------------------------------------------------------------------
set(0,'CurrentFigure',cF)
varargout = {F,struct(	'hConList',	hConList,...
			'hDesMtxAx',	hDesMtxAx,...
			'hParEstAx',	hParEstAx,...
			'hPrompt',	hPrompt,...
			'hSTATmode',	hSTATmode,...
			'hStatLin',	hStatLin,...
			'hNew',		hNew	)};


%=======================================================================
otherwise                                               %-unknown action
%=======================================================================
error(['Illegal Action string: ',varargin{1}])


%=======================================================================
end                                                            % - E N D
%=======================================================================
