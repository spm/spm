function varargout = spm_api_erp(varargin)
% SPM_API_ERP Application M-file for spm_api_erp.fig
%    FIG = SPM_API_ERP launch spm_api_erp GUI.
%    SPM_API_ERP('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 03-Sep-2004 14:08:52

if nargin == 0  % LAUNCH GUI

	fig     = openfig(mfilename,'reuse');
    Fgraph  = spm_figure('GetWin','Graphics');
	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    
    DCM.name = 'ERP';
    DCM.Y.xy = {};
    DCM.Y.y  = [];
    
    handles.DCM    = DCM;
    handles.Fgraph = Fgraph;
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.


% --------------------------------------------------------------------
function varargout = swd_Callback(h, eventdata, handles, varargin)
try
    cd(get(handles.swd,'String'));
    set(handles.swd,'String',pwd)
end
handles.DCM.swd = pwd;


% --------------------------------------------------------------------
function varargout = cd_Callback(h, eventdata, handles, varargin)

try
    cd(spm_get(-1,[],'select SWD'));
    set(handles.swd,'String',pwd)
end
handles.DCM.swd = pwd;


% --------------------------------------------------------------------
function varargout = name_Callback(h, eventdata, handles, varargin)
handles.DCM.name   = get(handles.name,'String');
set(handles.save,'enable','on')
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = Y_Callback(h, eventdata, handles, varargin)

% get y
%---------------------------------------------------------------------
[f,p] = uigetfile({'*.mat';'*.txt';'*.dat'},'please select ERP data matrix');
f     = fullfile(p,f);
try
    y = load(f,'-ascii');
end
try
    y = load(f,'-mat');
    while ~isnumeric(y)
        s = fieldnames(y);
        y = getfield(y,s{1});
    end
end

% check the number of channels is consistent with the lead field
%---------------------------------------------------------------------
try
    if size(y,2) ~= size(handles.DCM.L,1)
        warndlg('lead field and channel numbers are inconistent')
        return
    end
end

% check the number of channels is consistent previous ERP
%---------------------------------------------------------------------
try
    if size(y,2) ~= size(handles.DCM.Y.xy{end},2)
        warndlg('channels are inconistent')
        return
    end
end

% save in DCM.Y and update X
%---------------------------------------------------------------------
try

    handles.DCM.Y.xy{end + 1} = y;
    m                         = length(handles.DCM.Y.xy);
    spm_dcm_erp_results(handles.DCM,'Response');
    
    % default design
    %-----------------------------------------------------------------
    X               = eye(m,m);
    X(:,1)          = [];
    handles.DCM.U.X = X;
    guidata(h,handles);
    set(handles.X,'String',num2str(X))
    spm_api_erp('X_Callback',gcbo,[],handles)
    
catch
    warndlg('these data are not compatible')
end


% --------------------------------------------------------------------
function varargout = h_Callback(h, eventdata, handles, varargin)
handles.DCM.Y.h    = get(handles.h,'Value');
guidata(h,handles);


% --------------------------------------------------------------------
function varargout = dt_Callback(h, eventdata, handles, varargin)
try
    handles.DCM.Y.dt = str2num(get(handles.dt,'String'));
    guidata(h,handles);
    set(handles.dt,'String',num2str(handles.DCM.Y.dt))
catch
    set(handles.dt,'String',{})
end

% --------------------------------------------------------------------
function varargout = X_Callback(h, eventdata, handles, varargin)
try
    X     = get(handles.X,'String');
    X     = str2num(X);
    [n m] = size(X);
    if  m
        if n == length(handles.DCM.Y.xy) | ~n
            handles.DCM.U.X    = X;
            guidata(h,handles);
            Uname_Callback(h, eventdata, handles);
            reset_Callback(h, eventdata, handles);
        else
            warndlg('please ensure the number of rows and ERPs correspond')
            set(handles.X,'String',num2str(handles.DCM.U.X))
        end
    else
        set(handles.X,'String',num2str(handles.DCM.U.X))
    end
catch
    set(handles.X,'String',num2str(handles.DCM.U.X))
end


% --------------------------------------------------------------------
function varargout = Uname_Callback(h, eventdata, handles, varargin);
try
    name = get(handles.Uname,'String');
    name = name(1:size(handles.DCM.U.X,2));

catch
    name  = {};
    for i = 1:1:size(handles.DCM.U.X,2)
        name{i} = sprintf('input %i',i);
    end
end
set(handles.Uname,'String',name);
handles.DCM.U.name = name;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = Snames_Callback(h, eventdata, handles, varargin)
try
    name = get(handles.Sname,'String');
    name = name(1:size(handles.DCM.L,2));
catch
    name  = {};
    for i = 1:size(handles.DCM.L,2)
        name{i} = sprintf('source %i',i);
    end

end
set(handles.Sname,'String',name);
handles.DCM.Sname = name;
guidata(h,handles);

% --------------------------------------------------------------------
function varargout = L_Callback(h, eventdata, handles, varargin)
[f,p] = uigetfile({'*.mat';'*.txt';'*.dat'},'please select lead feild matrix');
f     = fullfile(p,f);
try
    L = load(f,'-ascii');
end
try
    L = load(f,'-mat');
    while ~isnumeric(L)
        s = fieldnames(L);
        L = getfield(L,s{1});
    end
end

% check the number of channels is consistent with the lead field
%---------------------------------------------------------------------
try
    if size(L,1) ~= size(handles.DCM.Y.xy{1},2)
        warndlg('lead field and channel numbers are inconistent')
        return
    end
end

% and save in DCM.L
%---------------------------------------------------------------------
try
    handles.DCM.L  = L;
    guidata(h,handles);
    Snames_Callback(h, eventdata, handles);
    reset_Callback(h, eventdata, handles);
catch
    warndlg('these data are not compatible')
end

% --------------------------------------------------------------------
function varargout = connections_Callback(h, eventdata, handles, varargin)
n = size(handles.DCM.L,2);    % number of sources
m = size(handles.DCM.U.X,2);  % number of inputs


% reset connection buttons
%---------------------------------------------------------------------
reset_Callback(h, eventdata, handles, varargin)

% connection buttons
%---------------------------------------------------------------------
set(handles.connections,'Units','Normalized')
p  = get(handles.connections,'Position');
x0 = p(1);
y0 = p(2);
sx = 1/32;
sy = 1/64;
for i = 1:n
    for j = 1:n
        for k = 1:3
            x          = x0 + (j - 1)*sx + (k - 1)*(n + 1)*sx;
            y          = y0 - (i + 4)*sy;
            str        = sprintf('data.DCM.A{%i}(%i,%i)',k,i,j);
            str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
            A{k}(i,j)  = uicontrol(handles.SPM,...
                      'Units','Normalized',...
                      'Position',[x y sx sy],...
                      'Style','radiobutton',...
                      'Tag','',...
                      'Callback',str);
            if i == j
                set(A{k}(i,j),'Enable','off')
            end
            try
                set(A{k}(i,j),'Value',handles.DCM.A{k}(i,j));
            catch
                handles.DCM.A{k}(i,j) = get(A{k}(i,j),'Value');
            end
        end
        for k = 1:m
            x          = x0 + (j - 1)*sx + (k - 1)*(n + 1)*sx;
            y          = y0 - (i + 4)*sy - (n + 1)*sy;
            str        = sprintf('data.DCM.B{%i}(%i,%i)',k,i,j);
            str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
            B{k}(i,j)  = uicontrol(handles.SPM,...
                      'Units','Normalized',...
                      'Position',[x y sx sy],...
                      'Style','radiobutton',...
                      'Tag','',...
                      'Callback',str);
            if i == j
                set(B{k}(i,j),'Enable','off')
            end
            try
                set(B{k}(i,j),'Value',handles.DCM.B{k}(i,j));
            catch
                handles.DCM.B{k}(i,j) = get(B{k}(i,j),'Value');
            end
        end
    end
end
for i = 1:n
    x          = x0 + (4 - 1)*(n + 1)*sx;
    y          = y0 - (i + 4)*sy;
    str        = sprintf('data.DCM.C(%i)',i);
    str        = ['data=guidata(gcbo);' str '=get(gcbo,''Value'');guidata(gcbo,data)'];
    C(i)       = uicontrol(handles.SPM,...
                      'Units','Normalized',...
                      'Position',[x y sx sy],...
                      'Style','radiobutton',...
                      'Tag','',...
                      'Callback',str);
    try
        set(C(i),'Value',handles.DCM.C(i));
    catch
        handles.DCM.C(i,1) = get(C(i),'Value');
    end
end

% string labels
%-------------------------------------------------------------------------
constr         = {'A forward' 'A backward' 'A lateral' 'C input'};
nsx   = (n + 1)*sx;
nsy   = 2*sy;
for k = 1:4
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 3*sy;
    str        = constr{k};
    S(k)       = uicontrol(handles.SPM,...
                'Units','Normalized',...
                'Position',[x y nsx nsy],...
                'HorizontalAlignment','left',...
                'Style','text',...
                'String',str);
end
constr         = get(handles.Uname,'String');
for k = 1:m
    x          = x0 + (k - 1)*nsx;
    y          = y0 - 6*sy - 2*(n + 1)*sy;
    str        = ['B ' constr{k}];
    S(4 + k)   = uicontrol(handles.SPM,...
                'Units','Normalized',...
                'Position',[x y nsx nsy],...
                'HorizontalAlignment','left',...
                'Style','text',...
                'String',str);
end
handles.S = S;
handles.A = A;
handles.B = B;
handles.C = C;
set(handles.connections,'Enable','off')
set(handles.reset,      'Enable','on' )
guidata(h,handles)


% --------------------------------------------------------------------
function varargout = save_Callback(h, eventdata, handles, varargin)

DCM         = handles.DCM;
[f,p]       = uiputfile(['DCM' DCM.name '.mat'],'save DCM');
save(fullfile(p,f),'DCM')
set(handles.estimate,'Enable','on')
set(handles.results, 'Enable','off')
set(handles.contrast,'Enable','off')
try
    handles.DCM.F;
    set(handles.estimate,'String','re-estimate')
end


% --------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)

[f,p]       = uigetfile('DCM*.mat','please select DCM');
DCM         = load(fullfile(p,f),'-mat');
DCM         = DCM.DCM;
handles.DCM = DCM;
guidata(h,handles)

try,set(handles.name, 'String',DCM.name),                end
try,set(handles.Y,    'Enable','off'),                   end
try,spm_dcm_erp_results(DCM,'Data'),                     end
try,set(handles.h,    'Value', DCM.Y.h),                 end
try,set(handles.dt,   'String',num2str(DCM.Y.dt)),       end
try,set(handles.X,    'String',num2str(DCM.U.X)),        end
try,set(handles.Uname,'String',DCM.U.name),              end
try,set(handles.Sname,'String',DCM.Sname),               end
try,spm_api_erp('connections_Callback',gcbo,[],handles), end

try,handles.DCM.F; set(handles.contrast,'Enable','on'),
                   set(handles.results, 'Enable','on'),
                   set(handles.estimate,'String',...
                   sprintf('F = %.3f',handles.DCM.F)),   end


% --------------------------------------------------------------------
function varargout = estimate_Callback(h, eventdata, handles, varargin)
set(handles.estimate,'String','Estimating','Foregroundcolor',[1 0 0])
handles.DCM = spm_dcm_erp(handles.DCM);

save_Callback(h, eventdata, handles, varargin)
set(handles.estimate,   'String',sprintf('F = %.3f',handles.DCM.F),...
                        'Callback',{})

set(handles.contrast,   'Enable','on' )
set(handles.results,    'Enable','on' )
set(handles.save,       'Enable','off')
set(handles.load,       'Enable','off')
set(handles.Y,          'Enable','off')
set(handles.L,          'Enable','off')
set(handles.name,       'Enable','off')
set(handles.swd,        'Enable','off')
set(handles.reset,      'Enable','off')

set(handles.connections,'String','specify contrast',...
                        'enable','on',...
                        'Foregroundcolor',[1 0 0],...
                        'Callback',{})

% reset and disable connections for contrast specification
%-------------------------------------------------------------------
for i = 1:length(handles.A)
    for j = 1:length(handles.A{i})
        for k = 1:length(handles.A{i})
            if ~handles.DCM.A{i}(j,k)
                set(handles.A{i}(j,k),'enable','off');
            end
            set(handles.A{i}(j,k),'Value',0);
            handles.DCM.A{i}(j,k)       = 0;
        end
    end
end
for i = 1:length(handles.B)
    for j = 1:length(handles.B{i})
        for k = 1:length(handles.B{i})
            if ~handles.DCM.B{i}(j,k)
                set(handles.B{i}(j,k),'enable','off');
            end
            set(handles.B{i}(j,k),'Value',0);
            handles.DCM.B{i}(j,k)       = 0;
        end
    end
end
for i = 1:length(handles.C)
    set(handles.C(i),'Value',0);
    set(handles.C(i),'enable','off');
    handles.DCM.C(i)       = 0;
end
guidata(h,handles)

% --------------------------------------------------------------------
function varargout = contrast_Callback(h, eventdata, handles, varargin)

% get contrast
%---------------------------------------------------------------------
[n m] = size(handles.DCM.C);
c     = handles.DCM.Ep*0;
c     = spm_erp_pack(c,m,n);
c.A   = handles.DCM.A;
c.B   = handles.DCM.B;
c.C   = handles.DCM.C;
c     = spm_vec(c);
c     = c/sum(c);

warning off
Ec    = full(c'*handles.DCM.Ep);
Cc    = full(c'*handles.DCM.Cp*c);
x     = Ec + [-32:32]*sqrt(Cc)*4/32;
pdf   =     spm_Npdf(x,Ec,Cc);
PP    = 1 - spm_Ncdf(0,abs(Ec),Cc);
warning on

% plot conditional density
%---------------------------------------------------------------------
figure(handles.Fgraph)
clf
subplot(2,1,1);
plot(x,pdf);
title({'Posterior density of contrast',...
        sprintf('P(contrast > 0) = %.1f%s',PP*100,'%')},...
       'FontSize',12)
xlabel('contrast')
ylabel('probability density')
axis square, grid on

% --------------------------------------------------------------------
function varargout = results_Callback(h, eventdata, handles, varargin)
Action  = get(handles.results,'String');
Action  = Action{get(handles.results,'Value')};
spm_dcm_erp_results(handles.DCM,Action);

% remove existing buttons
%---------------------------------------------------------------------
function varargout = reset_Callback(h, eventdata, handles, varargin)
try
    for i = 1:length(handles.A)
        for j = 1:length(handles.A{i})
            for k = 1:length(handles.A{i})
                delete(handles.A{i}(j,k));
            end
        end
    end
    for i = 1:length(handles.B)
        for j = 1:length(handles.B{i})
            for k = 1:length(handles.B{i})
                delete(handles.B{i}(j,k));
            end
        end
    end
    for i = 1:length(handles.C)
        delete(handles.C(i));
    end
    for i = 1:length(handles.S)
        delete(handles.S(i));
    end
    handles.DCM = rmfield(handles.DCM,{'A','B','C'});
    handles     = rmfield(handles,{'A','B','C','S'});
    guidata(h,handles)

end
set(handles.connections,'Enable','on')