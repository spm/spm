function [t,sts] = spm_select(varargin)
% File selector
% FORMAT [t,sts] = spm_select(n,typ,mesg,sel,wd)
%     n    - Number of files
%            A single value or a range.  e.g.
%            1       - Select one file
%            Inf     - Select any number of files
%            [1 Inf] - Select 1 to Inf files
%            [0 1]   - select 0 or 1 files
%            [10 12] - select from 10 to 12 files
%     typ  - file type
%           'any'   - all files
%           'image' - Image files (".img" and ".nii")
%                     Note that it gives the option to select
%                     individual volumes of the images.
%           'xml'   - XML files
%           'mat'   - Matlab .mat files
%           'batch' - SPM batch files (.mat and XML)
%           'dir'   - select a directory
%           Other strings act as a filter to regexp.  This means
%           that e.g. DCM*.mat files should have a typ of '^DCM.*\.mat$'
%      mesg - a prompt (default 'Select files...')
%      sel  - list of already selected files
%      wd   - Directory to start off in
%
%      t    - selected files
%      sts  - status (1 means OK, 0 means window quit)
%
% Files can be selected from disk, but "virtual" files can also be selected.
% Virtual filenames are passed by
%     spm_select('addvfiles',list)
%         where list is a cell array of filenames
% The list can be cleared by
%     spm_select('clearvfiles')
%
% FORMAT cpath = spm_select('CPath',path,cwd)
% function to canonicalise paths: Prepends cwd to relative paths, processes
% '..' & '.' directories embedded in path.
% path     - string matrix (or cell array of strings) containing path names
% cwd      - current working directory [defaut '.']
% cpath    - conditioned paths, in same format as input path argument
%____________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$


if nargin>1 && ischar(varargin{1}) && strcmpi(varargin{1},'addvfiles'),
    vfiles('add',varargin{2:end});
    return;
elseif nargin>0 && ischar(varargin{1}) && strcmpi(varargin{1},'clearvfiles'),
    vfiles('clear');
elseif nargin>0 && ischar(varargin{1}) && strcmpi(varargin{1},'vfiles'),
    t = vfiles('all');
    return;
elseif nargin>0 && ischar(varargin{1}) && strcmpi(varargin{1},'cpath'),
    t = cpath(varargin{2:end});
else
    [t,sts] = selector(varargin{:});
end;
return;
%=======================================================================

%=======================================================================
function [t,ok] = selector(n,typ,mesg,already,wd,varargin)
if nargin<5, wd   = pwd; end;
if nargin<4, already = {''}; end;
if nargin<3, mesg = 'Select files...'; end;
if nargin<2, typ  = 'any'; end;
if nargin<1, n    = [0 Inf]; end;
ok  = 0;
if numel(n)==1,   n    = [n n];    end;
if n(1)>n(2),     n    = n([2 1]); end;
if ~finite(n(1)), n(1) = 0;        end;
if numel(already)>n(2), already = already(1:n(2)); end
already = strvcat(already);

t = '';
switch lower(typ),
case {'any'},   code = 0; ext = {'.*'};
case {'image'}, code = 1; ext = {'.*\.nii$','.*\.img$','.*\.NII$','.*\.IMG$'};
case {'xml'},   code = 0; ext = {'.*\.xml$','.*\.XML$'};
case {'mat'},   code = 0; ext = {'.*\.mat$','.*\.MAT$'};
case {'batch'}, code = 0; ext = {'.*\.mat$','.*\.MAT$','.*\.xml$','.*\.XML$'};
case {'dir'},   code =-1; ext = {'.*'};
otherwise,      code = 0; ext = {typ};
end;

[col1,col2,col3,fs] = colours;

fg = figure('IntegerHandle','off',...
        'Tag','Select',...
        'Name',strvcat(mesg),...
        'NumberTitle','off',...
        'Color',[1 1 1],...
        'Units','Pixels',...
        'MenuBar','none',...
        'DefaultTextInterpreter','none',...
        'DefaultUicontrolInterruptible','on',...
        'ResizeFcn',@resize_fun,...
        'KeyPressFcn',@hitkey);

fh = 0.05;
%fs = 10;

h1 = (0.96-4*fh-5*0.01)/2;
if n(2)*fh<h1,
    h1 = n(2)*fh;
end;
h2 = 0.96-4*fh-5*0.01-h1;

SPMdir = fileparts(which(mfilename));
if str2double(version('-release'))==14 && isdeployed,
    ind = findstr(SPMdir,'_mcr')-1;
    [SPMdir,junk] = fileparts(SPMdir(1:ind(1)));
end;
prevdirs([SPMdir filesep]);
[pd,vl] = prevdirs([wd filesep]);

% Selected Files
hp = 0.02;
sel = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.02 hp 0.96 h1],...
    'FontSize',fs,...
    'Callback',@unselect,...
    'tag','selected',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'Max',10000,...
    'Min',0,...
    'String',already,...
    'Value',1);
c0 = uicontextmenu('Parent',fg);
set(sel,'uicontextmenu',c0);
uimenu('Label','Unselect All', 'Parent',c0,'Callback',@unselect_all);

% Messages
hp = hp+h1+0.01;
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 hp 0.96 fh],...
    'FontSize',fs,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'HorizontalAlignment','left',...
    'Tag','msg',...
    'String',mesg);

if strcmpi(typ,'image'),
    uicontrol(fg,...
        'style','edit',...
        'units','normalized',...
        'Position',[0.61 hp 0.37 fh],...
        'Callback',@update_frames,...
        'tag','frame',...
        'FontSize',fs,...
        'BackgroundColor',col1,...
        'String','1','UserData',1);
% 'ForegroundGolor',col3,...
end;

% Help
hp = hp+fh+0.01;
uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.02 hp fh fh],...
    'FontSize',fs,...
    'Callback',@heelp,...
    'tag','?',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','?',...
    'FontWeight','bold',...
    'ToolTipString','Show Help',...
    'FontSize',fs);

uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.03+fh hp fh fh],...
    'FontSize',fs,...
    'Callback',@editwin,...
    'tag','Ed',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Ed',...
    'FontWeight','bold',...
    'ToolTipString','Edit Selected Files',...
    'FontSize',fs);

% Done
dne = uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.04+2*fh hp 0.45-2*fh fh],...
    'FontSize',fs,...
    'Callback',@delete,...
    'tag','D',...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'String','Done',...
    'FontWeight','bold',...
    'FontSize',fs,...
    'Enable','off',...
    'DeleteFcn',@null);

if size(already,1)>=n(1) && size(already,1)<=n(2),
    set(dne,'Enable','on');
end;

% Filter Button
uicontrol(fg,...
    'Style','pushbutton',...
    'units','normalized',...
    'Position',[0.51 hp 0.1 fh],...
    'FontSize',fs,...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'Callback',@clearfilt,...
    'String','Filt',...
    'FontSize',fs);

% Filter
ud     = struct('ext',{ext},'code',code);
uicontrol(fg,...
    'style','edit',...
    'units','normalized',...
    'Position',[0.61 hp 0.37 fh],...
    'ForegroundColor',col3,...
    'BackgroundColor',col1,...
    'FontSize',fs,...
    'Callback',@update,...
    'tag','regexp',...
    'String','.*',...
    'UserData',ud);

% Directories
hp = hp + fh+0.01;
db = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.02 hp 0.47 h2],...
    'FontSize',fs,...
    'Callback',@click_dir_box,...
    'tag','dirs',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'Max',1,...
    'Min',0,...
    'String','',...
    'UserData',wd,...
    'Value',1);

% Files
tmp = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.51 hp 0.47 h2],...
    'FontSize',fs,...
    'Callback',@click_file_box,...
    'tag','files',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'UserData',n,...
    'Max',10240,...
    'Min',0,...
    'String','',...
    'Value',1);
c0 = uicontextmenu('Parent',fg);
set(tmp,'uicontextmenu',c0);
uimenu('Label','Select All', 'Parent',c0,'Callback',@select_all);

% Drives
if strcmpi(computer,'PCWIN'),
    dr  = spm_win32utils('drives');
    drivestr = cell(1,numel(dr));
    for i=1:numel(dr),
        drivestr{i} = [dr(i) ':'];
    end;
    %drivestr = {'A:','B:','C:','D:'};
    sz = get(db,'Position');
    sz(4) = sz(4)-fh-2*0.01;
    set(db,'Position',sz);
    uicontrol(fg,...
        'style','text',...
        'units','normalized',...
        'Position',[0.02 hp+h2-fh-0.01 0.10 fh],...
        'FontSize',fs,...
        'BackgroundColor',get(fg,'Color'),...
        'ForegroundColor',col3,...
        'String','Drive');
    uicontrol(fg,...
        'style','popupmenu',...
        'units','normalized',...
        'Position',[0.12 hp+h2-fh-0.01 0.37 fh],...
        'FontSize',fs,...
        'Callback',@setdrive,...
        'tag','drive',...
        'BackgroundColor',col1,...
        'ForegroundColor',col3,...
        'String',drivestr,...
        'Value',1);
end;

% Previous dirs
hp = hp+h2+0.01;
uicontrol(fg,...
    'style','popupmenu',...
    'units','normalized',...
    'Position',[0.12 hp 0.86 fh],...
    'FontSize',fs,...
    'Callback',@click_dir_list,...
    'tag','previous',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String',pd,...
    'Value',vl);
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 hp 0.10 fh],...
    'FontSize',fs,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Prev');

% Directory
hp = hp + fh+0.01;
uicontrol(fg,...
    'style','edit',...
    'units','normalized',...
    'Position',[0.12 hp 0.86 fh],...
    'FontSize',fs,...
    'Callback',@edit_dir,...
    'tag','edit',...
    'BackgroundColor',col1,...
    'ForegroundColor',col3,...
    'String','');
uicontrol(fg,...
    'style','text',...
    'units','normalized',...
    'Position',[0.02 hp 0.10 fh],...
    'FontSize',fs,...
    'BackgroundColor',get(fg,'Color'),...
    'ForegroundColor',col3,...
    'String','Dir');

resize_fun(fg);
update(sel,wd)

waitfor(dne);
if ishandle(sel),
    t  = get(sel,'String');
    ok = 1;
end;
if ishandle(fg),  delete(fg); end;
return;
%=======================================================================

%=======================================================================
function null(varargin)
%=======================================================================

%=======================================================================
function msg(ob,str)
ob = sib(ob,'msg');
set(ob,'String',str);
if nargin>=3,
    set(ob,'ForegroundColor',[1 0 0],'FontWeight','bold');
else
    set(ob,'ForegroundColor',[0 0 0],'FontWeight','normal');
end;
drawnow;
return;
%=======================================================================

%=======================================================================
function setdrive(ob,varargin)
st = get(ob,'String');
vl = get(ob,'Value');
update(ob,st{vl});
return;
%=======================================================================

%=======================================================================
function resize_fun(fg,varargin)
ob = findobj(fg,'String','Filt','Style','pushbutton');
if ~isempty(ob),
    ofs = get(ob,'FontSize');
    ex = get(ob,'Extent');
    ps = get(ob,'Position');
    fs = floor(ofs*min(ps(4)./ex(4))+1);
    fs = max(min(fs,30),4);
    ob = findobj(fg,'Fontsize',ofs);
    set(ob,'FontSize',fs);
end;
return;
%=======================================================================

%=======================================================================
function [d,mch] = prevdirs(d)
persistent pd
if ~iscell(pd), pd = {}; end;
d   = deblank(d);
mch = find(strcmp(d,pd));
if isempty(mch),
    pd  = {pd{:},d};
    mch = length(pd);
end;
d = pd;
return;
%=======================================================================

%=======================================================================
function clearfilt(ob,varargin)
set(sib(ob,'regexp'),'String','.*');
update(ob);
return;
%=======================================================================

%=======================================================================
function click_dir_list(ob,varargin)
vl = get(ob,'Value');
ls = get(ob,'String');
update(ob,deblank(ls{vl}));
return;
%=======================================================================

%=======================================================================
function edit_dir(ob,varargin)
update(ob,get(ob,'String'));
return;
%=======================================================================

%=======================================================================
function click_dir_box(lb,varargin)
update(lb,current_dir(lb));
return;
%=======================================================================

%=======================================================================
function dr = current_dir(lb,varargin)
pd  = get(sib(lb,'edit'),'String');
vl  = get(lb,'Value');
str = get(lb,'String');
dr = fullfile(pd,deblank(str(vl,:)));
if vl==1, % Current directory
    dr = fileparts(dr);
end;
if vl==2, % Parent directory
    dr = fileparts(dr);
    dr = fileparts(dr);
end;
return;
%=======================================================================

%=======================================================================
function re = getfilt(ob)
ob  = sib(ob,'regexp');
ud  = get(ob,'UserData');
re  = struct('code',ud.code,...
             'frames',get(sib(ob,'frame'),'UserData'),...
             'ext',{ud.ext},...
             'filt',get(sib(ob,'regexp'),'String'));
return;
%=======================================================================

%=======================================================================
function update(lb,dr)
lb = sib(lb,'dirs');
if nargin<2 || isempty(dr),
    dr = get(lb,'UserData');
end;
if ~strcmpi(computer,'PCWIN')
    dr    = [filesep dr filesep];
else
    dr    = [dr filesep];
end;
dr(findstr([filesep filesep],dr)) = [];
[f,d] = listfiles(dr,getfilt(lb));
if isempty(d),
    dr    = get(lb,'UserData');
    [f,d] = listfiles(dr,getfilt(lb));
else
    set(lb,'UserData',dr);
end;
set(lb,'Value',1,'String',d);
set(sib(lb,'files'),'Value',1,'String',f);
[ls,mch] = prevdirs(dr);
set(sib(lb,'previous'),'String',ls,'Value',mch);
set(sib(lb,'edit'),'String',dr);

if numel(dr)>1 && dr(2)==':',
    str = get(sib(lb,'drive'),'String');
    str = cat(1,char(str));
    mch = find(lower(str(:,1))==lower(dr(1)));
    if ~isempty(mch),
        set(sib(lb,'drive'),'Value',mch);
    end;
end;
return;
%=======================================================================

%=======================================================================
function update_frames(lb,varargin)
str = get(lb,'String');
%r   = get(lb,'UserData');
try
    r = eval(['[',str,']']);
catch
    msg(lb,['Failed to evaluate "' str '".'],'r');
    beep;
    return;
end;
if ~isnumeric(r),
    msg(lb,['Expression non-numeric "' str '".'],'r');
    beep;
else
    set(lb,'UserData',r);
    msg(lb,'');
    update(lb);
end;
%=======================================================================

%=======================================================================
function select_all(ob,varargin)
lb = findobj(get(get(ob,'Parent'),'Parent'),'Tag','files');
str  = get(lb,'String');
set(lb,'Value',1:size(str,1));
drawnow;
click_file_box(lb);
return;
%=======================================================================

%=======================================================================
function click_file_box(lb,varargin)
lim  = get(lb,'UserData');
ob   = sib(lb,'selected');
str3 = get(ob,'String');

str  = get(lb,'String');
vlo  = get(lb,'Value');
lim1  = min(max(lim(2)-size(str3,1),0),length(vlo));
if isempty(vlo),
    msg(lb,'Nothing selected');
    return;
end;
if lim1==0,
    msg(lb,['Selected ' num2str(size(str3,1)) '/' num2str(lim(2)) ' already.']);
    beep;
    set(sib(lb,'D'),'Enable','on');
    return;
end;

vl   = vlo(1:lim1);
msk  = false(size(str,1),1);
if vl>0, msk(vl) = true; else msk = []; end;
str1 = str( msk,:);
str2 = str(~msk,:);
dr   = [current_dir(sib(lb,'dirs')) filesep];
str1 = [repmat(dr,size(str1,1),1) str1];

set(lb,'Value',min(vl(1),size(str2,1)),'String',str2);
r    = (1:size(str1,1))+size(str3,1);
str3 = deblank(strvcat(str3,str1));
set(ob,'String',str3,'Value',r);
if length(vlo)>lim1,
    msg(lb,['Retained ' num2str(lim1) '/' num2str(length(vlo))...
        ' of selection.']);
    beep;
elseif finite(lim(2))
    if lim(1)==lim(2),
        msg(lb,['Selected ' num2str(size(str3,1)) '/' num2str(lim(2)) ' files.']);
    else
        msg(lb,['Selected ' num2str(size(str3,1)) '/' num2str(lim(1)) '-' num2str(lim(2)) ' files.']);
    end;
else
    if size(str3,1) == 1, ss = ''; else ss = 's'; end;
    msg(lb,['Selected ' num2str(size(str3,1)) ' file' ss '.']);
end;
if ~finite(lim(1)) || size(str3,1)>=lim(1),
    set(sib(lb,'D'),'Enable','on');
end;

return;
%=======================================================================

%=======================================================================
function obj = sib(ob,tag)
obj = findobj(get(ob,'Parent'),'Tag',tag);
return;
%if isempty(obj),
%    error(['Can''t find object with tag "' tag '".']);
%elseif length(obj)>1,
%    error(['Found ' num2str(length(obj)) ' objects with tag "' tag '".']);
%end;
%return;
%=======================================================================

%=======================================================================
function unselect(lb,varargin)
vl      = get(lb,'Value');
str     = get(lb,'String');
msk     = ones(size(str,1),1);
if vl~=0, msk(vl) = 0; end;
str2    = str(logical(msk),:);
set(lb,'Value',min(vl(1),size(str2,1)),'String',str2);
lim = get(sib(lb,'files'),'UserData');
if size(str2,1)>= lim(1) && size(str2,1)<= lim(2),
    set(sib(lb,'D'),'Enable','on');
else 
    set(sib(lb,'D'),'Enable','off');
end;

if size(str2,1) == 1, ss = ''; else ss = 's'; end;
msg(lb,['Unselected ' num2str(size(str2,1)) ' file' ss '.']);
return;
%=======================================================================

%=======================================================================
function unselect_all(ob,varargin)
lb = findobj(get(get(ob,'Parent'),'Parent'),'Tag','selected');
set(lb,'Value',[],'String','','ListBoxTop',1);
msg(lb,'Unselected all files.');
lim = get(sib(lb,'files'),'UserData');
if lim(1)>0, set(sib(lb,'D'),'Enable','off'); end;
return;
%=======================================================================

%=======================================================================
function varargout = vfiles(option,varargin)
persistent vfs
if isempty(vfs),
    vfs = newvfs;
end;

switch option,
case {'clear'}
    vfs = newvfs;
case {'add'}
    for j=1:numel(varargin),
        if ischar(varargin{j}),
            for i=1:size(varargin{j},1),
                fle = deblank(varargin{j}(i,:));
                vfs = addvfile(vfs,fle);
            end;
        elseif iscell(varargin{j}),
            for i=1:numel(varargin{j}),
                fle = deblank(varargin{j}{i});
                vfs = addvfile(vfs,fle);
            end;
        end;
    end;
case {'list'}
    [varargout{1:3}] = listvfiles(vfs,varargin{:});
case {'all'}
    varargout{1} = vfs;
otherwise
    error('Unknown option.');
end;
return;
%=======================================================================

%=======================================================================
function vfs = newvfs(nam)
if nargin==0, nam = ''; end;
vfs = struct('name',nam,'dirs',struct('name',{},'dirs',{},'files',{}),'files',struct('name',{},'ind',{}));
return;
%=======================================================================

%=======================================================================
function vfs = addvfile(vfs,fle)
ind = find(fle==filesep);
if any(ind==1),
    ind = ind(2:end)-1;
    fle = fle(2:end);
end;
if isempty(ind),
    [unused,nam,ext,num] = spm_fileparts(fle);
    if ~isempty(num),
        ind = [str2num(num) 1 1];
        ind = ind(1);
    else
        ind = [];
    end;
    fname = [nam ext];
    mch   = strcmp(fname,{vfs.files.name});
    if any(mch),
        mch                = find(mch);
        vfs.files(mch).ind = [vfs.files(mch).ind ind];
    else
        vfs.files(end+1).name = fname;
        vfs.files(end).ind    = ind;
    end;
else
    dr   = fle(1:(ind(1)-1));
    fle  = fle((ind(1)+1):end);
    mch  = strcmp(dr,{vfs.dirs.name});
    if any(mch),
        mch           = find(mch);
    else
        mch           = numel(vfs.dirs)+1;
        vfs.dirs(mch) = newvfs(dr);
    end;
    vfs.dirs(mch)     = addvfile(vfs.dirs(mch),fle);
end;
return;
%=======================================================================

%=======================================================================
function [f,d] = listfiles(dr,filt)
ob = gco;
msg(ob,'Listing directory...');
if nargin<2, filt = '';  end;
if nargin<1, dr   = '.'; end;
de      = dir(dr);
if ~isempty(de),
    d       = {de(find( cell2mat({de.isdir}))).name};
    if filt.code~=-1,
        f   = {de(find(~cell2mat({de.isdir}))).name};
    else
        % f = d(3:end);
        f   = d;
    end;
else
    d = {'.','..'};
    f = {};
end;

msg(ob,['Filtering ' num2str(numel(f)) ' files...']);
f  = do_filter(f,filt.ext);
f  = do_filter(f,{filt.filt});
ii = cell(1,numel(f));
if filt.code==1 && (numel(filt.frames)~=1 || filt.frames(1)~=1),
    msg(ob,['Reading headers of ' num2str(numel(f)) ' images...']);
    for i=1:numel(f),
        try
            ni = nifti(fullfile(dr,f{i}));
            dm = [ni.dat.dim 1 1 1 1 1];
            d4 = (1:dm(4))';
        catch
            d4 = 1;
        end;
        msk = false(size(filt.frames));
        for j=1:numel(msk), msk(j) = any(d4==filt.frames(j)); end;
        ii{i} = filt.frames(msk);
    end;
elseif filt.code==1 && (numel(filt.frames)==1 && filt.frames(1)==1),
    for i=1:numel(f),
        ii{i} = 1;
    end;
end;

msg(ob,'Listing virtual files...');
[fv,dv,iv] = vfiles('list',dr);
if filt.code==-1,
    fv = dv;
    iv = cell(size(fv));
end;
msg(ob,['Filtering ' num2str(numel(fv)) ' virtual files...']);
[fv,ind]   = do_filter(fv,filt.ext);
iv         = iv(ind);
[fv,ind]   = do_filter(fv,{filt.filt});
iv         = iv(ind);
if filt.code==1,
    for i=1:numel(iv),
        msk   = false(size(filt.frames));
        for j=1:numel(msk), msk(j) = any(iv{i}==filt.frames(j)); end;
        iv{i} = filt.frames(msk);
    end;
end;

d       = { d{:},dv{:}};
f       = { f{:},fv{:}};
ii      = {ii{:},iv{:}};

msg(ob,['Listing ' num2str(numel(f)) ' files...']);

[f,ind] = sortrows(f(:));
ii      = ii(ind);
msk     = true(1,numel(f));
if ~isempty(f), f{1} = deblank(f{1}); end;
for i=2:numel(f),
    f{i} = deblank(f{i});
    if strcmp(f{i-1},f{i}),
        if filt.code==1,
            tmp      = sort([ii{i}(:) ; ii{i-1}(:)]);
            tmp(~diff(tmp,1)) = [];
            ii{i}    = tmp;
        end;
        msk(i-1) = false;
    end;
end;
f        = f(msk);
if filt.code==1,
    ii       = ii(msk);
    c        = cell(size(f));
    for i=1:numel(f),
        c{i} = [repmat([f{i} ','],numel(ii{i}),1) num2str(ii{i}(:)) ];
    end;
    f        = strvcat(c{:});
elseif filt.code==-1,
    fs = filesep;
    for i=1:numel(f),
        f{i} = [f{i} fs];
    end;
    f        = strvcat(f{:});
else
    f        = strvcat(f{:});
end;

d        = sortrows(d(:));
d        = strvcat(d);
sam      = find(~any(diff(d+0,1),2));
d(sam,:) = [];
msg(ob,'');
return;
%=======================================================================

%=======================================================================
function [f,ind] = do_filter(f,filt)
t2 = false(numel(f),1);
for j=1:numel(filt),
    t1 = regexp(f,filt{j});
    if numel(f)==1, t1 = {t1}; end;
    for i=1:numel(t1),
        t2(i) = t2(i) || ~isempty(t1{i});
    end;
end;
ind = find(t2);
f   = f(ind);
return;
%=======================================================================

%=======================================================================
function [f,d,ii] = listvfiles(vfs,dr)
f  = {};
d  = {};
ii = {};
if isempty(dr),
    f  = {vfs.files.name};
    ii = {vfs.files.ind};
    d  = {vfs.dirs.name};
else
    if dr(1)==filesep, dr = dr(2:end); end;
    ind = find(dr==filesep);
    if isempty(ind),
        d1 = dr;
        d2 = '';
    else
        d1 = dr(1:(ind(1)-1));
        d2 = dr((ind(1)+1):end);
    end;
    for i=1:length(vfs.dirs),
        if strcmp(d1,vfs.dirs(i).name),
            [f,d,ii] = listvfiles(vfs.dirs(i),d2);
            break;
        end;
    end;
end;
return;
%=======================================================================

%=======================================================================
function heelp(ob,varargin)
[col1,col2,col3,fs] = colours;
fg = get(ob,'Parent');
t  = uicontrol(fg,...
    'style','listbox',...
    'units','normalized',...
    'Position',[0.01 0.01 0.98 0.98],...
    'FontSize',fs,...
    'FontName','FixedWidthFont',...
    'BackgroundColor',col2,...
    'ForegroundColor',col3,...
    'Max',0,...
    'Min',0,...
    'tag','HelpWin',...
    'String','                   ');
c0 = uicontextmenu('Parent',fg);
set(t,'uicontextmenu',c0);
uimenu('Label','Done', 'Parent',c0,'Callback',@helpclear);

ext = get(t,'Extent');
pw  = floor(0.98/ext(3)*20-4);
str  = spm_justify(pw,{[...
'File Selection help. You can return to selecting files via the right mouse button (the "Done" option). ',...
'Because of a bug in Matlab (on some machines), don''t resize this window when viewing the help.'],...
'',[...
'The panel at the bottom shows files that are already selected. ',...
'Clicking a selected file will un-select it. To un-select several, you can ',...
'drag the cursor over the files, and they will be gone on release. ',...
'You can use the right mouse button to un-select everything.'],...
'',[...
'Directories are navigated by editing the name of the current directory (where it says "Dir"), ',...
'by going to one of the previously entered directories ("Prev"), or by navigating around ',...
'the parent or subdirectories listed in the left side panel.'],...
'',[...
'Files matching the filter ("Filt") are shown in the panel on the right. ',...
'These can be selected by clicking or dragging.  Use the right mouse button if ',...
'you would like to select all files.  Note that when selected, the files disappear ',...
'from this panel.  They can be made to reappear by re-specifying the directory ',...
'or the filter. ',...
'Note that the syntax of the filter differs from that used by previous versions of ',...
'SPM.  The following is a list of symbols with special meaning for filtering the filenames:'],...
'    ^     start of string',...
'    $     end of string',...
'    .     any character',...
'    \     quote next character',...
'    *     match zero or more',...
'    +     match one or more',...
'    ?     match zero or one, or match minimally',...
'    {}    match a range of occurrances',...
'    []    set of characters',...
'    [^]   exclude a set of characters',...
'    ()    group subexpression',...
'    \w    match word [a-z_A-Z0-9]',...
'    \W    not a word [^a-z_A-Z0-9]',...
'    \d    match digit [0-9]',...
'    \D    not a digit [^0-9]',...
'    \s    match white space [ \t\r\n\f]',...
'    \S    not a white space [^ \t\r\n\f]',...
'    \<WORD\>    exact word match',...
'',[...
'Individual time frames of image files can also be selected.  The frame filter ',...
'allows specified frames to be shown, which is useful for image files that ',...
'contain multiple time points.  If your images are only single time point, then ',...
'reading all the image headers can be avoided by specifying a frame filter of "1". ',...
'The filter should contain a list of integers indicating the frames to be used. ',...
'This can be generated by e.g. "1:100", or "1:2:100".'],...
'',[...
'There is also an edit button (Ed), which allows you to edit your selection of files. ',...
'When you are done, then use the menu-button of your mouse to either cancel or accept your changes'],''});
pad = cellstr(char(zeros(max(0,floor(1.2/ext(4) - numel(str))),1)));
str = {str{:}, pad{:}};
set(t,'String',str);
return;
%=======================================================================

%=======================================================================
function helpclear(ob,varargin)
ob = get(ob,'Parent');
ob = get(ob,'Parent');
ob = findobj(ob,'Tag','HelpWin');
delete(ob);
%=======================================================================

%=======================================================================
function hitkey(fg,varargin)
ch = get(fg,'CurrentCharacter');
if isempty(ch), return; end;

ob = findobj(fg,'Tag','files');
if ~isempty(ob),
    f = get(ob,'String');
    f = f(:,1);
    fset = find(f>=ch);
    if ~isempty(fset),
        fset = fset(1);
        %cb = get(ob,'Callback');
        %set(ob,'Callback',[]);
        set(ob,'ListboxTop',fset);
        %set(ob,'Callback',cb);
    else
        set(ob,'ListboxTop',length(f));
    end;
end;
return;
%=======================================================================

%=======================================================================
function t = cpath(t,d)

switch spm_platform('filesys'),
case 'unx',
    mch = '^/';
    fs  = '/';
    fs1 = '/';
case 'win',
    mch = '^.:\\';
    fs  = '\';
    fs1 = '\\';
otherwise;
    error('What is this filesystem?');
end

if isempty(regexp(t,mch,'once')),
    if nargin<2, d = pwd; end;
    t = [d fs t];
end;

% Replace (most) occurences of '/./' by '/' (problems with e.g. /././././././')
for i=1:5, % Should really repeat until no more substitutions are made
    t  = regexprep(t,[fs1 '\.' fs1], fs);
end;
t  = regexprep(t,[fs1 '\.' '$'], fs);

% Replace (most) occurences of '/abc/../' by '/'
for i=1:5, % Should really repeat until no more substitutions are made
    t  = regexprep(t,[fs1 '[^' fs1 ']+' fs1 '\.\.' fs1],fs,'once');
end;
t  = regexprep(t,[fs1 '[^' fs1 ']+' fs1 '\.\.' '$'],fs,'once');

% Replace '//'
t  = regexprep(t,[fs1 '+'], fs);
%=======================================================================

%=======================================================================
function editwin(ob,varargin)
[col1,col2,col3,fs] = colours;
fg   = get(ob,'Parent');
lb   = findobj(fg,'Tag','selected');
str  = get(lb,'String');
str  = cellstr(str);
h    = uicontrol(fg,'Style','Edit',...
        'units','normalized',...
        'String',str,...
        'FontSize',16,...
        'Max',2,...
        'Tag','EditWindow',...
        'HorizontalAlignment','Left',...
        'ForegroundColor',col3,...
        'BackgroundColor',col1,...
        'Position',[0.01 0.01 0.98 0.98]);
c0 = uicontextmenu('Parent',fg);
set(h,'uicontextmenu',c0);
uimenu('Label','Cancel', 'Parent',c0,'Callback',@editclear);
uimenu('Label','Accept', 'Parent',c0,'Callback',@editdone);
%=======================================================================

%=======================================================================
function editdone(ob,varargin)
ob  = get(ob,'Parent');
fg  = get(ob,'Parent');
ob  = findobj(fg,'Tag','EditWindow');
str = deblank(get(ob,'String'));
str = strvcat(cellstr(str));
lb  = findobj(fg,'Tag','selected');
set(lb,'String',str,'Value',[]);
delete(ob);
%=======================================================================

%=======================================================================
function editclear(ob,varargin)
ob = get(ob,'Parent');
ob = get(ob,'Parent');
ob = findobj(ob,'Tag','EditWindow');
delete(ob);
%=======================================================================

%=======================================================================
function [c1,c2,c3,fs] = colours
global defaults
c1 = [1 1 1];
c2 = [1 1 1];
c3 = [0 0 0];
fs = 14;
if isfield(defaults,'ui'),
    ui = defaults.ui;
    if isfield(ui,'colour1'), c1 = ui.colour1; end;
    if isfield(ui,'colour2'), c2 = ui.colour2; end;
    if isfield(ui,'colour3'), c3 = ui.colour3; end;
    if isfield(ui,'fs'),      fs = ui.fs;      end;
end;
%=======================================================================

%=======================================================================

