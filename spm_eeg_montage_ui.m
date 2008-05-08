function montage = spm_eeg_montage_ui(montage)
%
% A montage is specified as a structure with the fields
%   montage.tra      = MxN matrix
%   montage.labelnew = Mx1 cell-array
%   montage.labelorg = Nx1 cell-array
%
% This function builds a GUI for montage editing.
% This is used for:
%   - visualization purposes
%   - classification purposes (not in use)
% IN:
%   - D: the user-data structure
% OUT:
%   - montage: a structure containing:
%       .M :the mxp matrix parametrizing the montage, where m is the number
%       of channels, and p is the number of sensors
%       .clab: a mx1 cell array of strings containing the labelso f the
%       channels (combinations of the labels of the sensors)

dbstop if error

if nargin < 1 || ~isfield(montage, 'labelorg') || isempty(montage.labelorg)
    help spm_eeg_montage_ui
    error('Insufficient or illegal input');
end

if ~isfield(montage, 'labelnew')
    montage.labelnew = montage.labelorg;
end

ns = length(montage.labelorg);
nc = length(montage.labelnew);

if ~isfield(montage, 'tra')
    montage.tra = sparse(zeros(nc, ns));
    [sel1 sel2] = spm_match_str(montage.labelnew, montage.labelorg);
    montage.tra(sub2ind([nc, ns], sel1, sel2)) = 1;
end


disp('______________________')
dispSTR = ['Current montage : ',num2str(nc),'original channels for ',num2str(ns),'new channels'];
disp(dispSTR)
fprintf('Building montage GUI ...')

hf = figure('menubar','none','numbertitle','off','name','Montage editing: building GUI...');
[posfig,pointers] = buildGUIperSe(montage, [], hf);
set(hf,'name','Montage editing');
pointers.hf = hf;

forceMontage(pointers, montage.tra);

handles.pointers = pointers;
set(hf,'userdata', handles);


waitfor(hf);

return

[M, str] = evalMontage(pointers);



pointers.handles = D.PSD.handles;
udb.M = D.PSD.VIZU.montage.M;
udb.str = str;
udb.pointers = pointers;
set(pointers.hf,'userdata',udb)

% ud.montage.pointers = pointers;

drawnow
fprintf(1,' OK.')
fprintf(1,'\n')
disp('______________________')




function doOK(ok_btn,evd)

dbstop if error
h = get(ok_btn,'parent');
hf=(get(gcbo,'userdata'));
udb = get(hf,'userdata');
[M,str] = feval(@evalMontage,udb.pointers);
udb.M = M;
udb.str = str;
% build montage structure variable
sM = sum(abs(M),2);
emptyChans = find(sM==0);
if isempty(emptyChans) ~= 1
    M(emptyChans,:) = [];
    fullChans = find(sM~=0);
    for i = 1:length(fullChans)
        clab{i} = str{fullChans(i)};
    end
else
    clab = str;
end

return



function doCheck(ok_btn,evd)
% reads edited montage and updates channel config texts
dbstop if error
hf = get(ok_btn,'parent');
udb = get(hf,'userdata');
[M,str] = feval(@evalMontage,udb.pointers);
udb.M = M;
udb.str = str;
set(hf,'userdata',udb)
return

function doAddChannel(ok_btn,evd)
% adds a line to the GUI
hf = get(ok_btn,'parent')
udb = get(hf,'userdata');
ud = get(udb.pointers.handles.axes,'userdata');
handles = udb.pointers.handles;
% add channel to relevant variables
udb.M = [udb.M;zeros(1,size(udb.M,2))];
udb.pointers.hce = [udb.pointers.hce;zeros(1,size(udb.M,2))];
udb.pointers.hct = [udb.pointers.hct;0];
udb.pointers.hct2 = [udb.pointers.hct2;0];
ud.montage.M = udb.M;
ud.handles.pointers = udb.pointers;
[posfig,udb.pointers] = feval(@buildGUIperSe,ud,hf);
feval(@forceMontage,udb.pointers,udb.M);
[udb.M,udb.str] = feval(@evalMontage,udb.pointers);
udb.pointers.handles = handles;
udb.pointers.hf = hf;
set(udb.pointers.hf,'userdata',udb)
return

function doSave(ok_btn,evd)
% save edited montage
dbstop if error
hf = get(ok_btn,'parent');
udb = get(hf,'userdata');
[M,str] = feval(@evalMontage,udb.pointers);
uisave('M','PSD_montage.mat');


function doLoad(ok_btn,evd)
% load montage and fill up GUI
dbstop if error
hf = get(ok_btn,'parent');
udb = get(hf,'userdata');
[FileName,PathName] = uigetfile('*.mat','Load montage file');
montage = load(fullfile(PathName,FileName));
name = fieldnames(montage);
try
    montage = getfield(montage, name{1});
    M = montage.tra;
    p = size(udb.pointers.hce,1);
    n = size(udb.pointers.hce,2);
    if size(M,1) < p
        ds = p - size(M,1);
        M = [M;zeros(ds,n)];
    elseif size(M,1) > p
        ud = get(udb.pointers.handles.axes,'userdata');
        handles = udb.pointers.handles;
        ud.montage.M = M;
        % create pointers for building GUI
        udb.pointers.hce = [udb.pointers.hce;zeros(size(M,1)-p,size(udb.M,2))];
        udb.pointers.hct = [udb.pointers.hct;zeros(size(M,1)-p,1)];
        udb.pointers.hct2 = [udb.pointers.hct2;zeros(size(M,1)-p,1)];
        ud.handles.pointers = udb.pointers;
        [posfig,udb.pointers] = feval(@buildGUIperSe, montage, ud, hf);
        udb.pointers.handles = handles;
        udb.pointers.hf = hf;
    end
    feval(@forceMontage,udb.pointers, M);
    [udb.M,udb.str] = feval(@evalMontage,udb.pointers);
    set(hf,'userdata',udb)
catch
    errordlg('Problem with the montage file');
end
return



function forceMontage(pointers,M)
% fill the editable boxes with montage M
p = size(M,1);
n = size(M,2);
for i = 1:p
    for j = 1:n
        if isequal(M(i,j),0) ~= 1
            set(pointers.hce(i,j),'string',num2str(M(i,j)));
        else
            set(pointers.hce(i,j),'string',[]);
        end
    end
end
return


function [M,str] = evalMontage(pointers)
% reads the editable boxes and return montage
dbstop if error
p = size(pointers.hce,1);
n = size(pointers.hce,2);
M = zeros(p,n);
for j = 1:p
    str{j} = [];
    non0 = 0;
    for i = 1:n
        tmp = get(pointers.hce(j,i),'string');
        if isempty(tmp) ~= 1
            M(j,i) = str2num(tmp);
            if M(j,i) ~= 0
                clab = get(pointers.hcc(i),'string');
                switch M(j,i)
                    case -1
                        str{j} = [str{j},' -',clab];
                    case 1
                        str{j} = [str{j},' +',clab];
                    otherwise
                        str{j} = [str{j},' ',tmp,'x ',clab];
                end
                non0 = non0 + 1;
            end
        end
    end
    if isempty(str{j}) ~= 1
        if isequal(str{j}(1:2),' +') == 1
            str{j}(1:2) = [];
        elseif isequal(str{j}(1),' ') == 1
            str{j}(1) = [];
        end
    end
    if  non0 ~= 0
        set(pointers.hct2(j),'string',str{j},...
            'tooltipstring',['channel ',num2str(j),': ',str{j}]);
        set(pointers.hct(j),'string',num2str(j));
    else
        set(pointers.hct2(j),'string',[],'tooltipstring',['channel ',num2str(j),': EMPTY CHANNEL']);
        set(pointers.hct(j),'string',[],'tooltipstring',['channel ',num2str(j),': EMPTY CHANNEL']);
    end
end
return




function [posfig,pointers] = buildGUIperSe(montage, handles, hf)
% creates the editable boxes from montage
p = size(montage.tra, 1);
n = size(montage.tra, 2);
mase = 20;              % max size of editable boxes
% fix the size of the GUI
scrsz = get(0,'ScreenSize');
H = 9*scrsz(4)/11;
se = (H-60)./p;
se = min([mase se]);
H = se.*p + 100;
L = se.*n + 100;
posfig = [scrsz(3)-max([L+40,540])-40 scrsz(4)/10 max([L+40,540]) max([H,500])];
hc = get(hf,'children');
if isfield(handles,'pointers')
    oldPos = get(hf,'position');
    posfig(1:2) = oldPos(1:2);
    hcc = ud.handles.pointers.hcc;
    hce = ud.handles.pointers.hce;
    hct = ud.handles.pointers.hct;
    hct2 = ud.handles.pointers.hct2;
    ud.handles.pointers
else
    hcc = zeros(n,1);       % sensor names
    hce = zeros(p,n);       % editable boxes
    hct = zeros(p,1);       % channel numbers
    hct2 = zeros(p,1);      % channel config
    % add buttons
    hch = uicontrol('style','pushbutton','position',[posfig(3)-100 40 90 20],...
        'string','Check montage','userdata',hf,'callback',{@doCheck},...
        'visible','off');
    hlo = uicontrol('style','pushbutton','position',[110 40 90 20],...
        'string','Load montage','userdata',hf,'callback',{@doLoad},...
        'visible','off');
    hsa = uicontrol('style','pushbutton','position',[10 40 90 20],...
        'string','Save montage','userdata',hf,'callback',{@doSave},...
        'visible','off');
    hco = uicontrol('style','pushbutton','position',[posfig(3)-200 40 90 20],...
        'string','Add channel','userdata',hf,'callback',{@doAddChannel},...
        'visible','off');
    hok = uicontrol('style','pushbutton','position',[posfig(3)./2-150 10 300 25],...
        'string','Apply','userdata',hf,'callback',{@doOK},...
        'visible','off');
    pointers.hch = hch;
    pointers.hlo = hlo;
    pointers.hsa = hsa;
    pointers.hok = hok;
end
set(hf,'position',posfig,'resize','off');
figure(hf);
% Build montage editing GUI
x0 = [30 50+p*se se se];
x = x0;
xt = x0;
xt(2) = x0(2) + 20;
xt2 = x0;
xt2(1) = x0(1) - 20;
xt3 = [posfig(3)-100 50+p*se 90 se];
for i = 1:p
    % create editable boxes and add sensor names
    x(2) = x0(2) + (1-i).*se;
    for j = 1:n
        x(1) = x0(1) + (j-1).*se;
        if hce(i,j) == 0
           
            ts = ['sensor ', montage.labelorg{i},' --> channel ',num2str(i)];
            hce(i,j) = uicontrol('style','edit','position',x,...
                'tooltipstring',ts,'backgroundcolor',[0.95 0.75 0.75]);
        else
            set(hce(i,j),'position',x);
        end
        % add sensor name
        if i == 1
            xt(1) = x0(1) + (j-1).*se;
            if hcc(j) == 0
                hcc(j) = uicontrol('style','text','position',xt,...
                    'string', montage.labelorg{j},...
                    'tooltipstring',['sensor ',num2str(j),': ', montage.labelorg{j}]);
            else
                set(hcc(j),'position',xt);
            end
        end
    end
    % add channel numbers and config
    xt2(2) = x0(2) + (1-i).*se;
    if hct(i) == 0
        hct(i) = uicontrol('style','text','position',xt2,...
            'string', montage.labelnew{i},...
            'tooltipstring',['channel ', montage.labelnew{i}]);
    else
        set(hct(i),'position',xt2);
    end
    xt3(2) = x0(2) + (1-i).*se;
    if hct2(i) == 0
        hct2(i) = uicontrol('style','text','position',xt3);
    else
        set(hct2(i),'position',xt3);
    end
end
pointers.hct = hct;
pointers.hce = hce;
pointers.hct2 = hct2;
pointers.hcc = hcc;

set(pointers.hch,'visible','on');
set(pointers.hlo,'visible','on');
set(pointers.hsa,'visible','on');
set(pointers.hok,'visible','on');

return


