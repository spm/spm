function montage = spm_eeg_montage_ui(montage)
% GUI for EEG montage (rereference EEG data to new reference channel(s))
% FORMAT montage = spm_eeg_montage_ui(montage)
%
% montage     - structure with fields:
%   tra       - MxN matrix
%   labelnew  - Mx1 cell-array - new labels
%   labelorg  - Nx1 cell-array - original labels
%
% Output is empty if the GUI is closed.
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_montage_ui.m 7119 2017-06-20 11:25:48Z guillaume $


% Create the figure
%--------------------------------------------------------------------------
fig  = figure;
S0   = spm('WinSize','0',1);
pos  = get(fig,'position');
pos2 = [40 70 pos(3)-60 pos(4)-100];
pos  = [S0(1) S0(2) 0 0] + [pos(1) pos(2) 1.8*pos(3) pos(4)];
set(fig,...
    'menubar',     'none',...
    'position',    pos,...
    'numbertitle', 'off',...
    'name',        'Montage edition');

addButtons(fig);

% Display the uitable component
%--------------------------------------------------------------------------
table    = cat(2,montage.labelnew(:),num2cell(montage.tra));
colnames = cat(2,'channel labels',montage.labelorg(:)');

pause(1e-1) % This is weird, but fixes java troubles.
ht       = my_uitable(table,colnames);
set(ht,'position',pos2, 'units','normalized');

% Display the matrix representation of the montage 
%--------------------------------------------------------------------------
ax = axes('position',[0.6 0.18 0.4 0.77]);
hi = imagesc(montage.tra,'parent',ax);
axis('image');
colormap('bone');
zoom(fig,'on');

% Store info in figure's userdata and wait for user interaction
%--------------------------------------------------------------------------
ud.hi      = hi;
ud.ht      = ht;
ud.montage = montage;
set(fig,'userdata',ud);
uiwait(fig);

% Get the montage from the GUI
%--------------------------------------------------------------------------
try
    ud = get(fig,'userdata');
    montage = ud.montage;
    close(fig);
catch % GUI was closed without 'OK' button
    montage = [];
end


%==========================================================================
function doAddRow(obj,evd,h)
%==========================================================================
% 'add row' button subfunction
ud = get(h,'userdata');
[M,newLabels] = getM(ud.ht);
M = [M;zeros(1,size(M,2))];
newLabels = cat(1,newLabels(:),num2str(size(M,1)));
set(ud.ht,'units','pixels');
pos = get(ud.ht,'Position');
delete(ud.ht);
table = cat(2,newLabels,num2cell(M));
colnames = cat(2,'channel labels',ud.montage.labelorg(:)');
pause(1) % This is weird, but fixes java troubles.
ht = my_uitable(table,colnames);
set(ht,'position',pos, 'units','normalized');
ud.ht = ht;
set(h,'userdata',ud);
doCheck(obj,evd,h);

%==========================================================================
function doLoad(obj,evd,h)
%==========================================================================
% 'load' button subfunction
[t,sts] = spm_select(1,'mat','Load montage file');
if sts
    montage = load(t);
    if ismember('montage', fieldnames(montage))
        montage    = montage.montage;
        ud         = get(h,'userdata');
        set(ud.ht,'units','pixels');
        pos        = get(ud.ht,'Position');
        delete(ud.ht);
        table      = cat(2,montage.labelnew(:),num2cell(montage.tra));
        colnames   = cat(2,'channel labels',montage.labelorg(:)');
        pause(1) % This is weird, but fixes java troubles.
        ht         = my_uitable(table,colnames);
        set(ht,'position',pos,...
            'units','normalized');
        ud.ht      = ht;
        ud.montage = montage;
        set(h,'userdata',ud);
        pause(1)
        doCheck(obj,evd,h);
    else
        spm('alert!','File did not contain any montage!','Montage edition');
    end
end

%==========================================================================
function doSave(obj,evd,h)
%==========================================================================
% 'save as' button subfunction
doCheck(obj,evd,h);
ud = get(h,'userdata');
[M,newLabels] = getM(ud.ht);
% delete row if empty:
ind              = ~any(M,2);
M(ind,:)         = [];
newLabels(ind)   = [];
montage.tra      = M;
montage.labelorg = ud.montage.labelorg;
montage.labelnew = newLabels;
[filename, pathname] = uiputfile({'*.mat','MAT-files (*.mat)'}, ...
    'Save montage', 'SPMeeg_montage.mat');
if ~isequal(filename, 0)
	save(fullfile(pathname, filename), 'montage', spm_get_defaults('mat.format'));
end

%==========================================================================
function doOK(obj,evd,h)
%==========================================================================
% 'OK' button subfunction
doCheck(obj,evd,h);
ud               = get(h,'userdata');
[M, newLabels]   = getM(ud.ht);
% delete row if empty:
ind              = ~any(M,2);
M(ind,:)         = [];
newLabels(ind)   = [];
montage.tra      = M;
montage.labelorg = ud.montage.labelorg(:);
montage.labelnew = newLabels(:);
ud.montage = montage;
set(h,'userdata',ud);
uiresume(h);

%==========================================================================
function doCheck(obj,evd,h)
%==========================================================================
% Update the montage display
ud = get(h,'userdata');
M  = getM(ud.ht);
set(ud.hi,'cdata',M);
set(gca,'xlim',[0.5 size(M,1)]);
set(gca,'ylim',[0.5 size(M,2)]);
axis('image');
drawnow;

%==========================================================================
function [M,newLabels] = getM(ht)
%==========================================================================
% extracting montage from java object
nnew = get(ht,'NumRows');
nold = get(ht,'NumColumns')-1;
M    = zeros(nnew,nold);
data = get(ht,'data');
for i =1:nnew
    if ~isempty(data(i,1))
        newLabels{i} = data(i,1);
    else
        newLabels{i} = [];
    end
    for j =1:nold
        if ~isempty(data(i,j+1))
            if ~ischar(data(i,j+1))
                M(i,j) = data(i,j+1);
            else
                M0 = str2double(data(i,j+1));
                if ~isempty(M0)
                    M(i,j) = M0;
                else
                    M(i,j) = 0;
                end
            end
        else
            M(i,j) = 0;
        end
    end
end

%==========================================================================
function addButtons(h)
%==========================================================================
% adding buttons to the montage GUI
hAdd = uicontrol('style','pushbutton',...
    'string','Add row','callback',{@doAddRow,h},...
    'position',[60 20 80 20]);
set(hAdd,'units','normalized');
hLoad = uicontrol('style','pushbutton',...
    'string','Load file','callback',{@doLoad,h},...
    'position',[180 20 80 20]);
set(hLoad,'units','normalized');
hSave = uicontrol('style','pushbutton',...
    'string','Save as','callback',{@doSave,h},...
    'position',[280 20 80 20]);
set(hSave,'units','normalized');
hOK = uicontrol('style','pushbutton',...
    'string',' OK ','callback',{@doOK,h},...
    'position',[400 20 80 20]);
set(hOK,'units','normalized');
hCheck = uicontrol('style','pushbutton',...
    'string',' Refresh display ','callback',{@doCheck,h},...
    'position',[760 20 120 20]);
set(hCheck,'units','normalized');

%==========================================================================
function [ht,hc] = my_uitable(varargin)
%==========================================================================
% conversion layer for various MATLAB versions
persistent runOnce
try
    if spm_check_version('matlab','8.4') >= 0
        if isempty(runOnce)
            warning('Consider migrating to the new uitable component.');
            runOnce = true;
        end
        [ht,hc] = uitable('v0',varargin{:});
    else
        [ht,hc] = spm_uitable(varargin{:});
    end
catch
    [ht,hc]     = deal([]);
end
