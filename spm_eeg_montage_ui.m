function montage = spm_eeg_montage_ui(montage)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_eeg_montage_ui.m 2258 2008-09-30 16:50:29Z jean $

%% creates GUI from spm_uitable
fig = figure;
pos = get(fig,'position');
pos2 = [40 70 pos(3)-60 pos(4) - 100];
pos = [pos(1) pos(2) 1.8*pos(3) pos(4)];
set(gcf,'menubar','none','position',pos,...
    'numbertitle','off','name','Montage edition');

ha = gca;
set(ha,'position',pos2);

addButtons;

table = cat(2,montage.labelnew(:),num2cell(montage.tra));
colnames = cat(2,'channel labels',montage.labelorg(:)');

[ht,hc] = spm_uitable(table,colnames);
set(ht,'position',pos2,...
    'units','normalized');

ud.hi = imagesc(montage.tra);
set(ha,'position',[0.6 0.18 0.4 0.77])
axis image
colormap('bone')
zoom

ud.b4MontageEditing = get(0,'userdata');
ud.ht = ht;
ud.fig = fig;
ud.montage = montage;
set(0,'userdata',ud);



%% gets the montage from the GUI
uiwait(fig)

ud = get(0,'userdata');
if isfield(ud,'b4MontageEditing')
    montage = ud.montage;
    set(0,'userdata',ud.b4MontageEditing);
else % GUI was closed without 'OK' button
    montage = [];
end




%% 'add row' button subfunction
function [ha] = doAddRow(o1,o2)
ud = get(0,'userdata');
[M,newLabels] = getM(ud.ht);
M = [M;zeros(1,size(M,2))];
newLabels = cat(1,newLabels(:),num2str(size(M,1)));
set(ud.ht,'units','pixels');
pos = get(ud.ht,'Position');
delete(ud.ht);
table = cat(2,newLabels,num2cell(M));
colnames = cat(2,'channel labels',ud.montage.labelorg(:)');
pause(1) % This is weird, but fixes java troubles.
ht = spm_uitable(table,colnames);
set(ht,'position',pos,...
    'units','normalized');
ud.ht = ht;
set(0,'userdata',ud);
doCheck

%% 'load' button subfunction
function [ha] = doLoad(o1,o2)
[FileName,PathName] = uigetfile('*.mat','Load montage file');
if ~isempty(FileName) || ~isequal(FileName,0)
    montage = load(fullfile(PathName,FileName));
    name = fieldnames(montage);
    if isequal(name{1},'montage')
        montage = getfield(montage, name{1});
        ud = get(0,'userdata');
        set(ud.ht,'units','pixels');
        pos = get(ud.ht,'Position');
        delete(ud.ht);
        table = cat(2,montage.labelnew(:),num2cell(montage.tra));
        colnames = cat(2,'channel labels',montage.labelorg(:)');
        pause(1) % This is weird, but fixes java troubles.
        ht = spm_uitable(table,colnames);
        set(ht,'position',pos,...
            'units','normalized');
        ud.ht = ht;
        ud.montage = montage;
        set(0,'userdata',ud);
        pause(1)
        doCheck
    else
        warndlg('File did not contain any montage!')
    end
end


%% 'save as' button subfunction
function [ha] = doSave(o1,o2)
ud = get(0,'userdata');
[M,newLabels] = getM(ud.ht);
montage.tra = M;
montage.labelorg = ud.montage.labelorg;
montage.labelnew = newLabels;
uisave('montage','SPMeeg_montage.mat');


%% 'OK' button subfunction
function [] = doOK(o1,o2)
ud = get(0,'userdata');
[M,newLabels] = getM(ud.ht);
% delete row if empty:
ind = find(sum(M,2)==0);
M(ind,:) = [];
newLabels = {newLabels{setdiff(1:length(newLabels),ind)}};
montage.tra = M;
montage.labelorg = ud.montage.labelorg(:);
montage.labelnew = newLabels(:);
ud.montage = montage;
set(0,'userdata',ud);
pause(1)
close(gcf)


function [] = doCheck(o1,o2)
ud = get(0,'userdata');
[M,newLabels] = getM(ud.ht);
set(ud.hi,'cdata',M)
set(gca,'xlim',[0.5 size(M,1)])
set(gca,'ylim',[0.5 size(M,2)])
axis image


%% extracting montage from java object
function [M,newLabels] = getM(ht)
nnew = get(ht,'NumRows');
nold = get(ht,'NumColumns')-1;
M = zeros(nnew,nold);
data = get(ht,'data');
for i =1:nnew
    if ~isempty(data(i,1))
        newLabels{i} = data(i,1);
    else
        newLabels{i} = [];
    end
    for j =1:nold
        if ~isempty(data(i,j+1))
            if ~isstr(data(i,j+1))
                M(i,j) = data(i,j+1);
            else
                M0 = str2num(data(i,j+1));
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


%% adding buttons to the montage GUI
function [] = addButtons(ha)
hAdd = uicontrol('style','pushbutton',...
    'string','Add row','callback',{@doAddRow},...
    'position',[60 20 80 20]);
set(hAdd,'units','normalized')
hLoad = uicontrol('style','pushbutton',...
    'string','load file','callback',{@doLoad},...
    'position',[180 20 80 20]);
set(hLoad,'units','normalized')
hSave = uicontrol('style','pushbutton',...
    'string','save as','callback',{@doSave},...
    'position',[280 20 80 20]);
set(hSave,'units','normalized')
hOK = uicontrol('style','pushbutton',...
    'string',' OK ','callback',{@doOK},...
    'position',[400 20 80 20]);
set(hOK,'units','normalized')
hCheck = uicontrol('style','pushbutton',...
    'string',' Check montage ','callback',{@doCheck},...
    'position',[760 20 120 20]);
set(hCheck,'units','normalized')

