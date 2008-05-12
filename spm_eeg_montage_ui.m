function montage = spm_eeg_montage_ui(montage)

dbstop if error

fig = figure;
pos = [520   627   559   494];
set(gcf,'menubar','none',...
    'numbertitle','off','name','Montage edition');
set(gcf,'position',pos);
pos2 = [40 70 pos(3)-60 pos(4) - 100];
ha = gca;
set(ha,'position',pos2);


addButtons;

ht = spm_uitable(num2cell(montage.tra),montage.labelorg);
set(ht,'position',pos2,...
    'units','normalized');

ud.ob4MontageEditing = get(0,'userdata');
ud.ht = ht;
ud.montage = montage;
set(0,'userdata',ud);


uiwait(fig)


ud = get(0,'userdata');
if isfield(ud,'b4MontageEdition')
    montage = ud.montage;
    for i=1:size(montage.tra,1)
        if i <= size(montage.tra,2)
            montage.labelnew{i} = montage.labelorg{i}
        else
            montage.labelnew{i} = num2str(i);
        end
    end
    set(0,'userdata',ud.b4MontageEdition);
else
    montage = [];
end






function [ha] = doAddRow(o1,o2)
dbstop if error
ud = get(0,'userdata');
M = getM(ud.ht);
M = [M;zeros(1,size(M,2))];
set(ud.ht,'units','pixels');
pos = get(ud.ht,'Position');
delete(ud.ht);
ht = spm_uitable(num2cell(M),ud.montage.labelorg);
set(ht,'position',pos,...
    'units','normalized');
ud.ht = ht;
set(0,'userdata',ud);



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
        ht = spm_uitable(num2cell(montage.tra),montage.labelorg);
        set(ht,'position',pos,...
            'units','normalized');
        ud.ht = ht;
        set(0,'userdata',ud);
    end
end

function [ha] = doSave(o1,o2)
dbstop if error
ud = get(0,'userdata');
M = getM(ud.ht);
montage.tra = M;
montage.labelorg = ud.montage.labelorg;
for i=1:size(montage.tra,1)
    if i <= size(montage.tra,2)
        montage.labelnew{i} = montage.labelorg{i};
    else
        montage.labelnew{i} = num2str(i);
    end
end
uisave('montage','SPMeeg_montage.mat');

function [] = doOK(o1,o2)
dbstop if error
ud = get(0,'userdata');
M = getM(ud.ht);
montage.tra = M;
montage.labelorg = ud.montage.labelorg;
ud.b4MontageEdition = get(0,'userdata');
ud.montage = montage;
set(0,'userdata',ud);
close(gcf)


function M = getM(ht)
M = zeros(get(ht,'NumRows'),get(ht,'NumColumns'));
data = get(ht,'data');
for i =1:get(ht,'NumRows')
    for j =1:get(ht,'NumColumns')
        M(i,j) = data(i,j);
    end
end



function [] = addButtons(ha)
dbstop if error
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

