function spm_pulldown
% Create a pulldown for individual jobs
%_______________________________________________________________________
% %W% John Ashburner %E%

c  = spm_config;
fg = spm_figure('findwin','Graphics');
if isempty(fg), return; end;
set(0,'ShowHiddenHandles','on');
delete(findobj(fg,'tag','jobs'));
set(0,'ShowHiddenHandles','off');
f0 = uimenu(fg,'Label','TASKS','HandleVisibility','off','tag','jobs');
rec(f0,c,c.tag);
f1 = uimenu(f0,'Label','Batch','CallBack','spm_jobman;','Separator','on');
f1 = uimenu(f0,'Label','Defaults','CallBack','spm_jobman(''defaults'');','Separator','off');
return;

function rec(f0,c0,tag0)
for i=1:length(c0.values),
    c1 = c0.values{i};
    if isfield(c1,'tag'),
        tag1 = [tag0 '.' c1.tag];
    end;
    if isfield(c1,'prog'),
        str = ['spm_jobman(''interactive'','''',''' tag1 ''');'];
        f1 = uimenu(f0,'Label',c1.name,'CallBack',str);
    else
        f1 = uimenu(f0,'Label',c1.name);
        rec(f1,c1,tag1);
    end;
end;
return;
