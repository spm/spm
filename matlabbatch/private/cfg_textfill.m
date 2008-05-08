function str = cfg_textfill(obj, left, right, tflag)

% function str = cfg_textfill(obj, left, right)
% Fill a text object, so that the left part is left justified and the
% right part right justified. If tflag is set, try to fit text in widget
% by truncating right until at least 5 characters are displayed.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: copyright_cfg.m 218 2008-04-17 10:34:11Z glauche $

rev = '$Rev: 218 $';

if ~ishandle(obj)
    error('matlabbatch:textfill:argchk',...
          'First input must be a graphics handle.');
else
    try
        get(obj,'string');
    catch
        error('matlabbatch:textfill:argchk',...
              'Input object must have a ''string'' property.');
    end;
end;
if ~iscellstr(left)
    if ischar(left)
        left = cellstr(left);
    else
        error('matlabbatch:textfill:argchk',...
              'Second input must be a string array or cellstr.');
    end;
    % add one space as delimiter
    left = strcat(left, ' ');
end;
if ~iscellstr(right)
    if ischar(right)
        right = cellstr(right);
    else
        error('matlabbatch:textfill:argchk',...
              'Third input must be a string array or cellstr.');
    end;
end;
if numel(left) ~= numel(right)
    error('matlabbatch:textfill:argchk',...
          'Second and third input must have the same number of lines.');
end;    

TempObj=copyobj(obj,get(obj,'Parent'));
set(TempObj,'Visible','off','Max',100);

% Find max extent of left string
lext = maxextent(TempObj, left);
mlext = max(lext);

% Find max extent of right string
rext = maxextent(TempObj, right);
mrext = max(rext);

% Find extent of single space
% Work around MATLAB inaccuracy by measuring 100 spaces and dividing by 100
spext = maxextent(TempObj, {repmat(' ',1,100)})/100;

% try to work out slider size
pos = get(TempObj,'Position');
oldun = get(TempObj,'units');
set(TempObj,'units','points');
ppos = get(TempObj,'Position');
set(TempObj,'units',oldun);
sc = pos(3)/ppos(3);
% assume slider width of 15 points
swidth=15*sc;

pos = get(obj,'Position');
width = max(mlext+mrext,pos(3)-swidth);
if tflag && mlext+mrext > pos(3)-swidth
    newrext = max(5*spext, pos(3)-swidth-mlext);
    trind = find(rext > newrext);
    for k = 1:numel(trind)
        tright = right{trind(k)};
        set(TempObj,'String',['...' tright]);
        ext = get(TempObj,'Extent');
        while ext(3) > newrext
            tright = tright(2:end);
            set(TempObj,'String',['...' tright]);
            ext = get(TempObj,'Extent');
        end;
        right{trind(k)} = ['...' tright];
        rext(trind(k)) = ext(3);
    end;
    width = mlext+newrext;
end;

str = cell(size(left));
for k = 1:numel(left)
    fillstr = repmat(' ',1,floor((width-(lext(k)+rext(k)))/spext));
    str{k} = [left{k} fillstr right{k}];
end;

function ext = maxextent(obj, str)
ext = zeros(size(str));
for k = 1:numel(str)
    set(obj,'String',str{k});
    next = get(obj, 'Extent');
    ext(k) = next(3);
end;
