function spec = cfg_findspec(cellspec)

% function spec = cfg_findspec(cellspec)
% Convert a shorthand for a find spec into a struct/cell array. cellspec
% should contain a cell array of cells, each of them containing name/value
% pairs that will be combined into a struct suitable for match/find.
% Name/value pairs within a cell will be OR concatenated, while cells
% will be AND concatenated.
% A cellspec
%  {{'field1','val1','field2','val2'},{'field3','val3'}}
% matches an item if 
%  (item.field1==val1 || item.field2==val2) && item.field3==val3
% If the field name is 'class', an item matches, if its class name is equal to
% spec.value.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_findspec.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok

if nargin == 0 || isempty(cellspec)
    spec = {};
    return;
end;
for k = 1:numel(cellspec)
    spec{k} = struct('name',{}, 'value',{});
    for l = 1:2:numel(cellspec{k})
        spec{k}(end+1).name = cellspec{k}{l};
        spec{k}(end).value = cellspec{k}{l+1};
    end;
end;