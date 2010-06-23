function spec = cfg_findspec(cellspec)

% function spec = cfg_findspec(cellspec)
% Create a find specification. cellspec should contain a cell array of
% cells, each of them containing name/value pairs that will be combined
% into a struct suitable for e.g. @cfg_item/match and @cfg_item/list. 
% These methods will be used to e.g. select items in a configuration tree or
% to match dependencies and input items.
%
% Name/value pairs within a cell will be OR concatenated, while cells
% will be AND concatenated.
%
% A cellspec
%  {{'field1','val1','field2','val2'},{'field3','val3'}}
%
% matches an item if 
%  (item.field1==val1 || item.field2==val2) && item.field3==val3
%
% If the field name is 'class', an item matches, if its class name is equal to
% spec.value.
%
% For class specific matching rules, see the help for the
% resp. @cfg_.../match method.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_findspec.m 3944 2010-06-23 08:53:40Z volkmar $

rev = '$Rev: 3944 $'; %#ok

if nargin == 0 || isempty(cellspec)
    spec = {};
    return;
end;
spec = cell(size(cellspec));
for k = 1:numel(cellspec)
    spec{k} = struct('name',{}, 'value',{});
    for l = 1:2:numel(cellspec{k})
        spec{k}(end+1).name = cellspec{k}{l};
        spec{k}(end).value = cellspec{k}{l+1};
    end;
end;