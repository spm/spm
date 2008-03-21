function obj = putlabel(obj, ind, label)
% Method for putting new labels for all trials ind
% FORMAT res = putlabel(obj, ind, label)
%
% ind: vector of trial indices
% label: either cell-vector of strings or character array
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

if iscell(label)
    [obj.trials(ind).label] = deal(label);
else
    [obj.trials(ind).label] = deal(cellstr(label));
end