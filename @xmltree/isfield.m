function F = isfield(tree,uid,parameter)
% XMLTREE/ISFIELD Is parameter a field of tree{uid} ?
% FORMAT F = isfield(tree,uid,parameter)
%
% tree      - a tree
% uid       - uid of the element
% parameter - a field of the root tree
% F         - 1 if present, 0 otherwise
%__________________________________________________________________________
%
% Is parameter a field of tree{uid} ?
%__________________________________________________________________________

% Copyright (C) 2002-2022 Guillaume Flandin


F = false(1,length(uid));
for i=1:length(uid)
    if isfield(tree.tree{uid(i)},parameter)
        F(i) = true;
    end
end
