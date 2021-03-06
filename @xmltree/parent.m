function p = parent(tree,uid)
% XMLTREE/PARENT Parent Method
% FORMAT uid = parent(tree,uid)
% 
% tree   - XMLTree object
% uid    - UID of the lonely child
% p      - UID of the parent ([] if root is the child)
%__________________________________________________________________________
%
% Return the uid of the parent of a node.
%__________________________________________________________________________

% Copyright (C) 2002-2022 Guillaume Flandin


p = tree.tree{uid}.parent;
