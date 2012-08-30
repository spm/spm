function jsubs = cfg2jobsubs(item, subs)

% function jsubs = cfg2jobsubs(item, subs)
% Return the subscript into the job tree for a given subscript vector into
% the val part of the cfg tree. In a cfg_choice, this is a struct reference
% to a field with the name of the tag of the corresponding child node.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg2jobsubs.m 4867 2012-08-30 13:04:51Z volkmar $

rev = '$Rev: 4867 $'; %#ok

% Only de-reference subscripts into .val{1} field
if isequal(subs(1:2), substruct('.','val','{}',{1})) 
    % return possible indices into .val{1} field
    citem = subsref(item, subs(1:2));
    if numel(subs) > 2
        jsubs1 = cfg2jobsubs(citem, subs(3:end));
    else
        jsubs1 = [];
    end;
    jsubs = [substruct('.', gettag(citem)) jsubs1];
else
    cfg_message('matlabbatch:cfg2jobsubs:wrongsubs', 'Inappropriate subscript reference in item ''%s''.', item.tag);
    jsubs = struct('type',{},'subs',{});
end;