function jsubs = cfg2jobsubs(item, subs)

% function jsubs = cfg2jobsubs(item, subs)
% Return the subscript into the job tree for a given subscript vector into
% the val part of the cfg tree. This generic function should be called only
% for leafs (cfg_entry, cfg_file, cfg_menu) of the cfg tree. It returns the
% subscripts that remain after item.val{1} has been addressed.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg2jobsubs.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

% Only de-reference subscripts into .val{1} field
if isequal(subs(1:2), substruct('.','val','{}',{1}))
    % return possible indices into .val{1} field
    jsubs = subs(3:end);
else
    warning('matlabbatch:cfg2jobsubs:wrongsubs', 'Inappropriate subscript reference in item ''%s''.', item.tag);
    jsubs = struct('type',{},'subs',{});
end;