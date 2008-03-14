function [sts val] = subsasgn_check(item,subs,val)

% function [sts val] = subsasgn_check(item,subs,val)
% Perform assignment checks for .num field. Checks for .val field could
% include filtering and numel checks, if inputs are not passed as
% reference.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: subsasgn_check.m 1215 2008-03-14 21:32:25Z volkmar $

rev = '$Rev: 1215 $';

sts = true;
switch subs(1).subs
    case {'num'}
	sts = subsasgn_check_num(val);
    case {'val'}
        if ~isa(val, 'cfg_dep')
            sts = iscell(val) && (isempty(val) || iscellstr(val{1}));
        end;
end;
