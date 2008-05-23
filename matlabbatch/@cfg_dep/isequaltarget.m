function sts = isequaltarget(dep1, dep2)

% function sts = isequaltarget(dep1, dep2)
% Compare source references of two dependencies and return true if both
% point to the same object. If multiple dependencies are given, the
% number and order of dependencies must match.
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: isequaltarget.m 1716 2008-05-23 08:18:45Z volkmar $

rev = '$Rev: 1716 $'; %#ok
sts = numel(dep1) == numel(dep2);
if sts
    for k = 1:numel(dep1)
        sts = sts && isequal(dep1.tgt_exbranch, dep2.tgt_exbranch) && ...
              isequal(dep1.tgt_input, dep2.tgt_input);
        if ~sts
            break;
        end;
    end;
end;
