function sts = isequalsource(dep1, dep2)

% function sts = isequalsource(dep1, dep2)
% Compare source references of two dependencies and return true if both
% point to the same object. If multiple dependencies are given, the
% number and order of dependencies must match.
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: isequalsource.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';
sts = numel(dep1) == numel(dep2);
if sts
    for k = 1:numel(dep1)
        sts = sts && isequal(dep1(k).src_exbranch, dep2(k).src_exbranch) && ...
              isequal(dep1(k).src_output, dep2(k).src_output);
        if ~sts
            break;
        end;
    end;
end;
