function dep = dep_add(cdep, dep, ntgt_input, njtsubs)

% augment cdep tsubs references, and add them to dependency list
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: dep_add.m 2512 2008-12-01 13:21:29Z volkmar $

rev = '$Rev: 2512 $'; %#ok

for k = 1:numel(cdep)
    cdep(k).tgt_input = [ntgt_input cdep(k).tgt_input];
    cdep(k).jtsubs = [njtsubs cdep(k).jtsubs];
end;
if isempty(dep)
    dep = cdep(:);
else
    dep = [dep(:); cdep(:)];
end;
