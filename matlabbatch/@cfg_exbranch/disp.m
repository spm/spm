function disp(obj)

% function disp(obj)
% Disp a cfg_exbranch object.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: disp.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

sz = size(obj);
fprintf('%s object: ', class(obj));
if length(sz)>4,
    fprintf('%d-D\n',length(sz));
else
    for i=1:(length(sz)-1),
        fprintf('%d-by-',sz(i));
    end;
    fprintf('%d\n',sz(end));
end;
if prod(sz)==1,
    so = struct(obj);
    % display branch fields
    disp(so.cfg_branch);
    % display exbranch additional fields
    disp(rmfield(so,'cfg_branch'));
end;
return;
