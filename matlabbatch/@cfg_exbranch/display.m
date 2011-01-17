function display(item)

% function display(item)
% Display a configuration object
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: display.m 4166 2011-01-17 15:06:41Z volkmar $

rev = '$Rev: 4166 $'; %#ok

disp(' ');
disp([inputname(1),' = ']);
disp(' ');
disp(item);
disp(' ');
