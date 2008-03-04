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
% $Id: display.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

disp(' ');
disp([inputname(1),' = ']);
disp(' ');
disp(item);
disp(' ');
