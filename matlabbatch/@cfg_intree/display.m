function display(varargin)

% This class should not display any information about its structure.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: display.m 1184 2008-03-04 16:27:57Z volkmar $

rev = '$Rev: 1184 $';

disp(class(varargin{1}))