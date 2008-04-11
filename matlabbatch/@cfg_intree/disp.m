function disp(varargin)

% function disp(varargin)
% This class should not display any information about its structure.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: disp.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

disp(class(varargin{1}))