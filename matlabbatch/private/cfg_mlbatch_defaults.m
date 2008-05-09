function cfg_defaults = cfg_mlbatch_defaults

% function cfg_defaults = cfg_mlbatch_defaults
% This file contains defaults that control the behaviour and appearance 
% of matlabbatch.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: copyright_cfg.m 218 2008-04-17 10:34:11Z glauche $

rev = '$Rev: 218 $';

% Font definition for cfg_ui user interface
% cfg_defaults.cfg_ui.font is a font struct as returned by uisetfont
cfg_defaults.cfg_ui.font.FontName   = 'FixedWidth';
cfg_defaults.cfg_ui.font.FontWeight = 'normal';
cfg_defaults.cfg_ui.font.FontAngle  = 'normal';
cfg_defaults.cfg_ui.font.FontSize   = 10;
cfg_defaults.cfg_ui.font.FontUnits  = 'points';
% Toggle ExpertEdit mode. Value can be 'on' or 'off'
cfg_defaults.cfg_ui.ExpertEdit = 'off';