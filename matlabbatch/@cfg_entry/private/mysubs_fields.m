function [fnames defaults] = mysubs_fields

% Additional fields for class cfg_entry. See help of
% @cfg_item/subs_fields for general help about this function.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: mysubs_fields.m 1246 2008-03-26 10:45:13Z volkmar $

rev = '$Rev: 1246 $';

fnames = {'strtype','num','def','extras'};
defaults = {'e',[],{},{}};