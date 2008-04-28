function [fnames, defaults] = mysubs_fields

% function [fnames, defaults] = mysubs_fields
% This function returns a cell string of names containing the fields
% implemented by a derived class and their default values. It is called
% from the class constructor directly and indirectly for subsasgn/subsref
% via the subs_fields public function of each class.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: mysubs_fields.m 1489 2008-04-28 06:32:55Z volkmar $

rev = '$Rev: 1489 $';

fnames = {'tag','name','val','check','help','expanded','hidden','def'};
defaults = {'generic','Generic',{{}},'',{''},true,false,{{}}};