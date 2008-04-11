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
% $Id: mysubs_fields.m 1366 2008-04-11 10:24:17Z volkmar $

rev = '$Rev: 1366 $';

fnames = {'tag','name','val','check','help','id','expanded','hidden'};
defaults = {'generic','Generic',{},'',{''},[],true,false};