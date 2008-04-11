function obj = cfg_inv_out(varargin)
    
% function obj = cfg_inv_out(varargin)
% Auxiliary class to mark invalid (i.e. not yet available) outputs of
% cfg_exbranch'es. An object of this type will be assigned automatically
% to a cfg_exbranch'es .jout field. resolve_deps will not resolve outputs
% that consist of a cfg_inv_out object.
%
% This code is part of a batch job configuration system for MATLAB. See 
%      help matlabbatch
% for a general overview.
%_______________________________________________________________________
% Copyright (C) 2007 Freiburg Brain Imaging

% Volkmar Glauche
% $Id: cfg_inv_out.m 1373 2008-04-11 14:24:03Z spm $

rev = '$Rev: 1373 $';

% Dummy field - nothing should be stored here
obj.invalid = true;
obj = class(obj, 'cfg_inv_out');