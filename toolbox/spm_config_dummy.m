function ob = spm_config_dummy(varargin)
% Dummy configuration file
%
% This file should remain in the toolbox directory, so that other
% toolboxes can be more easily added by including spm_config_*.m files.
% It ensures that the data structure is of the form jobs{*}.tools{*}.*
%
% If you write a toolbox, then you can include it in this directory -
% but remember to try to keep the function names unique (to reduce
% clashes with other toolboxes.  See spm_config.m for information about
% the form of SPM's configuration files.
% 
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_config_dummy.m 245 2005-09-27 14:16:41Z guillaume $

ob = struct('type','branch','tag','unused','name','Toolboxes','val',{{}});
ob.help = {[...
'Toolbox configuration files should be placed in the toolbox directory, ',...
'with their own spm_config_*.m files. ',...
'If you write a toolbox, then you can include it in this directory - ',...
'but remember to try to keep the function names unique (to reduce ',...
'clashes with other toolboxes.  See spm_config.m for information about ',...
'the form of SPM''s configuration files.']};
