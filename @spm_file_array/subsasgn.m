function a=subsasgn(a,b,c)
% Overloaded subsasgn function for spm_file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: subsasgn.m 253 2005-10-13 15:31:34Z guillaume $

error('Sorry, but the elements of spm_file_array objects are read-only.');
