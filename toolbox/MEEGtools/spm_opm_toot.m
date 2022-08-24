function spm_opm_toot
% A helper function to audibly notify the end of a script, or the working
% day perhaps?
%__________________________________________________________________________
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% George O'Neill
% $Id: spm_opm_toot.m 8306 2022-08-24 15:41:21Z george $

str = 'bG9hZCB0cmFpbjsgc291bmQoeSxGcyk=';
eval(native2unicode(matlab.net.base64decode(str)));