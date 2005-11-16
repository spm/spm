function varargout = spm_list_files(varargin)
% Compiled routine that lists files and directories
% FORMAT [Files,Dirs] = spm_list_files(Dir,Filter)
% Dir    - directory to list
% Filter - e.g. '*.img'
% Files  - filenames
% Dirs   - directories
%_______________________________________________________________________
%
% See also: spm_get.m
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_list_files.m 300 2005-11-16 21:05:24Z guillaume $


warning('spm_list_files is a bit old now. Please try to use spm_select.');

error(nargchk(2,2,nargin));

filt = varargin{2};
filt = strrep(filt, '.', '\.'  );
filt = strrep(filt, '?', '.{1}');
filt = strrep(filt, '*', '.*'  );

[varargout{1},varargout{2}] = spm_select('List',varargin{1},filt);
