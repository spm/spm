function varargout = spm_file_ext(ext)
% Return or set file extension for SPM images
% FORMAT ext = spm_file_ext
% ext   - file extension (e.g. '.img' or '.nii' for NIfTI images)
%
% FORMAT spm_file_ext(ext)
% ext   - file extension (e.g. '.img' or '.nii' for NIfTI images)
%__________________________________________________________________________
%
% The file extension returned by this function is defined in spm_defaults.m
% in field 'defaults.images.format'.
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_file_ext.m 4422 2011-08-04 16:18:05Z guillaume $ 

if ~nargin
    defaults = spm_get_defaults;
    % faster access than spm_get_defaults('images.format')
    varargout = {['.' defaults.images.format]};
else
    if ~isempty(ext) && ext(1)=='.', ext = ext(2:end); end
    spm_get_defaults('images.format', ext);
end
