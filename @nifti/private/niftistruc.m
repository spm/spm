function o = niftistruc(fmt)
% Create a data structure describing NIFTI headers
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: niftistruc.m 4959 2012-09-24 18:26:31Z guillaume $


if ~nargin, fmt = 'nifti1'; end
switch lower(fmt)
    case 'nifti1'
        o = nifti1struc;
    case 'nifti2'
        o = nifti2struc;
    otherwise
        error('Unknown format.');
end