function o = niftistruc(fmt)
% Create a data structure describing NIFTI headers
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


if ~nargin, fmt = 'nifti1'; end
switch lower(fmt)
    case {'nifti1','ni1','n+1'}
        o = nifti1struc;
    case {'nifti2','ni2','n+2'}
        o = nifti2struc;
    otherwise
        error('Unknown format.');
end
