function o = niftistruc(fmt)
% Create a data structure describing NIFTI headers
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: niftistruc.m 7147 2017-08-03 14:07:01Z spm $


if ~nargin, fmt = 'nifti1'; end
switch lower(fmt)
    case {'nifti1','ni1','n+1'}
        o = nifti1struc;
    case {'nifti2','ni2','n+2'}
        o = nifti2struc;
    otherwise
        error('Unknown format.');
end
