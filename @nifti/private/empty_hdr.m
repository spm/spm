function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

%
% $Id: empty_hdr.m 4959 2012-09-24 18:26:31Z guillaume $


org = niftistruc;
for i=1:length(org)
    hdr.(org(i).label) = feval(org(i).dtype.conv,org(i).def);
end
