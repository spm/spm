function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id$


org = niftistruc;
for i=1:length(org),
    hdr.(org(i).label) = org(i).def;
end;

