function hdr = empty_hdr
% Create an empty NIFTI-1 header
% FORMAT hdr = empty_hdr
% _______________________________________________________________________
% @(#)empty_hdr.m	1.1 John Ashburner 04/11/26

org = niftistruc;
for i=1:length(org),
    hdr.(org(i).label) = org(i).def;
end;

