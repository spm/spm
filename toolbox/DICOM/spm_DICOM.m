function spm_DICOM
% Convert DICOM files into images for SPM2
%_______________________________________________________________________
% @(#)spm_DICOM.m	1.1 John Ashburner 02/08/12

P   = spm_get(Inf,'*','Select DICOM files');
hdr = spm_dicom_headers(P);
spm_dicom_convert(hdr)

