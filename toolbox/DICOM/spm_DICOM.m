function spm_DICOM
% Convert DICOM files into images for SPM2
%_______________________________________________________________________
% %W% John Ashburner %E%

P   = spm_get(Inf,'*','Select DICOM files');
hdr = spm_dicom_headers(P);
spm_dicom_convert(hdr)

