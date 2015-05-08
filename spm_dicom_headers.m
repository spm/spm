function hdr = spm_dicom_headers(P, essentials)
% Read header information from DICOM files
% FORMAT hdr = spm_dicom_headers(P [,essentials])
% P          - array of filenames
% essentials - if true, then only save the essential parts of the header
%
% hdr        - cell array of headers, one element for each file.
%
% Contents of headers are approximately explained in:
% http://medical.nema.org/standard.html
%
% This code may not work for all cases of DICOM data, as DICOM is an
% extremely complicated "standard".
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_headers.m 6431 2015-05-08 18:24:28Z john $

if nargin<2, essentials = false; end

dict = readdict;
j    = 0;
hdr  = {};
if size(P,1)>1, spm_progress_bar('Init',size(P,1),'Reading DICOM headers','Files complete'); end
for i=1:size(P,1)
    tmp = spm_dicom_header(P(i,:),dict);
    if ~isempty(tmp)
        if isa(essentials,'function_handle')
            tmp = feval(essentials,tmp);
        elseif essentials
            tmp = spm_dicom_essentials(tmp);
        end
        if ~isempty(tmp)
            j      = j + 1;
            hdr{j} = tmp;
        end
    end
    if size(P,1)>1, spm_progress_bar('Set',i); end
end
if size(P,1)>1, spm_progress_bar('Clear'); end


%==========================================================================
% function dict = readdict(P)
%==========================================================================
function dict = readdict(P)
if nargin<1, P = 'spm_dicom_dict.mat'; end
try
    dict = load(P);
catch
    fprintf('\nUnable to load the file "%s".\n', P);
    rethrow(lasterror);
end

