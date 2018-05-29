function Headers = spm_dicom_headers(DicomFilenames, Essentials)
% Read header information from DICOM files
% FORMAT Headers = spm_dicom_headers(DicomFilenames [,Essentials])
% DicomFilenames - array of filenames
% Essentials     - if true, then only save the essential parts of the header
%
% Headers        - cell array of headers, one element for each file.
%
% Contents of headers are approximately explained in:
% http://medical.nema.org/standard.html
%
% This code may not work for all cases of DICOM data, as DICOM is an
% extremely complicated "standard".
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_headers.m 7320 2018-05-29 10:19:49Z john $

if nargin<2, Essentials = false; end

DicomDictionary = ReadDicomDictionary;
j        = 0;
Headers  = {};
if size(DicomFilenames, 1)>1, spm_progress_bar('Init',size(DicomFilenames, 1), 'Reading DICOM headers', 'Files complete'); end
for i=1:size(DicomFilenames,1)
    Header = spm_dicom_header(DicomFilenames(i,:), DicomDictionary);
    if ~isempty(Header)
        if isa(Essentials,'function_handle')
            Header = feval(Essentials, Header);
        elseif Essentials
            Header = spm_dicom_essentials(Header);
        end
        if ~isempty(Header)
            j          = j + 1;
            Headers{j} = Header;
        end
    end
    if size(DicomFilenames, 1)>1, spm_progress_bar('Set',i); end
end
if size(DicomFilenames, 1)>1, spm_progress_bar('Clear'); end


%==========================================================================
% function DicomDictionary = ReadDicomDictionary(DictionaryFile)
%==========================================================================
function DicomDictionary = ReadDicomDictionary(DictionaryFile)
if nargin<1, DictionaryFile = 'spm_dicom_dict.mat'; end
try
    DicomDictionary = load(DictionaryFile);
catch problem
    fprintf('\nUnable to load the file "%s".\n', DictionaryFile);
    rethrow(problem);
end

