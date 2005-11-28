function scales = spm_eeg_write(fpout, d, direc, dtype)
% write EEG data in a scaled format to file
% FORMAT scales = spm_eeg_write(fpout, d, direc, dtype)
%
% fpout  - file descriptor to write to
% d      - data matrix
% direc  - Dimension along which the data is scaled
% dtype  - data type with which data is written
% 
% scales - scaling vector
%_______________________________________________________________________
% 
% spm_eeg_write is used to write data to different data types. If the data
% type is different from float, the data is compressed using the vector
% scales. The only reason for writing to a format like int16 is to save
% disk space.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_write.m 317 2005-11-28 18:31:24Z stefan $


if strcmp(dtype, 'float32')
    scales = ones(size(d, 3-direc), 1);
elseif strcmp(dtype, 'int16')
    scales = max(abs(d), [], direc)./32767;
    scales = scales(:);
    ind = find(scales == 0);
    scales(ind) = ones(size(ind));
    d = int16(round(d./repmat(scales, 1, size(d, 2))));
elseif strcmp(dtype, 'int32')
    scales = max(abs(d), [], direc)./2147483647;
    scales = scales(:);
    ind = find(scales == 0);
    scales(ind) = ones(size(ind));
    d = int32(d./repmat(scales, 1, size(d, 2)));
end

if fpout ~= - 1
    fwrite(fpout, d, dtype);
end