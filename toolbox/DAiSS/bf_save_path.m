function bf_save_path(BF,path, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save_path.m 7703 2019-11-22 12:06:29Z guillaume $

if nargin == 2 && exist(path, 'file')
    save(path, '-struct', 'BF', '-append');
else
    save(path, '-struct', 'BF', '-v7.3');
end