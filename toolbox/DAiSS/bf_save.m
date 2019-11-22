function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save.m 7703 2019-11-22 12:06:29Z guillaume $

if nargin == 1 && exist(fullfile(pwd, 'BF.mat'), 'file')
    save('BF.mat', '-struct', 'BF', '-append');
else
    save('BF.mat', '-struct', 'BF', '-v7');
end