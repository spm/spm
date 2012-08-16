function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save.m 4847 2012-08-16 17:29:23Z vladimir $

if nargin == 1 && exist(fullfile(pwd, 'BF.mat'), 'file')
    save('BF.mat', '-struct', 'BF', '-append');
else
    save('BF.mat', '-struct', 'BF');
end