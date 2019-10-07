function spm_eeg_display_tf(files)
% Display TF images saved as NIfTI
% FORMAT spm_eeg_display_tf(files)
% files  -  list of images to display (as char or cell aray of strings)
%           Up to 6 images are supported 
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_display_tf.m 7673 2019-10-07 09:07:14Z guillaume $


if ~nargin
    [files, sts] = spm_select([1 6], 'image', 'Select TF image files');
    if ~sts, return; end
end
files = cellstr(files);


Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);

switch numel(files)
    case 1
        dim = [1 1];
    case 2
        dim = [2 1];
    case 3
        dim = [3 1];
    case 4
        dim = [2 2];
    case {5, 6}
        dim = [3 2];
end

for f = 1:numel(files)
    subplot(dim(1), dim(2), f);
    N = nifti(files{f});
    
    if numel(N.dat.dim) ~= 2 || ~all(size(N.dat) > 1)
        error('2D image expected.');
    end
    
    fr = (1:size(N.dat, 1)) * N.mat(1,1) + N.mat(1, 4);
    t  = (1:size(N.dat, 2)) * N.mat(2,2) + N.mat(2, 4);
    
    imagesc(t, fr, N.dat(:,:));
    axis xy;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    set(gca, 'FontSize', 15);
end
