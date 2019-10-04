function spm_eeg_display_tf(files)
% Display of TF images saved as NIfTI
% Up to 6 images are supported 
% FORMAT spm_eeg_display_tf
% FORMAT spm_eeg_display_tf(files)
%        files - list of images to display 
%                (as char or cell aray of strings)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% Vladimir Litvak
% $Id: 

if nargin == 0
    files = cellstr(spm_select([1 6], 'image', 'Select TF image files'));%
elseif isa(files, 'char')
    files = cellstr(files);
end
    
%%
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph);

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
    n = nifti(files{f});
    
    if ~ismatrix(n.dat) || ~all(size(n.dat)>1)
        error('2D image expected');
    end
    
    f = (1:size(n.dat, 1))*n.mat(1,1)+n.mat(1, 4);
    t = (1:size(n.dat, 2))*n.mat(2,2)+n.mat(2, 4);
    
    imagesc(t, f, n.dat(:,:));
    axis xy;
    xlabel('time (ms)');
    ylabel('frequency (Hz)');
    set(gca, 'FontSize', 15);
end