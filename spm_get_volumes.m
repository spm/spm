function gl = spm_get_volumes(P)
% Compute total volumes from tissue segmentations
% FORMAT gl = spm_get_volumes(P)
% P  - a matrix of image filenames
% gl - a vector of volumes (in litres)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2006-2022 Wellcome Centre for Human Neuroimaging


warning('spm:deprecated', ...
    ['spm_get_volumes will be removed in the future, please use ' ...
    'the Tissue Volumes Utility in the Batch interface, or:\n\t' ...
    'spm_summarise(P, ''all'', ''litres'')']);

if ~nargin
    [P,sts] = spm_select(Inf,'image','Select images');
    if ~sts, gl = []; return; end
end

V = spm_vol(P);
if spm_check_orientations(V, false)
    gl = spm_summarise(V, 'all', 'litres');
else
    N = numel(V);
    gl = nan(N, 1);
    for n = 1:N
        gl(n) = spm_summarise(V(n), 'all', 'litres');
    end
end
