function gl = spm_get_volumes(P)
% Compute total volumes from tissue segmentations
% FORMAT gl = spm_get_volumes(P)
% P  - a matrix of image filenames
% gl - a vector of volumes (in litres)
%__________________________________________________________________________
% Copyright (C) 2006-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_volumes.m 5647 2013-09-20 13:03:44Z ged $

warning('spm:deprecated', ...
    ['spm_get_volumes will be removed in the future, please use ' ...
    'the Tissue Volumes Utility in the Batch interface, or:\n\t' ...
    'spm_summarise(P, ''all'', ''litres'')']);

if ~nargin
    [P,sts] = spm_select(Inf,'image','Select images');
    if ~sts, gl = []; return; end
end

gl = spm_summarise(P, 'all', 'litres');
