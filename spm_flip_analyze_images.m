function flip = spm_flip_analyze_images
% Do Analyze format images need to be left-right flipped?
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_flip_analyze_images.m 1131 2008-02-06 11:17:09Z spm $


global defaults
if isempty(defaults) | ~isfield(defaults,'analyze') |...
         ~isfield(defaults.analyze,'flip')
    warning('Cant get default Analyze orientation - assuming flipped');
    flip = 1;
    return;
end;
flip = defaults.analyze.flip;
