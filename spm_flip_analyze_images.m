function flip = spm_flip_analyze_images
% Do Analyze format images need to be left-right flipped?
%-----------------------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: spm_flip_analyze_images.m 112 2005-05-04 18:20:52Z john $


global defaults
if isempty(defaults) | ~isfield(defaults,'analyze') |...
         ~isfield(defaults.analyze,'flip')
	warning('Cant get default Analyze orientation - assuming flipped');
	flip = 1;
	return;
end;
flip = defaults.analyze.flip;
