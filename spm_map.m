function V = spm_map(P,ID)
% memory map of a volume image - an obsolete function
% FORMAT V = spm_map(P,[ID])
% P - filename
% ID - Image descriptors passed to spm_map_vol
%____________________________________________________________________________
%
% spm_map returns a vector V identifying a memory mapped image
% volumne on disk.  Memory mapping avoids having very large objects
% in working memory.
% Information about the image is read from *.hdr if is exists, otherwise
% the global default variables are used.  
%
% SPM_MAP is obsolete and will be eliminated.  Please use spm_vol.m.
%
% see also spm_map_vol.m
%
%__________________________________________________________________________
% %W% %E%

warning('SPM_MAP is obsolete and will be eliminated.');

% ensure correct suffix for header filename
%-----------------------------------------------------------------------
P        = deblank(P);
Filename = P;
q        = length(P);
if P(q - 3) == '.'; P = P(1:(q - 4)); end

% get image descriptors and memory map
%-----------------------------------------------------------------------
if (nargin == 1)
	[DIM VOX SCALE TYPE OFFSET] = spm_hread([P '.hdr']);
	ID = [DIM(1:3) VOX(1:3) SCALE TYPE OFFSET];
end

err = [...
'f=spm_figure(''findwin'',''Graphics'');' ...
'if ~isempty(f),' ...
'	spm_figure(''Clear'',''Graphics'');' ...
'	spm_figure(''Clear'',''Interactive'');' ...
'	figure(f);' ...
'	ax=axes(''Visible'',''off'');' ...
'	text(0,0.5,[''There is a problem with: '' spm_str_manip(P,''a16'') ],''FontSize'',20);' ...
'end;' ...
'error([''There is a problem with: '' P]);'];

eval('V = spm_map_vol(Filename,ID);',err);
