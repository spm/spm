function V = spm_map(P,ID)
% memory map of a volume image
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
% see also spm_map_vol.m
%
%__________________________________________________________________________
% %W% %E%

% ensure correct suffix for header filename
%-----------------------------------------------------------------------
P        = deblank(P);
Filename = P;
q        = length(P);
if P(q - 3) == '.'; P = P(1:(q - 4)); end

% get image descriptors and memory map
%-----------------------------------------------------------------------
if (nargin = 1)
	[DIM VOX SCALE TYPE OFFSET] = spm_hread([P '.hdr']);
	ID = [DIM(1:3) VOX(1:3) SCALE TYPE OFFSET];
end

V        = spm_map_vol(Filename,ID);
