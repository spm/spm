function V = spm_write_plane(V,dat,n)
% Write a transverse plane of image data.
% FORMAT VO = spm_write_plane(V,dat,n)
% V   - data structure containing image information.
%       - see spm_vol for a description.
% dat - the two dimensional image to write.
% n   - the plane number (beginning from 1).
%
% VO  - (possibly) modified data structure containing image information.
%       It is possible that future versions of spm_write_plane may
%       modify scalefactors (for example).
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id$


if isfield(V,'n'), n = num2cell([n V.n]); else, n = {n}; end;
S     = struct('type','()','subs',{{':',':',n{:}}});
V.private.dat = subsasgn(V.private.dat,S,dat);
