function spm_check_orientations(V)
% Check the dimensions and orientations of the images
% FORMAT spm_check_orientations(V)
%
% V - a structure as returned by spm_vol.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_check_orientations.m 253 2005-10-13 15:31:34Z guillaume $

if numel(V)<=1, return; end;

dims = cat(1,V.dim);
if any(any(diff(dims,1,1),1)),
    fprintf('The images do not all have the same dimensions.\n\n');
    for i=1:numel(V),
        fprintf('[%d %d %d]  %s\n',V(i).dim, V(i).fname);
    end;
    fprintf('\n');
    error('The dimensions must be identical for this procedure.');
end;

matx = reshape(cat(3,V.mat),[16,numel(V)]);
if any(any(abs(diff(matx,1,2))>1e-4)),
    fprintf('The images do not all have same orientation and/or voxel sizes.\n\n');
    for i=1:numel(V),
        fprintf('[%g %g %g %g; %g %g %g %g; %g %g %g %g]  %s\n',...
                V(i).mat(1:3,:)', V(i).fname);
    end;
    fprintf('\n');
    error('The orientations etc must be identical for this procedure.');
end;

