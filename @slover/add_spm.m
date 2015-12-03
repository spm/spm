function obj = add_spm(obj,xSPM)
% Add SPM blobs as new img to object, split effect, 'hot' colormap
% FORMAT obj = add_spm(obj)
%
% SPM results are fetched from the workspace
%__________________________________________________________________________

% Matthew Brett
% $Id: add_spm.m 6623 2015-12-03 18:38:08Z guillaume $

if nargin == 2
    XYZ = xSPM.XYZ;
    Z   = xSPM.Z;
    M   = xSPM.M;
else
    [XYZ,Z,M] = pr_get_spm_results;
end
if isempty(XYZ)
    warning('slover:noSPM', 'No SPM results to add');
    return
end

newimg = length(obj.img)+1;
obj.img(newimg).vol   = pr_blobs2vol(XYZ, Z, M);
obj.img(newimg).type  = 'split';
obj.img(newimg).prop  = 1;
obj.img(newimg).cmap  = hot;
obj.img(newimg).range = [0 max(Z)];
obj.cbar = [obj.cbar newimg];
