function obj = add_spm(obj)
% Adds SPM blobs as new img to object, split effect, 'hot' colormap
% FORMAT obj = add_spm(obj)
% 
% SPM results are fetched from the workspace
% 
% $Id: add_spm.m,v 1.1 2005/04/20 15:05:36 matthewbrett Exp $
  
[XYZ Z M] = pr_get_spm_results;
if isempty(XYZ)
  warning('slover:noSPM', 'No SPM results to add');
  return
end

newimg = length(obj.img)+1;
obj.img(newimg).vol = pr_blobs2vol(XYZ, Z, M);
obj.img(newimg).type = 'split';
obj.img(newimg).prop = 1;
obj.img(newimg).cmap = hot;
obj.img(newimg).range = [0 max(Z)];
obj.cbar = [obj.cbar newimg];
