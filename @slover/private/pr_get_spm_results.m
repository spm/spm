function [XYZ, Z, M] = pr_get_spm_results;
% fetches SPM results, returns as point list
% FORMAT [XYZ, Z, M] = pr_get_spm_results;
% 
% Outputs
% XYZ    - XYZ point list in voxels (empty if not found)
% Z      - values at points in XYZ
% M      - 4x4 voxel -> world transformation matrix 
% 
% $Id: pr_get_spm_results.m,v 1.2 2005/05/06 22:57:40 matthewbrett Exp $ 
  
errstr = '''Cannot find SPM results in workspace''';
[XYZ Z M] = deal([]);

V = spm('ver');
switch V(4:end)
 case '99'
  have_res = evalin('base', 'exist(''SPM'', ''var'')');
  if ~have_res, return, end
  SPM = evalin('base', 'SPM', ['error(' errstr ')']);
  XYZ = SPM.XYZ;
  Z   = SPM.Z;
  M   = evalin('base', 'VOL.M', ['error(' errstr ')']);
 otherwise
  have_res = evalin('base', 'exist(''xSPM'', ''var'')');
  if ~have_res, return, end
  xSPM = evalin('base', 'xSPM', ['error(' errstr ')']);
  XYZ = xSPM.XYZ;
  Z   = xSPM.Z;
  M   = xSPM.M;
end
