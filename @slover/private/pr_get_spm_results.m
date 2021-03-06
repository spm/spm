function [XYZ, Z, M] = pr_get_spm_results
% Fetch SPM results and return as point list
% FORMAT [XYZ, Z, M] = pr_get_spm_results
%
% Outputs
% XYZ    - XYZ point list in voxels (empty if not found)
% Z      - values at points in XYZ
% M      - 4x4 voxel -> world transformation matrix
%__________________________________________________________________________

% Copyright (C) 2005-2022 Matthew Brett


errstr = '''Cannot find SPM results in workspace''';
[XYZ,Z,M] = deal([]);

have_res = evalin('base', 'exist(''xSPM'', ''var'')');
if ~have_res, return, end
xSPM = evalin('base', 'xSPM', ['error(' errstr ')']);
XYZ = xSPM.XYZ;
Z   = xSPM.Z;
M   = xSPM.M;
