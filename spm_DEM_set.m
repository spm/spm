function [DEM] = spm_DEM_set(DEM)
% Performs checks on DEM structures
% FORMAT [DEM] = spm_DEM_set(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - inputs or data
% DEM.U  - prior expectation of causes
% DEM.X  - observation confounds
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id$
 
% check model and data
% -------------------------------------------------------------------------
try
    DEM.M = spm_M_set(DEM.M);
catch
    msgbox('please check your model')
    error(' ')
end
try
    N     = size(DEM.Y,2);
catch
    msgbox('please specify data')
    error(' ')
end
try
    DEM.class;
catch
    DEM.class = 'unknown';
end
 
% ensure model and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.Y,1) ~= DEM.M(1).l
    errordlg('Outputs and data are not compatible')
end
 
% Default priors and confounds
% -------------------------------------------------------------------------
try
    if ~isfield(DEM,'U')
        DEM.U = sparse(DEM.M(end).l,N);
    end
end
try
    if ~isfield(DEM,'X')
        DEM.X = sparse(0,N);
    end
end
 
% ensure Inputs and cause dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,1) ~= DEM.M(end).l
    errordlg('Inputs and causes are not compatible')
end
 
% ensure causes and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,2) < N
    errordlg('causes and data not compatible')
end
 
% ensure confounds and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.X,2) < N
    errordlg('confounds and data not compatible')
end
