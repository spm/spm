function [varargout] = spm_DEM_set(DEM)
% Performs checks on DEM structures
% FORMAT [DEM] = spm_DEM_set(DEM)
%
% DEM.M  - hierarchical model
% DEM.Y  - response variable, output or data
% DEM.U  - explanatory variables, inputs or prior expectation of causes
% DEM.X  - confounds
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_set.m 3695 2010-01-22 14:18:14Z karl $
 
% check recognition model
% -------------------------------------------------------------------------
try
    DEM.M = spm_DEM_M_set(DEM.M);
catch
    errordlg('please check your model')
end


% check format of inputs and data
% -------------------------------------------------------------------------
try, DEM.Y = DEM.Y.y'; end
try, DEM.U = DEM.U.u'; end


% check whether data are specified explicitly or with a generative model
% -------------------------------------------------------------------------
try
    N = size(DEM.Y,2);
catch
    try
        DEM.G = spm_DEM_M_set(DEM.G);
        N     = size(DEM.C,2);
    catch
        errordlg('please specify data or inputs')
    end
end
try
    DEM.class;
catch
    DEM.class = 'unknown';
end
 
% ensure model and data dimensions check
% -------------------------------------------------------------------------
try
    if size(DEM.Y,1) ~= DEM.M(1).l
        errordlg('DCM and data are incompatible')
    end
catch
    if size(DEM.C,1) ~= DEM.M(end).l
        errordlg('DCM and causes are incompatible')
    end
end

% Default priors and confounds
% -------------------------------------------------------------------------
n  = DEM.M(end).l;
if ~isfield(DEM,'U')
    DEM.U = sparse(n,N);
end
if ~isfield(DEM,'X')
    DEM.X = sparse(0,N);
end

% transpose causes and confounds, if specified in conventional fashion
%--------------------------------------------------------------------------
if size(DEM.U,2) < N, DEM.U = DEM.U';    end
if size(DEM.X,2) < N, DEM.X = DEM.X';    end

% check prior expectation of causes (at level n) and confounds
%--------------------------------------------------------------------------
if ~nnz(DEM.U), DEM.U = sparse(n,N); end
if ~nnz(DEM.X), DEM.X = sparse(0,N); end
 
% ensure inputs and cause dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,1) ~= DEM.M(end).l
    errordlg('DCM inputs and priors are not compatible')
end
 
% ensure causes and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,2) < N
    errordlg('priors and data have different lengths')
end
 
% ensure confounds and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.X,2) < N
    errordlg('confounds and data have different lengths')
end

% check length of time-series
%--------------------------------------------------------------------------
if N < DEM.M(1).E.n
    errordlg('Please ensure time-series is longer than embedding order')
    return
end

% unpack DEM if necessary
% -------------------------------------------------------------------------
if nargout > 1
    varargout{1} = DEM.M;
    varargout{2} = DEM.Y;
    varargout{3} = DEM.U;
    varargout{4} = DEM.X;
else
    varargout{1} = DEM;
end

