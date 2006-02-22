function [varargout] = spm_DEM_set(DEM)
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
    DEM.M = spm_DEM_M_set(DEM.M);
catch
    errordlg('please check your model')
end
try
    N     = size(DEM.Y,2);
catch
    errordlg('please specify data')
end
try
    DEM.class;
catch
    DEM.class = 'unknown';
end
 
% ensure model and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.Y,1) ~= DEM.M(1).l
    errordlg('DCM outputs and data are incompatible')
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

% ensure parameters and casues are not all zero
% -------------------------------------------------------------------------
if ~any(spm_vec(DEM.M.P)) && ~any(spm_vec(DEM.U))
    warndlg('please intialise parameters or inputs')
end

% unpack DEM if necessary
% -------------------------------------------------------------------------
if nargout == 4
    varargout{1} = DEM.M;
    varargout{2} = DEM.Y;
    varargout{3} = DEM.U;
    varargout{4} = DEM.X;
else
    varargout{1} = DEM;
end

