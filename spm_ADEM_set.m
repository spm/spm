function [DEM] = spm_ADEM_set(DEM)
% Perform checks on DEM structures for active inversion
% FORMAT [DEM] = spm_ADEM_set(DEM)
%
% DEM.G  - generative model
% DEM.M  - recognition model
% DEM.C  - exogenous causes
% DEM.U  - prior expectation of causes
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ADEM_set.m 5962 2014-04-17 12:47:43Z spm $


% check recognition model
% -------------------------------------------------------------------------
try
    DEM.M = spm_DEM_M_set(DEM.M);
catch
    errordlg('please check your generative model')
end
try
    DEM.G = spm_ADEM_M_set(DEM.G);
catch
    errordlg('please check your generative process')
end

% check data or generative model
% -------------------------------------------------------------------------
try
    N     = length(DEM.C);
catch
    errordlg('please specify causes (e.g., sparse(0,N)')
end
try
    DEM.class;
catch
    DEM.class = 'active';
end
 
% Default priors
% -------------------------------------------------------------------------
if ~isfield(DEM,'U')
    DEM.U = sparse(DEM.M(end).l,N);
end

% transpose causes and confounds, if specified in conventional fashion
%--------------------------------------------------------------------------
try, if size(DEM.U,2) < N, DEM.U = DEM.U';    end, end
try, if size(DEM.C,2) < N, DEM.C = DEM.C';    end, end


% ensure model and input dimensions check
% -------------------------------------------------------------------------
if size(DEM.C,1) ~= DEM.G(end).l
    errordlg('model (G) and causes (C) are incompatible')
end
if size(DEM.U,1) ~= DEM.M(end).l
    errordlg('model (M) and priors (U) are incompatible')
end

% check prior expectation of causes (at level n) and confounds
%--------------------------------------------------------------------------
if ~nnz(DEM.U), DEM.U = sparse(DEM.M(end).l,N); end
if ~nnz(DEM.C), DEM.C = sparse(DEM.G(end).l,N); end
 

% ensure causes and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,2) < N
    errordlg('priors (U) and causes (C) have different lengths')
end

% check length of time-series
%--------------------------------------------------------------------------
if N < DEM.M(1).E.n
    errordlg('Please ensure time-series is longer than embedding order')
    return
end

