function [DEM] = spm_ADEM_set(DEM)
% Performs checks on DEM structures for active inversion
% FORMAT [DEM] = spm_ADEM_set(DEM)
%
% DEM.G  - generative model
% DEM.M  - recognition model
% DEM.C  - causes
% DEM.U  - explanatory variables, inputs or prior expectation of causes
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Karl Friston
% $Id: spm_ADEM_set.m 2521 2008-12-02 19:49:39Z karl $
 
% check recognition model
% -------------------------------------------------------------------------
try
    DEM.M = spm_DEM_M_set(DEM.M);
catch
    errordlg('please check your inversion model')
end
try
    DEM.G = spm_ADEM_M_set(DEM.G);
catch
    errordlg('please check your generative model')
end

% check data or generative model
% -------------------------------------------------------------------------
try
    N     = length(DEM.C);
catch
    errordlg('please specify causes')
end
try
    DEM.class;
catch
    DEM.class = 'active';
end
 

% transpose causes and confounds, if specified in conventional fashion
%--------------------------------------------------------------------------
try, if size(DEM.U,2) < N, DEM.U = DEM.U';    end, end
try, if size(DEM.C,2) < N, DEM.C = DEM.C';    end, end


% ensure model and data dimensions check
% -------------------------------------------------------------------------
if DEM.G(1).l ~= DEM.M(1).l
    errordlg('models are incompatible in terms of data')
end
if size(DEM.C,1) ~= DEM.G(end).l
    errordlg('DEM.G and causes (C) are incompatible')
end

% Default priors and confounds
% -------------------------------------------------------------------------
n  = DEM.M(end).l;
if ~isfield(DEM,'U')
    DEM.U = sparse(n,N);
end



% check prior expectation of causes (at level n) and confounds
%--------------------------------------------------------------------------
if ~nnz(DEM.U), DEM.U = sparse(n,N); end
if ~nnz(DEM.C), DEM.C = sparse(n,N); end
 
% ensure inputs and cause dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,1) ~= DEM.M(end).l
    errordlg('DEM.M and priors (U) are incompatible')
end
 
% ensure causes and data dimensions check
% -------------------------------------------------------------------------
if size(DEM.U,2) < N
    errordlg('priors and data have different lengths')
end

% check length of time-series
%--------------------------------------------------------------------------
if N < DEM.M(1).E.n
    errordlg('Please ensure time-series is longer than embedding order')
    return
end

