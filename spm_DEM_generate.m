function [DEM] = spm_DEM_generate(M,U,P,h)
% Generates data for a HDM
% FORMAT [DEM] = spm_DEM_generate(M,N,[P,h]): N-samples using z
% FORMAT [DEM] = spm_DEM_generate(M,U,[P,h]): size(U,2) samples using U
%
% M(i)     - HDM
% U(n x N} - causes or number of casues
% P{i}     - model-parameters for level i (defaults to M.P)
% H{i}     - hyper-parameters for level i (defaults to M.h)
%
% generates
% DEM.M    - hierarchical model (with M.P and M.h removed)
% DEM.Y    - responses or data
%
% and true causes NB: v{end} = U or z{end} (last level innovations)
% DEM.pU.v 
% DEM.pU.x
% DEM.pU.e
% DEM.pP.P
% DEM.pH.h
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$
 
% sequence length specified by priors on causes
%--------------------------------------------------------------------------
M  = spm_DEM_M_set(M);
try
    if length(U) > 1
        N = size(U,2);
    else
        N = U;
    end
catch
    warndlg('Please specify model inputs U or their number')
    return
end
 
% initialise model-parameters if specified
%--------------------------------------------------------------------------
try
    if ~iscell(P)
        errordlg('please ensure parameters are a cell array')
        return
    end
    for i = 1:length(P)
        M(i).P = spm_unvec(P{i},M(i).P);
    end
end
 
% initialise hyper-parameters if specified
%--------------------------------------------------------------------------
try
    for i = 1:length(h)
        M(i).h = h{i};
    end
end
 
% create innovations & add causes
%--------------------------------------------------------------------------
z        = spm_DEM_z(M,N);
if length(U) > 1
    z{end} = U;
end
 
% integrate HDM to obtain causal (v) and hidden states (x)
%--------------------------------------------------------------------------
[v,x]    = spm_DEM_int(M,z);
 
% Fill in DEM with response and its causes
%--------------------------------------------------------------------------
DEM.Y    = v{1};
DEM.pU.v = v;
DEM.pU.x = x;
DEM.pU.e = z;
DEM.pP.P = {M.P};
DEM.pH.h = {M.h};
 
% remove true values from M
%--------------------------------------------------------------------------
M        = rmfield(M,'P');
M        = rmfield(M,'h');
DEM.M    = M;
 


