function [DEM] = spm_DEM_generate(M,U,P,H)
% Generates data for a HDM
% FORMAT [DEM] = spm_DEM_generate(M,U,P,H)
%
% M      - HDM
% U      - causes
% P{i}   - model-parameters for level i (defaults t0 M.P)
% H{i}   - hyper-parameters for level i (defaults t0 M.h)
%
% generates
% DEM.M  - hierarchical model (with M.P and M.h removed)
% DEM.Y  - responses or data
%
% and true causes NB: v{end} = e{end} + U; e = z = innovations
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
M       = spm_M_set(M);
try
    N   = size(U,2);
catch
    warndlg('Please specify model inputs U')
    return
end
 
% initialise model-parameters if specified
%--------------------------------------------------------------------------
try
    for i = 1:length(P)
        M(i).P = P{i};
    end
end
 
% initialise hyper-parameters if specified
%--------------------------------------------------------------------------
try
    for i = 1:length(H)
        M(i).h = H{i};
    end
end
 
% create innovations & add causes
%--------------------------------------------------------------------------
z         = spm_DEM_z(M,N);
z{end}    = z{end} + U;
 
% integrate HDM to obtain causal (v) and hidden states (x)
%--------------------------------------------------------------------------
[v,x]     = spm_DEM_int(M,z);
z{end}    = z{end} - U;
 
% Fill in DEM with response and its causes
%--------------------------------------------------------------------------
DEM.Y     = v{1};
DEM.pU.v  = v;
DEM.pU.x  = x;
DEM.pU.e  = z;
DEM.pP.P = {M.P};
DEM.pH.h = {M.h};
 
% remove true values from M
%--------------------------------------------------------------------------
M         = rmfield(M,'P');
M         = rmfield(M,'h');
DEM.M     = M;
 


