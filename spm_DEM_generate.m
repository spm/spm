function [DEM] = spm_DEM_generate(M,U,P,h,g)
% Generates data for a Hierarchical Dynamic Model (HDM)
% FORMAT [DEM] = spm_DEM_generate(M,N,P,h,g): N-samples using z
% FORMAT [DEM] = spm_DEM_generate(M,U,P,h,g): size(U,2) samples using U
%
% M(i)     - HDM
% U(n x N} - causes or N number of causes
% P{i}     - model-parameters for level i (defaults to M.pE)
% h{i}     - log-precisions   for level i (defaults to 32 - no noise)
% g{i}     - log-precisions   for level i (defaults to 32 - no noise)
%
% generates
% DEM.M    - hierarchical model (checked)
% DEM.Y    - responses or data
%
% and true causes NB: v{end} = U or z{end} (last level innovations)
% DEM.pU.v 
% DEM.pU.x
% DEM.pU.e
% DEM.pP.P
% DEM.pH.h
%
% NB: [lower bound on] random fluctuations will default to unit variance if
% not specified in M(i).V and M(i).W
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_generate.m 7322 2018-05-31 09:47:15Z karl $
 
% please parametric form in model if necessary
%--------------------------------------------------------------------------
% m     = numel(M);
% for i = 1:m
%     if ~isfield(M(i),'pE')
%         try, M(i).pE = P{i}; end
%     end
% end

% sequence length specified by priors on causes
%--------------------------------------------------------------------------
M     = spm_DEM_M_set(M);
DEM.M = M;
try
    if length(U) > 1
        N = size(U,2);
    else
        N = U;
        U = sparse(M(end).l,N);
    end
catch
    warndlg('Please specify model inputs U or number of samples')
    return
end
 
% initialize model-parameters if specified
%--------------------------------------------------------------------------
try, P; if ~iscell(P), P = {P}; end, catch, P = {M.pE}; end
try, h; if ~iscell(h), h = {h}; end, catch, h = {}; end
try, g; if ~iscell(P), g = {g}; end, catch, g = {}; end
 
% transcribe parameters and hyperparameters into prior expectations
%--------------------------------------------------------------------------
m     = numel(M);
for i = 1:m
    try
        M(i).pE = spm_unvec(P{i},M(i).pE);
    catch
        try
            M(i).pE = P{i};
        end
    end
end
for i = 1:m
    try
        M(i).hE = h{i};
    catch
        M(i).hE = (M(i).hE - M(i).hE) + 32;
    end
end
for i = 1:m
    try
        M(i).gE = g{i};
    catch
        M(i).gE = (M(i).gE - M(i).gE) + 32;
    end
end
 
% re-set M and create innovations
%--------------------------------------------------------------------------
M     = spm_DEM_M_set(M);
DEM.G = M;
[z,w] = spm_DEM_z(M,N);
 
% place exogenous causes in cell array
%--------------------------------------------------------------------------
for i = 1:m
    u{i} = sparse(M(i).l,N);
end
u{m}  = U;
 
 
% integrate HDM to obtain causal (v) and hidden states (x)
%--------------------------------------------------------------------------
[v,x,z,w] = spm_DEM_int(M,z,w,u);
 
% Fill in DEM with response and its causes
%--------------------------------------------------------------------------
DEM.Y     = v{1};
DEM.pU.v  = v;
DEM.pU.x  = x;
DEM.pU.z  = z;
DEM.pU.w  = w;
DEM.pP.P  = {M.pE};
DEM.pH.h  = {M.hE};
DEM.pH.g  = {M.gE};
