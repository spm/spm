function [S L] = spm_salience_map(M,n)
% creates a salience map
% FORMAT [S L] = spm_salience_map(M,n)
%
% S  - Salience (n x n,1)
% L  - list of (fictive) hidden control states (range of S)
%
% M  - generative model (with M(2).v and M(1).xo encoding location (L)
% n  - dimension of map (S)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_salience_map.m 4595 2011-12-19 13:06:22Z karl $
 
 
% Defaults
%==========================================================================
try, n; catch, n = 32; end             % dimension of salience map
global STIM                            % get stimulus hypotheses
 
% Salience map
%--------------------------------------------------------------------------
M     = spm_DEM_M_set(M);
DIM   = STIM.V.dim;
X     = linspace(-DIM(1)/2,DIM(1)/2,n);
Y     = linspace(-DIM(2)/2,DIM(2)/2,n);
[Y,X] = meshgrid(X,Y);
L     = [X(:) Y(:)]'/16;
 
% compute salience
%==========================================================================
for i = 1:length(L)
    M(1).x.o = L(:,i);                 % fictive hidden state
    qP       = spm_DEM_qC(M);
    S(i,1)   = spm_logdet(qP)/2;
end
