function [F] = spm_DEM_t_fmin(n,r,s,dt,pt)
% FORMAT [F] = spm_DEM_t_fmin(n,r,s,dt,pt)
%__________________________________________________________________________
% n    - embedding order 
% d    - restriction order
% s    - temporal smoothness - s.d. of kernel {seconds}
% dt   - time interval {seconds}
% pt   - prediction interval
%
% F    - KL{P(0,n) P(pt,r)}
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

P    =  spm_DEM_P(n,r,s,dt,pt);
L    =  svd(full(P));
F    = -log(prod(L(1:r)))/2;
