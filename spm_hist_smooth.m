function x = spm_hist_smooth(x,s)
% histogram smoothing
% FORMAT x = spm_hist_smooth(x,s)
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Karl Friston
% $Id: spm_hist_smooth.m 7934 2020-08-19 09:34:35Z karl $

% remove negative values
%--------------------------------------------------------------------------
x    = x(:);
i    = x < 0;
x(i) = 0;

% remove spikes
%--------------------------------------------------------------------------
dx   = gradient(x);
i    = abs(dx) > 4*std(dx);
x(i) = 0;

% graph Laplacian smoothing
%--------------------------------------------------------------------------
n    = numel(x);
K    = spm_speye(n,n,-1) - 2*spm_speye(n,n,0) + spm_speye(n,n,1);
K(1) = -1; K(end) = -1;
K    = spm_speye(n,n,0) + K/4;
K    = K^(s*4);
x    = K*x;