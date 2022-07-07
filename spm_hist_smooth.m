function x = spm_hist_smooth(x,s)
% Histogram smoothing (graph Laplacian)
% FORMAT x = spm_hist_smooth(x,s)
% x   - data vector
% s   - smoothing
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


% remove negative values
%--------------------------------------------------------------------------
x    = x(:);
i    = x < 0;
x(i) = 0;

% remove spikes
%--------------------------------------------------------------------------
dx   = gradient(x);
i    = abs(dx) > 16*std(dx);
x(i) = 0;

% graph Laplacian smoothing
%--------------------------------------------------------------------------
n    = numel(x);
K    = spm_speye(n,n,-1) - 2*spm_speye(n,n,0) + spm_speye(n,n,1);
K(1) = -1; K(end) = -1;
K    = spm_speye(n,n,0) + K/4;
K    = K^(s*4);
x    = K*x;
