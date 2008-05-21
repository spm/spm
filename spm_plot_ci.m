function spm_plot_ci(varargin)
% plots mean and conditional confidence intervals
% FORMAT spm_plot_ci(t,E,C)
% FORMAT spm_plot_ci(E,C)
%
% t - domain
% E - expectation
% C - variance or covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_plot_ci.m 1704 2008-05-21 14:00:09Z karl $
 
% unpack
%--------------------------------------------------------------------------
try
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
catch
    E = varargin{1};
    C = varargin{2};
    t = [1:size(E,2)];
end
 
 
% order and length of sequence
%--------------------------------------------------------------------------
[n N] = size(E);
 
% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);
try
    for i = 1:N
        c(:,i) = ci*sqrt(diag(C{i}));
    end
catch
    c = ci*sqrt(C);
end
 
% conditional covariances
%--------------------------------------------------------------------------
if N > 1
    fill([t fliplr(t)],[full(E + c) fliplr(full(E - c))],...
        [1 1 1]*.8,'EdgeColor',[1 1 1]*.8),hold on
    plot(t,E)
else
    
    % plot in state-space
    %--------------------------------------------------------------------
    try,  C = C{1};  end
 
    [x y] = ellipsoid(E(1),E(2),0,c(1),c(2),0,32);
    fill(x',y',[1 1 1]*.8,'EdgeColor',[1 1 1]*.8),hold on
    plot(E(1,1),E(2,1),'.','MarkerSize',16)
 
end
hold off
drawnow
