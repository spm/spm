function spm_plot_ci(varargin)
% plots mean and conditional confidence intervals
% FORMAT spm_plot_ci(t,E,C,j,s)
% FORMAT spm_plot_ci(t,E,C,j)
% FORMAT spm_plot_ci(t,E,C)
% FORMAT spm_plot_ci(E,C)
%
% t - domain
% E - expectation
% C - variance or covariance
% j - indices of spm_cat(E(:)) to plot
% s - string to specify plot type
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_plot_ci.m 3715 2010-02-08 13:57:26Z karl $
 
% unpack
%--------------------------------------------------------------------------
if nargin == 5
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
    j = varargin{4};
    s = varargin{5};
elseif nargin == 4
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
    j = varargin{4};
elseif nargin == 3
    t = varargin{1};
    E = varargin{2};
    C = varargin{3};
else
    E = varargin{1};
    C = varargin{2};
end
 
if iscell(E)
    E = spm_cat(E(:));
end
try, j; catch, j = 1:size(E,1); end
try, t; catch, t = 1:size(E,2); end
try, s; catch, s = '';          end

 
% order and length of sequence
%--------------------------------------------------------------------------
E     = E(j,:);
[n N] = size(E);
 
% unpack conditional covariances
%--------------------------------------------------------------------------
ci    = spm_invNcdf(1 - 0.05);
try
    for i = 1:N
        c(:,i) = ci*sqrt(diag(C{i}(j,j)));
    end
catch
    c = ci*sqrt(C(j,j));
end

 
% conditional covariances
%--------------------------------------------------------------------------
if N > 1
    fill([t fliplr(t)],[full(E + c) fliplr(full(E - c))],...
        [1 1 1]*.8,'EdgeColor',[1 1 1]*.5),hold on
    plot(t,E,s)
else
    
    % plot in state-space
    %--------------------------------------------------------------------
    try,  C = C{1};  end
 
    [x y] = ellipsoid(E(1),E(2),0,c(1),c(2),0,32);
    fill(x',y',[1 1 1]*.8,'EdgeColor',[1 1 1]*.5),hold on
    plot(E(1,1),E(2,1),'.','MarkerSize',16)
 
end
hold off
drawnow
