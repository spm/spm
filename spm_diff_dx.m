function [DX] = spm_diff_dx(varargin)
% optimisation of finite difference for numerical differentiation
% FORMAT [dx] = spm_diff_dx(f,x,...,n)
% FORMAT [dx] = spm_diff_dx(f,x,...,n,V)
% FORMAT [dx] = spm_diff_dx(f,x,...,n,'q')
%
% f      - [inline] function f(x{1},...)
% x      - input argument[s]
% n      - arguments to differentiate w.r.t.
%
% dx     - 'best' step size
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_diff_dx.m 7143 2017-07-29 18:50:38Z karl $


% Stability of numerical gradients
%==========================================================================
global GLOBAL_DX

% line search of step sizes
%--------------------------------------------------------------------------
dh    = 1;
dx    = -8:dh:0;
nd    = numel(dx);
for i = 1:nd
    GLOBAL_DX = exp(dx(i));
    dxdp{i}   = spm_diff(varargin{:});
end
np    = size(dxdp{1},2);
for i = 1:(nd - 1)
    for j = 1:np
        dgdh     = spm_vec(dxdp{i}{j}) - spm_vec(dxdp{i + 1}{j});
        ssd(i,j) = mean(abs(dgdh/dh).^2);
    end
end

% graphics
%--------------------------------------------------------------------------
dx    = dx(1:end - 1);
mss   = mean(log(ssd),2);
[~,j] = min(mss);
DX    = dx(j);

% graphics
%--------------------------------------------------------------------------
if nargout, return, end

str   = sprintf('Log stability: dx = exp(%.1f)',DX);

subplot(2,2,1), plot(dx,mss,dx,log(ssd),':')
title('Log stability','Fontsize',16), xlabel('log(dx)')
axis square, spm_axis tight

subplot(2,2,2), imagesc(1:np,dx,log(ssd))
title('Log stability','Fontsize',16), xlabel('Paramter mode')
axis square, ylabel('log(dx)')

subplot(2,2,3), plot(log(ssd(j,:)))
title(str,'Fontsize',16), xlabel('Paramter mode')
axis square, spm_axis tight

if iscell(varargin{end})
    V = varargin{end};
    subplot(2,2,4), imagesc(V{1})
    title('Paramter modes','Fontsize',16), xlabel('Paramter mode')
    axis square, ylabel('parameter')
end



