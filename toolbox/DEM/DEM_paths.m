function DEM_paths
%--------------------------------------------------------------------------
% From densities to dynamics. This routine provides a numerical
% illustration of how to construct a random dynamical system from a
% probability density over some metric state-space. This involves sampling
% a set of points from an inhomogeneous density and creating an ordered
% list as follows: Add, without replacement, the point closest to the last
% point added. Start with the point furthest from the origin and repeat,
% until all points have been appended to the list. The ordered list can now
% be treated as a time-series, describing a unique path through
% state-space; namely, a time-homogenous Markov process. The implicit
% movement through state space can now be decomposed into a part that can
% be predicted given its position in state space, and a residual component,
% using a least squares estimator of the movement (i.e., flow) at each
% position. By construction, the residuals will have a zero mean Gaussian
% distribution by central limit theorem. This follows because the predicted
% flow is a linear mixture of some basis functions used for the least
% squares estimation.
%
% The ensuing decomposition can now be read as a random dynamical system,
% apt for analysis in terms of pullback attractors (Arnold, 2003; Crauel
% and Flandoli, 1994; Zaou et al., 2024) or, in a physics setting, a
% Langevin formulation (Girolami and Calderhead, 2011; Kerr and Graham,
% 2000; Pavliotis, 2014; Seifert, 2012; Sekimoto, 1998). In the current
% example, 2^12 points are sampled from a distribution over a
% two-dimensional state-space—with a flat Euclidean geometry—supplied. The
% functional form of the probability density, over the first state,
% comprises two Gaussian distributions. The conditional density of the
% second state was a quadratic function of the first.
% 
% The figure produced by this routine shows the sample points in the upper
% right panel, and the accompanying timeseries are shown on the upper
% right. These timeseries are rendered as a path--in the left middle
% panel--as a succession of vectors that we can associate with motion.
% Following least squares estimation of this motion--using fourth-order
% polynomial basis functions of the two states--the predicted motion (i.e.,
% flow) can be evaluated at each point (middle panel) and subtracted from
% the motion to identify the residuals (right middle panel), that can be
% treated as random fluctuations. Centring the fluctuations on their
% expectation of zero discloses their Gaussian form (lower panels). This
% construction associates time with some unique ordering of points sampled
% from a metric state-space. The ensuing dynamics can have complicated
% structure, as illustrated in the middle panel depicting flow as a vector
% field.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% set up
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
fontname(gcf,'Times New Roman');
rng('default')

% 2*N points in state-space
%--------------------------------------------------------------------------
N     = 2^12;
for i = 1:N
    x(1,i) = 1/2 + randn*sqrt(1/32);
    x(2,i) = x(1,i).^2 + randn*sqrt(1/8);

end
for i = (1:N) + N 
    x(1,i) = randn*sqrt(1/16) - 1/2;
    x(2,i) = x(1,i).^2 + randn*sqrt(1/16);
end
N     = size(x,2);

% plot
%--------------------------------------------------------------------------
subplot(3,2,1)
plot(x(1,:),x(2,:),'k.','MarkerSize',1/8)
title('Inhomogeneous points','FontSize',16)
xlabel('first state'), ylabel('second states'), axis square

% order points to minimise path length
%--------------------------------------------------------------------------
dx     = sum(x.^2);
[~,j]  = max(dx);
y(:,1) = x(:,j);
x(:,j) = [];
for i = 2:N
    dx     = sum((x - y(:,i - 1)).^2);
    [~,j]  = min(dx);
    y(:,i) = x(:,j);
    x(:,j) = [];
end

% plot
%--------------------------------------------------------------------------
subplot(3,2,2)
plot(y')
title('Time-homogeneous process','FontSize',16)
xlabel('time'), ylabel('states'), axis square, spm_axis tight


% decompose into a vector field and residuals
%==========================================================================
X     = y';
Y     = diff(X);
X     = X(1:(end - 2),:);
Y     = Y(1:(end - 1),:);
D     = [];
for i = 0:4
    D = [D, X.^i];
end

% polynomial basis functions
%--------------------------------------------------------------------------
d     = size(D,2);
for i = 1:d
    for j = 1:d
        D = [D, D(:,i).*D(:,j)];
    end
end

% least squares solution
%--------------------------------------------------------------------------
F     = D*(pinv(D)*Y);
W     = Y - F;

% plot
%--------------------------------------------------------------------------
subplot(3,3,4)
quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),'off','k')
title('Motion','FontSize',16)
xlabel('first state'), ylabel('Second states')
axis square, axis xy, a = axis;

% plot
%--------------------------------------------------------------------------
subplot(3,3,5)
quiver(X(:,1),X(:,2),F(:,1),F(:,2),'off','k')
title('Flow','FontSize',16)
xlabel('first state'), ylabel('Second states')
axis square, axis(a)

% plot
%--------------------------------------------------------------------------
subplot(3,3,6)
quiver(X(:,1),X(:,2),W(:,1),W(:,2),'off','k')
title('fluctuations','FontSize',16)
xlabel('first state'), ylabel('Second states')
axis square, axis(a)

% plot
%--------------------------------------------------------------------------
subplot(3,2,5)
sd  = std(W(:))*3;
plot(W(:,1),W(:,2),'k.','MarkerSize',1/8), hold on
title('Random fluctuations','FontSize',16)
xlabel('first dimension'), ylabel('Second dimension')
axis([-1,1,-1,1]*sd);
axis square, grid on

% plot
%--------------------------------------------------------------------------
subplot(3,2,6)
histogram(W,linspace(-sd,sd,64),'FaceColor',[1,1,1]/2)
title('Sample density','FontSize',16)
xlabel('fluctuations'), ylabel('Frequency')
axis square

return
