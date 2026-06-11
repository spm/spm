function DEM_reverse_paths
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

%% dynamics and parameters of a Lorentz system (with Jacobian)
%==========================================================================
% dxdt = f(x) + w:  see notes at the end of this script
%--------------------------------------------------------------------------
f    = @(x,v,P,M) [P(1)*x(2) - P(1)*x(1);
                   P(3)*x(1) - x(2) - x(1)*x(3);
                  -P(2)*x(3) + x(1)*x(2)]/64;
J    = @(x,v,P,M) [[     -P(1),P(1),     0];
                   [P(3) - x(3), -1, -x(1)];
                   [     x(2), x(1), -P(2)]]/64;
P    = [10; 8/3; 32];                  % parameters [sig, beta, rho]
x0   = [1; 1; 24];                     % initial state
W    = 1/16;                           % precision of random fluctuations

disp('Lyapunov dimnesion')
disp(3 - 2*(P(1) + P(2) + 1)/(P(1) + 1 + sqrt((P(1) - 1)^2 + 4*P(1)*P(3))))

% state-space model (for SPM integrators)
%--------------------------------------------------------------------------
M.f  = f;
M.J  = J;
M.g  = @(x,v,P,M) x;
M.x  = x0;
M.m  = 0;
M.pE = P;
M.W  = W;

% solution of stochastic and underlying ordinary differential equation
% -------------------------------------------------------------------------
spm_figure('GetWin','NESS (0)'); clf
N    = 2^10;
U.u  = sparse(N,M(1).m);
U.dt = 1;
t    = spm_int_ode(M.pE,M,U);
r    = spm_int_sde(M.pE,M,U);

subplot(3,1,1)
plot(t), hold on, set(gca,'ColorOrderIndex',1)
plot(r,':'), hold off
title('Trajectory','Fontsize',16)
xlabel('time'), ylabel('states')
spm_axis tight, box off

subplot(3,2,3)
% plot3(t(:,1),t(:,2),t(:,3)), hold on, set(gca,'ColorOrderIndex',1)
plot3(r(:,1),r(:,2),r(:,3),':'), hold off
title('State-space','Fontsize',16)
xlabel('1st state'), ylabel('2nd state'), zlabel('3rd state')
axis square, box off


% plot
%--------------------------------------------------------------------------
% subplot(3,2,1)
% plot(x(1,:),x(2,:),'k.','MarkerSize',1/8)
% title('Inhomogeneous points','FontSize',16)
% xlabel('first state'), ylabel('second states'), axis square

% re-order points to increase path length
%--------------------------------------------------------------------------
% dx     = sum(x.^2);
% [~,j]  = max(dx);
% y(:,1) = x(:,j);
% x(:,j) = [];
% for i = 2:N
%     dx     = sum((x - y(:,i - 1)).^2);
%     [~,j]  = min(dx);
%     y(:,i) = x(:,j);
%     x(:,j) = [];
% end
% 
% % plot
% %--------------------------------------------------------------------------
% subplot(3,2,2)
% plot(y')
% title('Time-homogeneous process','FontSize',16)
% xlabel('time'), ylabel('states'), axis square, spm_axis tight


% decompose into a vector field and residuals
%==========================================================================
x     = t';
L     = [];
B     = 0;
for t = 1:1280

    % order points to minimise path length
    %--------------------------------------------------------------------------
    j      = 1;
    y(:,1) = x(:,j);
    x(:,j) = [];
    for i = 2:N
        if t > 1
            D      = spm_basis(y(:,i - 1)');
            dx     = D*B;
        else
            dx = 0;
        end
        p      = y(:,i - 1) + dx(:);
        dx     = sum((x - p).^2);
        [~,j]  = min(dx);
        y(:,i) = x(:,j);
        if i < N
            x(:,j) = [];
        end
    end

    x     = y;
    X     = y';
    Y     = diff(X);
    X     = X(1:(end - 2),:);
    Y     = Y(1:(end - 1),:);
    D     = spm_basis(X);

    % least squares solution
    %----------------------------------------------------------------------
    B     = pinv(D)*Y;
    F     = D*B;
    W     = Y - F;
    L     = [L, sum(sqrt(sum(W.^2,2)))]


    
%     % plot
%     %--------------------------------------------------------------------
%     subplot(3,3,4)
%     quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),'off','k')
%     title('Motion','FontSize',16)
%     xlabel('first state'), ylabel('Second states')
%     axis square, axis xy, a = axis;
%     drawnow

end

return

function D = spm_basis(X)

D     = [];
for i = 0:4
    D = [D, X.^i];
end

% polynomial basis functions
%----------------------------------------------------------------------
d     = size(D,2);
for i = 1:d
    for j = 1:d
        D = [D, D(:,i).*D(:,j)];
    end
end

return

% decompose into a vector field and residuals
%==========================================================================
Z     = r;
L     = [];
for i = 1:128

    for j = 1:8
        r = randi(N,2,1);
        T = Z;
        T(r(2),:) = Z(r(1),:);
        T(r(1),:) = Z(r(2),:);
        Z = T;
    end

    X     = Z;
    Y     = diff(X);
    X     = X(1:(end - 2),:);
    Y     = Y(1:(end - 1),:);
    D     = [];
    for d = 0:4
        D = [D, X.^d];
    end

    % polynomial basis functions
    %----------------------------------------------------------------------
    d     = size(D,2);
    for i = 1:d
        for j = 1:d
            D = [D, D(:,i).*D(:,j)];
        end
    end

    % least squares solution
    %----------------------------------------------------------------------
    F     = D*(pinv(D)*Y);
    W     = Y - F;
    L     = [L, sum(sqrt(sum(W.^2,2)))];


    
%     % plot
%     %----------------------------------------------------------------------
%     subplot(3,3,4)
%     quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),'off','k')
%     title('Motion','FontSize',16)
%     xlabel('first state'), ylabel('Second states')
%     axis square, axis xy, a = axis;
%     drawnow

end


plot(L)

return
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
histogram(W(:,1),linspace(-sd,sd,64),'FaceColor',[1,1,1]/2)
title('Sample density','FontSize',16)
xlabel('fluctuations'), ylabel('Frequency')
axis square

return
