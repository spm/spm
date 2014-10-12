function DEM_demo_Lagrangian
% Demo to illustrate divergence and curl free flow specified by a 
% Lagrangian and antisymmetric matrices. This example uses a double well 
% potential and Newtonian dynamics.

 
% Lagrangian and Q
%==========================================================================
L = @(x) (1/4)*x(1,:).^2 + (1/2)*(x(2,:)).^2;
L = @(x) ((x(1,:).^2 - 8).^2)/64 + (x(2,:).^2)/2;
Q = [0 1; -1 0];
R = [1 0;
     0 1]/32;
 
% create surface
%--------------------------------------------------------------------------
N     = 64;
x     = -linspace(-8,8,N);
for i = 1:N
    for j = 1:N
        G(i,j) = L([x(i); x(j)]);
        P(i,j) = exp(-L([x(i); x(j)]));
    end
end
 

% and render
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
set(gcf,'Renderer','OpenGL')
colormap pink
 
hL = subplot(2,2,1);
surfc(x,x,G)
xlabel('velocity')
ylabel('state')
title('Lagrangian','FontSize',16)
shading interp
axis tight square
alpha(.5)
hold on
 
hP = subplot(2,2,2);
surfl(x,x,P)
shading interp
xlabel('velocity')
ylabel('state')
title('Ergodic density','FontSize',16)
axis tight square
hold on
 
hC = subplot(2,1,2);
surfc(x,x,G)
xlabel('velocity')
ylabel('state')
title('Flow','FontSize',16)
shading interp
axis tight square
view([0 0 1])
alpha(.5)
hold on
 
 
% path
%--------------------------------------------------------------------------
N     = 4098;
dt    = 1/64;
x     = [0; -4];
for i = 1:N
    f        = (Q - R)*spm_diff(L,x,1)';
    x        = x + f*dt;
    p(:,i)   = x;
end
 
% Superimpose on sufaces
%--------------------------------------------------------------------------
subplot(hL)
plot3(p(2,:),p(1,:),L(p),'r')
subplot(hC)
plot3(p(2,:),p(1,:),L(p),'r')
subplot(hP)
plot3(p(2,:),p(1,:),exp(-L(p)),'r')


