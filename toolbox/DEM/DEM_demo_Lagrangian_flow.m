function [Ep,Cp,Eh,F] = DEM_demo_Lagrangian_flow
% Demo to illustrate divergence and curl free flow specified by a 
% Lagrangian and antisymmetric matrices. This example uses a double well 
% potential and Newtonian dynamics.


% flow (to be matched)
%==========================================================================
f = @(x,u,P,M) [x(2,:); -(1/2)*x(1,:)];


% target flow
%--------------------------------------------------------------------------
x    = 8*rand(2,1024) - 4;
U    = x;
Y.y  = f(x);

% Lagrangian and Q
%==========================================================================
P.Q    = [0 1; 0 0];
P.P{1} = 0;
P.P{2} = [0; 0];
P.P{3} = [1/2 0; 0 1];

L    = @(x,P) P.P(1)*x(1,:).^2 + P.P(2)*x(2,:).^2;
dLdx = @(x,P) [2*P.P(1)*x(1,:);
               2*P.P(2)*x(2,:)];
           
% estimate parameters of Lagrangian and Q
%==========================================================================
M.IS = @(P,M,U)([0 P.Q;-P.Q 0]*dLdx(U,P));
M.IS = @(P,M,U) spm_fx_Lagrangian(P,M,U);
M.FS = @(y) spm_vec(y);
M.pE = P;
M.pC = diag(~~spm_vec(P));
M.hE = 8;
M.hC = 1/128;

[Ep,Cp,Eh,F] = spm_nlsi_GN(M,U,Y);


% Show results
%==========================================================================
 
% create surface
%--------------------------------------------------------------------------
N     = 64;
x     = -linspace(-8,8,N);
for i = 1:N
    for j = 1:N
        Lx(i,j) = L([x(i); x(j)],Ep);
        Px(i,j) = exp(-Lx(i,j));
    end
end
 

% and render
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
set(gcf,'Renderer','OpenGL')
colormap pink
 
hL = subplot(2,2,1);
surfc(x,x,Lx)
xlabel('velocity')
ylabel('state')
title('Lagrangian','FontSize',16)
shading interp
axis tight square
alpha(.5)
hold on
 
hP = subplot(2,2,2);
surfl(x,x,Px)
shading interp
xlabel('velocity')
ylabel('state')
title('Ergodic density','FontSize',16)
axis tight square
hold on
 
hC = subplot(2,1,2);
surfc(x,x,Lx)
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
G.x  = [0; -4];
U.dt = 1/32;
U.u  = sparse(512,0);

G.f  = @(x,u,P,M)([0 P.Q;-P.Q 0]*dLdx(x,P));
y    = spm_int_J(P,G,U);

 
% Superimpose on sufaces
%--------------------------------------------------------------------------
subplot(hL)
plot3(y(:,2),y(:,1),L(y',Ep),'r')
subplot(hC)
plot3(y(:,2),y(:,1),L(y',Ep),'r')
subplot(hP)
plot3(y(:,2),y(:,1),exp(-L(y',Ep)),'r')

% path
%--------------------------------------------------------------------------
G.f  = f;
y    = spm_int_J(P,G,U);

% Superimpose on sufaces
%--------------------------------------------------------------------------
subplot(hL)
plot3(y(:,2),y(:,1),L(y',Ep),'b')
subplot(hC)
plot3(y(:,2),y(:,1),L(y',Ep),'b')
subplot(hP)
plot3(y(:,2),y(:,1),exp(-L(y',Ep)),'b')

