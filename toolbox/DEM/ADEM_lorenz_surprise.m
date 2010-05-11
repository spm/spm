% This demo computes the cost-function (negative reward) for a Lorenz
% system; to show cost can be computed easily from value (negative 
% surprise or sojourn time). However, value is not a Lyapunov function
% because the flow is not curl-free (i.e., is not irrotational).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_lorenz_surprise.m 3884 2010-05-11 19:43:53Z karl $
 

% dynamics and parameters
%-------------------------------------------------------------------------
spm_figure('GetWin','Cost and surprise');
clear

DEMO     = 0;                          % switch for demo
LOR      = 1;                          % Lorenz vs a linear system
 
% generative model
%==========================================================================                       % switch for demo
G(1).E.s = 1/4;                        % smoothness
G(1).E.n = 6;                          % smoothness
G(1).E.d = 2;                          % smoothness
 
% dynamics and parameters
%--------------------------------------------------------------------------
if LOR
    fL  = '[v; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
    PL  = [10; -8/3; 32];
    x0  = [1; 1; 24];
    W   = exp(8);
else
    fL  = '[v; 0; 0] + [-P(1) 0 0; 0 -P(2) 0; 0 0 -P(3)]*x/64';
    PL  = [1; 1; 8];
    x0  = [-16; 16; 0];
    W   = diag([32; 32; exp(8)]);
end
 
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = x0;
G(1).f  = inline(fL ,'x','v','P');
G(1).g  = inline('x','x','v','P');
G(1).pE = PL;
G(1).V  = exp(8);                           % error precision
G(1).W  = W;                                % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = [0;0;0];                          % action variables
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_DEM_M_set(G);
 
% space
%--------------------------------------------------------------------------
N        = 66;
if LOR
    x{1} = linspace(-32,32,N);
    x{2} = linspace(-32,32,N);
    x{3} = linspace(  4,64,N);
else
    x{1} = linspace(-8,8,N);
    x{2} = linspace(-8,8,N);
    x{3} = linspace(-1,1,3);
end
 
% equilibrium density (q0), loss (L) and value (V) functions
%==========================================================================
if DEMO
    
    % Fokker-Planck operator and equilibrium density
    %----------------------------------------------------------------------
    [M0,q0] = spm_fp(G,x);
 
    % loss-function and negative surprise (value)
    %----------------------------------------------------------------------
    V    = log(q0);
    L    = spm_unvec(spm_vec(V)'*M0,q0);
 
    % trim
    %----------------------------------------------------------------------
    q0   = q0(2:end - 1,2:end - 1,2:end - 1);
    L    =  L(2:end - 1,2:end - 1,2:end - 1);
    V    =  V(2:end - 1,2:end - 1,2:end - 1);
    x{1} = x{1}(2:end - 1);
    x{2} = x{2}(2:end - 1);
    x{3} = x{3}(2:end - 1);
 
    if LOR
        save DEM_lorenz_suprise q0 L V x
    end
else
    load DEM_lorenz_suprise
end
 
 
% axes and trajectory
%--------------------------------------------------------------------------
if LOR
    z = 24;
else
    z = 1;
end
i    = 3;
j    = 1:3;
j(i) = [];
T    = 1024;
U.u  = sparse(T,G(1).m);
t    = spm_int_J(G(1).pE,G,U);
 
% surprise
%--------------------------------------------------------------------------
a    = [x{j(2)}(1) x{j(2)}(end) x{j(1)}(1) x{j(1)}(end)];
subplot(3,2,1)
imagesc(x{j(2)},x{j(1)},V(:,:,z))
hold on, plot(t(:,j(2)),t(:,j(1)),'r'), hold off
axis(a)
axis square xy
title('value','Fontsize',16)
 
% cost function
%--------------------------------------------------------------------------
subplot(3,2,2)
imagesc(x{j(2)},x{j(1)},L(:,:,z))
hold on, plot(t(:,j(2)),t(:,j(1)),'r'), hold off
axis(a)
axis square xy
title('cost','Fontsize',16)
 
% equilibrium density
%--------------------------------------------------------------------------
subplot(3,2,3)
imagesc(x{j(2)},x{j(1)},squeeze(max(q0,[],i)))
hold on, plot(t(:,j(2)),t(:,j(1)),'r'), hold off
axis(a)
axis square xy
title('density','Fontsize',16)
 
% exemplar trajectory
%--------------------------------------------------------------------------
subplot(3,2,4)
plot(t(:,j(2)),t(:,j(1)),'k')
axis(a)
axis square xy
title('trajectory','Fontsize',16)
 
% evaluate V(t)
%--------------------------------------------------------------------------
for i = 1:T
    [q ii] = min(abs(t(i,1) - x{1}));
    [q ij] = min(abs(t(i,2) - x{2}));
    [q ik] = min(abs(t(i,3) - x{3}));
    v(i)   = V(ii,ij,ik);
end
 
% and plot
%--------------------------------------------------------------------------
maxV = max(V(:));
subplot(3,1,3)
plot(v,'k'), hold on
plot([1 T],[maxV maxV],'r:'), hold off
title('value','Fontsize',16)
set(gca,'XLim',[1 T]);
box off
