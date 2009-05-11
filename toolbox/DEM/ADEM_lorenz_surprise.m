% This demo simply computes the loss-function (negative reward for a Lorenz
% system) to show the mapping from value (expected reward) to reward is easy
% to compute
 
 
% generative model
%==========================================================================                        % switch for demo
G(1).E.s = 1/4;                        % smoothness
G(1).E.n = 6;                          % smoothness
G(1).E.d = 2;                          % smoothness
 
% dynamics
%--------------------------------------------------------------------------
fL      = '[v; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';
 
% parameters
%--------------------------------------------------------------------------
PL      = [10; -8/3; 32];
PM      = [0;     0;  0];
 
% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number
 
% level 1
%--------------------------------------------------------------------------
G(1).x  = [1; 1; 24];
G(1).f  = inline(fL ,'x','v','a','P');
G(1).g  = inline('x','x','v','a','P');
G(1).pE = PL;
G(1).V  = exp(8);                           % error precision
G(1).W  = exp(8);                           % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a  = [0;0;0];                          % action variables
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);
 
% space
%--------------------------------------------------------------------------
x{1}    = linspace(-32,32,64);
x{2}    = linspace(-32,32,64);
x{3}    = linspace(  4,64,64);
 
% Fokker-Planck operator and equilibrium density
%==========================================================================
[M0,q0] = spm_fp(G,x);
U.u     = sparse(1024,G(1).m);
t       = spm_int_J(G(1).pE,G,U);

 
% loss-function and negative surprise (value)
%--------------------------------------------------------------------------
L    = -spm_unvec(spm_vec(log(q0))'*M0,q0);
V    = log(q0);

% trim
%--------------------------------------------------------------------------
q0   = q0(2:end - 1,2:end - 1,2:end - 1);
L    =  L(2:end - 1,2:end - 1,2:end - 1);
V    =  V(2:end - 1,2:end - 1,2:end - 1);
x{1} = x{1}(2:end - 1);
x{2} = x{2}(2:end - 1);
x{3} = x{3}(2:end - 1);

% axes
%--------------------------------------------------------------------------
i    = 3;
j    = 1:3;
j(i) = [];
a    = [x{j(2)}(1) x{j(2)}(end) x{j(1)}(1) x{j(1)}(end)];
 
% surprise
%--------------------------------------------------------------------------
subplot(3,2,1)
imagesc(x{j(2)},x{j(1)},V(:,:,34))
hold on, plot(t(:,j(2)),t(:,j(1)),'r'), hold off
axis(a)
axis square xy
title('surprise','Fontsize',16)
 
% cost function
%--------------------------------------------------------------------------
subplot(3,2,2)
imagesc(x{j(2)},x{j(1)},L(:,:,34))
hold on, plot(t(:,j(2)),t(:,j(1)),'r'), hold off
axis(a)
axis square xy
title('loss','Fontsize',16)
 
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
