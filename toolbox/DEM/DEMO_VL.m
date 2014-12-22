% function DEMO_VL
% evaluates the free energy landscape around the posterior
% FORMAT: spm_dcm_local_minima(DCM)
% DCM - (invert) model structure
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_VL.m 6291 2014-12-22 11:15:19Z karl $

% Model specification
%==========================================================================
clear; rng('default');

% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = 0;
A{2}  = 0;
A{3}  = 0;
B     = {};
C     = 1;

M.dipfit.model = 'CMC';
M.dipfit.type  = 'LFP';
M.dipfit.Nc    = 1;
M.dipfit.Ns    = 1;

% get priors
%--------------------------------------------------------------------------
[pE pC] = spm_dcm_neural_priors(A,B,C,M.dipfit.model);  % neuronal priors
[pE pC] = spm_L_priors(M.dipfit,pE,pC);                 % spatial  priors
[pE pC] = spm_ssr_priors(pE,pC);                        % spectral priors

[x,f]   = spm_dcm_x_neural(pE,M.dipfit.model);

rC      = spm_zeros(pC);
rC.G(1) = 1/4;
rC.G(2) = pC.G(2);

% create LFP model
%--------------------------------------------------------------------------
M.IS  = @spm_gen_erp;
M.f   = f;
M.g   = @spm_gx_erp;
M.x   = x;
M.pE  = pE;
M.pC  = rC;
M.hE  = 6;
M.hC  = 1/4;
M.m   = 1;
M.n   = numel(x);
M.l   = 1;



% Integrate system to see response (time-frequency and hemodynamic)
%==========================================================================
spm_figure('GetWin','Figure 1');clf

PP      = pE;
PP.G(1) = 2;
PP.G(2) = 0.5;


M.ons  = 64;
M.dur  = 16;
N      = 64;
U.dt   = 4/1000;
t      = (1:N)*U.dt;
pst    = t*1000;
U.u    = spm_erp_u(t,PP,M);
Yy     = spm_gen_erp(PP,M,U);
Y.y{1} = Yy{1} + mean(std(Yy{1}))*randn(size(Yy{1}))/4;
Y.dt   = U.dt;

% input
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(pst,U.u)
axis square
xlabel('time (ms)')
title('Exogenous input')

% LFP
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst,Y.y{1},':',pst,Yy{1})
axis square
xlabel('time (ms)')
title('LFP response')
drawnow

% return

% Stability analysis (over excitatory and inhibitory time constants)
%==========================================================================
fprintf('Stability analysis - please wait\n')
n     = 32;
P1    = linspace(0,3,n);
P2    = linspace(0,3,n);
for i = 1:length(P1)
    for j = 1:length(P2)
        P        = M.pE;
        P.G(:,1) = P1(i);
        P.G(:,2) = P2(j);
        S        = spm_ssm2s(P,M,0);
        S        = S(abs(imag(S)) > 2*pi & abs(imag(S)) < 128*2*pi);
        try
            [s k]   = max(real(S));
            HZ(i,j) = abs(imag(S(k)))/(2*pi);
            LE(i,j) = s;
        end
        
        M.Nmax         = 1;
        M.nograph      = 1;
        M.P            = P;
        [Qp,Cp,Eh,F,L] = spm_nlsi_GN(M,U,Y);
        FE(i,j)        = F;
        LL(i,j)        = L(1);
        Dp(i,j)        = L(2);
        Dh(i,j)        = L(3);
    end
end

% occam's window
%--------------------------------------------------------------------------
WIN  = 64;
fe   = FE - max(max(FE)); fe(fe < -WIN) = -WIN;
ll   = LL - max(max(LL)); ll(ll < -WIN) = -WIN;
dp   = Dp - max(max(Dp)); dp(dp < -WIN) = -WIN;
dh   = Dh - max(max(Dh)); dh(dh < -WIN) = -WIN;

[i j] = find(FE == max(FE(:)));
qP(1) = P1(i);
qP(2) = P2(j);



% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');clf

subplot(3,2,1)
imagesc(P1,P2,HZ), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold on
plot(qP(2),qP(1),'.g','MarkerSize',32), hold off
axis square
ylabel('inhibitory time constant (ms)')
xlabel('excitatory time constant (ms)')
title('Frequency','FontSize',16)

subplot(3,2,2)
imagesc(P1,P2,LE), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold off
axis square
ylabel('inhibitory time constant (ms)')
xlabel('excitatory time constant (ms)')
title('Stability','FontSize',16)

subplot(3,2,3)
imagesc(P1,P2,fe), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold off
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('Free energy','FontSize',16)

subplot(3,2,4)
imagesc(P1,P2,ll), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold off
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('Log likelihood','FontSize',16)

subplot(3,2,5)
imagesc(P1,P2,dp), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold off
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('Parameters','FontSize',16)

subplot(3,2,6)
imagesc(P1,P2,dh), hold on, 
plot(PP.G(2),PP.G(1),'.r','MarkerSize',32), hold off
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('Hyperparameters','FontSize',16)



