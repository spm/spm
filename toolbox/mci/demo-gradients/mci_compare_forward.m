function [els_sun,els_ode,els_spm] = mci_compare_forward (model)
% Compare integration methods 
% FORMAT [els_sun,els_ode,els_spm] = mci_compare_forward (model)
%
% model     'phase', 'nmm-r2p2'
%
% Run integration 9 times - compare speed and accuracy
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_forward.m 6548 2015-09-11 12:39:47Z will $

[P,M,U,Y] = mci_compare_setup (model);

M.reltol=1e-2;
M.abstol=1e-4;

M = spm_mci_minit (M);
if isstruct(M.pC)
    pC=full(diag(spm_vec(M.pC)));
else
    pC = M.pC;
end

UI.u=U';
UI.dt=M.T/M.N;
        
Nsamp=9;
for s=1:Nsamp,
    % New point
    P=spm_normrnd(M.vpE,pC,1);
    
    % Parameters in reduced space
    Pr = M.V'*(P-M.vpE);
    
    M.int='sundials';
    tic;
    y = spm_mci_fwd (P,M,U);
    ysun(:,s)=y(:,1);
    els_sun(s)=toc;
    
    M.int='ode15';
    tic;
    y = spm_mci_fwd (P,M,U);
    yode(:,s)=y(:,1);
    els_ode(s)=toc;

    % DCM for fMRI uses spm_int
    % DCM for ERP uses spm_int_L
    tic;
    y=spm_int_L(P,M,UI);
    yspm(:,s)=y(:,1);
    els_spm(s)=toc;
    
end

hs=figure;
set(hs,'Name','First Time Series');
lw=2;
rN=ceil(sqrt(Nsamp));
for s=1:Nsamp,
    subplot(rN,rN,s);
    plot(ysun(:,s),'k');
    hold on
    grid on
    plot(yode(:,s),'b');
    plot(yspm(:,s),'r');
end
legend('Sundials','ODE15','spm-int-L');

hs=figure;
set(hs,'Name','Integration Speed');
k=1; lw=2;
plot(els_ode,els_sun,'kx','MarkerSize',10);
hold on
mo=min(els_ode);
ma=max(els_ode);
plot([mo ma],[mo ma],'k','LineWidth',lw);
set(gca,'FontSize',18);
grid on
xlabel('ODE15');
ylabel('Sundials');

hs=figure;
set(hs,'Name','Integration Speed');
k=1; lw=2;
plot(els_spm,els_sun,'kx','MarkerSize',10);
hold on
mo=min(els_spm);
ma=max(els_spm);
plot([mo ma],[mo ma],'k','LineWidth',lw);
set(gca,'FontSize',18);
grid on
xlabel('spm-int-L');
ylabel('Sundials');
