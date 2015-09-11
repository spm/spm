function [Fx,Fp,FxFD,FpFD] = mci_compare_jacobians (model)
% Compare user supplied and finite difference methods
% FORMAT [Fx,Fp,FxFD,FpFD] = mci_compare_jacobians (model)
%
% model     'phase'
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_jacobians.m 6548 2015-09-11 12:39:47Z will $

switch model,
    case 'phase',
        d=3;
        disp(sprintf('Weakly coupled oscillator network : d=%d regions',d));
        [P,M,u,Y] = mci_phase_init(d);
        M.dfdx='mci_phase_dfdx';
        M.dfdp='mci_phase_dfdp';
        
    otherwise
        disp('Unknown model');
        return
end

M.x0=randn(M.n,1);

disp('User supplied Jacobian (blue)');
disp('Finite difference Jacobian (red)');
disp(' ');

Fx = feval(M.dfdx,M.x0,u,P,M);

FxFD = spm_diff(M.f,M.x0,u,P,M,1);

figure
for i=1:M.n,
    subplot(M.n,1,i);
    plot(Fx(i,:));
    hold on
    plot(FxFD(i,:),'r');
    ylabel(sprintf('df(%d)/dx(j)',i));
    grid on
end
xlabel('j');
disp('Percentage discrepancy in dfdx:');
100*sum(abs(Fx(:)-FxFD(:)))/sum(abs(FxFD(:)))

Fp = feval(M.dfdp,M.x0,u,P,M);
FpFD = spm_diff(M.f,M.x0,u,P,M,3);

figure
for i=1:M.n,
    subplot(M.n,1,i);
    plot(Fp(i,:));
    hold on
    plot(FpFD(i,:),'r');
    ylabel(sprintf('df(%d)/dp(j)',i));
    grid on
end
xlabel('j');
disp('Percentage discrepancy in dfdp:');
100*sum(abs(Fp(:)-FpFD(:)))/sum(abs(FpFD(:)))
