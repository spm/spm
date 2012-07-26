function DFP_demo_double_well
% DEMO comparing Variational filtering with particle filtering in the 
% context of a bimodal conditional density.  This demonstrates that the
% variational filter can not only represent free-form densities on the 
% states but also the causes of responses.
 
% phase diagram - to show bimodal energy function
%==========================================================================
spm_figure('GetWin','Figure 1');
 
x       = -32:1/16:32;
dxdt    = -x/2 + 16*x./(1 + x.^2);
V       =  1/4*x.^2 - 8*log(1 + x.^2);
 
subplot(2,1,1)                                 
plot(x,x*0,':',x,dxdt/8)
axis square
xlabel('state','FontSize',12)
ylabel('velocity','FontSize',12)
title('phase diagram','FontSize',14)
 
subplot(2,1,2)
plot(x,V/8)
axis square
xlabel('state','FontSize',12)
ylabel('potential','FontSize',12)
title('double well','FontSize',14)
 
 
% get nonlinear state-space model
%==========================================================================
spm_figure('GetWin','DEM');
M        = spm_DEM_M('ssm');
 
 
% generate data (output) and graph
%--------------------------------------------------------------------------
M(1).E.N = 32;                                        % number of particles
T        = 32;                                        % number of time bins
U        = 8*sin(pi*[1:T]/16);
DEM      = spm_DEM_generate(M,U);
 
spm_DEM_qU(DEM.pU);
 
% DFP
%--------------------------------------------------------------------------
DFP           = spm_DFP(DEM);
 
% PF
%--------------------------------------------------------------------------
[pf_x,P,Q,xQ] = spm_pf(M,DEM.Y);
 
% DEM
%--------------------------------------------------------------------------
DEM           = spm_DEM(DEM);
 
% Graphical comparison of DEM and true states
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
% Graphical comparison of DFP and true states
%--------------------------------------------------------------------------
spm_figure('GetWin','DFP');
spm_DFP_plot(DFP.QU,DFP.pU)
 
 
% Graphical comparison of DFP and PF
%==========================================================================
spm_figure('GetWin','Figure 2');
 
% show density on state (PF)
%--------------------------------------------------------------------------
subplot(2,2,2)
imagesc([1:T],xQ,Q)
xlabel('time (bins)')
ylabel('hidden state')
axis xy square
title('sample density (PF)')
 
% density on state (VF)
%--------------------------------------------------------------------------
for i = 1:T
    x = {DFP.QU{i}.x};
    for j = 1:length(x)
       qx(j) = x{j}{1};
    end
    q       = hist(qx,xQ);
    Qx(:,i) = q(:);
end
 
subplot(2,2,1)
imagesc([1:T],xQ,Qx)
xlabel('time (bins)')
ylabel('hidden state')
axis xy square
title('sample density (VF)')
 
 
% show density on cause
%==========================================================================
for i = 1:T
    x = {DFP.QU{i}.v};
    for j = 1:length(x)
       qv(j) = x{j}{1};
    end
    q       = hist(qv,xQ);
    Qv(:,i) = q(:);
end
 
spm_DFP_plot(DFP.QU,DFP.pU)
subplot(2,1,2)
imagesc([1:T],xQ,Qv)
xlabel('time (bins)')
ylabel('hidden state')
axis xy square
title('sample density (VF)')
