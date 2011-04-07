% DEM demo for linear deconvolution:  This demo considers the deconvolution
% of the responses of a single-input-multiple output input-state-output
% model (DCM) to disclose the input or causes.  It focuses on estimating the
% causes and hidden states: The notes provide a comparative evaluation with 
% extended Kalman filtering.
 
 
% get a simple convolution model
%==========================================================================
f       = spm_figure('GetWin','Graphics');
M       = spm_DEM_M('convolution model');

M(1).V  = exp(8);                             % error precision
M(1).W  = exp(8);                             % error precision
 
% and generate data
%==========================================================================
N       = 32;                                 % length of data sequence
U       = exp(-([1:N] - 12).^2/(2.^2));       % this is the Gaussian cause;
DEM     = spm_DEM_generate(M,U,{},{[] 16},{16});
 
% display
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.pU)
 

% invert model
%==========================================================================
DEM     = spm_DEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
return
 
 
% explore dimensions (n)
%==========================================================================
for i = 1:12
    clear functions
    DEM.M(1).E.n = i;
    DEM.M(1).E.d = 3;
    
    D{i}   = spm_DEM(DEM);
    F(i)   = D{i}.F;
    Sx(i)  = sum(sum((D{i}.qU.x{1} - DEM.pU.x{1}).^2));
    Sv(i)  = sum(sum((D{i}.qU.v{2} - DEM.pU.v{2}).^2));
end
 
% plot
%--------------------------------------------------------------------------
figure(f)
subplot(2,1,1)
bar(Sv)
xlabel('n - 1 (d = 2)')
ylabel('sum squared error (causal states)')
axis square
 
d     = [1 8];
for i = 1:length(d)
    subplot(2,length(d),i + length(d))
    hold on
    plot([1:N],D{d(i)}.qU.v{2})
    plot([1:N],DEM.pU.v{2},'linewidth',2,'color',[1 1 1]/2)
    hold off
    axis square
    xlabel('time')
    title(sprintf('n = %i',d(i)))
    if i == 1, a = axis; else, axis(a); end
end
 
% and d
%--------------------------------------------------------------------------
clf; clear D F Sx Sv
for i = 1:6
    clear functions
    DEM.M(1).E.n = 15;
    DEM.M(1).E.d = i - 1;
    
    D{i}   = spm_DEM(DEM);
    F(i)   = D{i}.F;
    Sx(i)  = sum(sum((D{i}.qU.x{1} - DEM.pU.x{1}).^2));
    Sv(i)  = sum(sum((D{i}.qU.v{2} - DEM.pU.v{2}).^2));
    
end
 
% plot
%--------------------------------------------------------------------------
clf
subplot(2,2,1)
bar(Sx)
xlabel('d (n = 15)')
ylabel('sum squared error (hidden states)')
axis square
set(gca,'XTickLabel',[0:7])

subplot(2,2,2)
bar(F - min(F) + 32)
xlabel('d (n = 15)')
ylabel('log-evidence')
axis square
set(gca,'XTickLabel',[0:7])
 
d     = [1 4];
for i = 1:length(d)
    subplot(2,length(d),i + length(d))
    hold on
    plot([1:N],D{d(i)}.qU.x{1})
    plot([1:N],DEM.pU.x{1},'linewidth',2,'color',[1 1 1]/2)
    hold off
    axis square
    xlabel('time')
    title(sprintf('d = %i',d(i) - 1))
    if i == 1, a = axis; else, axis(a); end
end
 
 
% Comparison with EKF
%==========================================================================
clear SSE
for i = 1:8
 
    % i.i.d.
    %----------------------------------------------------------------------
    DEM      = spm_DEM_generate(M,U,{},{[] 16});
    DEM      = spm_DEM(DEM);
 
    % serial correlations
    %----------------------------------------------------------------------
    D0          = DEM;
    D0.M(1).E.s = 0;
    D0          = spm_DEM(D0);
 
    % extended kalman filter - i.i.d.
    %----------------------------------------------------------------------
    e_x      = spm_ekf(DEM.M,DEM.Y);
    d_x      = DEM.qU.x{1};
    t_x      = DEM.pU.x{1};
    i_x      = D0.qU.x{1};
 
    SSE(i,1) = sum(sum((e_x - t_x).^2));
    SSE(i,2) = sum(sum((i_x - t_x).^2));
    SSE(i,3) = sum(sum((d_x - t_x).^2));
 
end
 
 
subplot(2,1,1)
plot(1:N,i_x,'g',1:N,d_x,'r',1:N,e_x,'b',1:N,t_x,'k')
legend('DEM(0)',' ','DEM(0.5)',' ','EKF',' ','true',' ')
xlabel('time')
title('hidden states')
axis square
 
subplot(2,1,2)
plot(SSE','k:'), hold on
plot(SSE','k.','Markersize',16), hold off
set(gca,'Xtick',[1 2 3],'XLim',[0 4])
set(gca,'Xticklabel',{'EKF','DEM(0)','DEM(0.5)'})
title('sum of squared error (hidden states)')
axis square
 
 
% Show equivalence when causes are not structured
%==========================================================================
M(1).E.s  = 0;
M(1).pE.h = eye(2);
M(1).pC   = [];

M(2).v    = [0;0];
M(2).V    = eye(2);
U         = [U; 0*U];

DEM       = spm_DEM_generate(M,U,{},{[] 16});
DEM       = spm_DEM(DEM);
e_x       = spm_ekf(DEM.M,DEM.Y);
d_x       = DEM.qU.x{1};
    
clf
subplot(2,2,1)
plot(1:N,d_x,'r',1:N,e_x,'b')
legend('DEM(0)',' ','EKF',' ')
xlabel('time')
title('hidden states')
axis square

% Show equivalence when causes are removed
%==========================================================================
DEM.M(1).pE.h = sparse(2,2);
DEM.M(1).W    = speye(2);

DEM           = spm_DEM(DEM);
e_x           = spm_ekf(DEM.M,DEM.Y);
d_x           = DEM.qU.x{1};
    
clf
subplot(2,2,1)
plot(1:N,d_x,'r',1:N,e_x,'b')
legend('DEM(0)',' ','EKF',' ')
xlabel('time')
title('hidden states')
axis square
