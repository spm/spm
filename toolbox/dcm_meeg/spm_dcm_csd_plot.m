function spm_dcm_csd_plot(DCM,i1,i2,C)
% Plots the conditional density of coherence etc for a given connection
% FORMAT spm_dcm_csd_plot(DCM,i,j,C)
%
% DCM - inverted DCM structure for CSD models
% i   - target source (or channel mode)
% j   - source source (or channel mode)
% C   - flag for channels (as opposed to sources
%
% This routine is a graphics routine that plots the Bayesian confidence 
% tubes around cross-covariance, coherence and phase delays as functions 
% of lag and frequency. It also plots the conditional density over the 
% delay. The confidence tubes (Bayesian confidence intervals) are 
% approximated by sampling the underlying parameters from the 
% [approximate] conditional density.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_csd_plot.m 4095 2010-10-22 19:37:51Z karl $
 
 
% check options
%==========================================================================
try, C; catch, C = 0; end
    
 
% Compute conditional density from N samples (cf, numerical integration)
%--------------------------------------------------------------------------
N     = 32;
D     = exp(DCM.Ep.D(i1,i2))*16;
Cp    = spm_sqrtm(DCM.Cp);
Qp    = spm_vec(DCM.Ep);
Vp    = spm_unvec(diag(DCM.Cp),DCM.Ep);

for i = 1:N
    
    % predictions (csd) and error (sensor space)
    %----------------------------------------------------------------------
    qp  = Qp + Cp*randn(size(Qp));
    qp  = spm_unvec(qp,DCM.Ep);

    % if source space is required
    %----------------------------------------------------------------------
    if ~C
        DCM.M.dipfit.type = 'LFP';
        qp.L  = qp.L - qp.L;
        qp.b  = qp.b - 32;
        qp.c  = qp.c - 32;
    end
    hc        = spm_csd_mtf(qp,DCM.M,DCM.xU);
    
    % functions
    %----------------------------------------------------------------------
    [ccf pst] = spm_csd2ccf(hc,DCM.M.Hz);
    [coh fsd] = spm_csd2coh(hc,DCM.M.Hz);
    
    CCF(:,i)  = ccf{1}(:,i1,i2);
    COH(:,i)  = coh{1}(:,i1,i2);
    FSD(:,i)  = fsd{1}(:,i1,i2);
    
    fprintf('sample: %i(%i)\n',i,N);
end
 
 
% Plot
%==========================================================================
HZ  = DCM.Hz(:);
PST = DCM.pst(:)*1000;
 

% Covariance function
%-------------------------------------------------------------------------- 
subplot(2,2,1); hold off
i   = abs(PST) < 128;
CI  = 1.6449*std(CCF')';
CC  = mean(CCF')';
CI  = CI(i);
CC  = CC(i);
PST = PST(i);
 
fill([PST; flipud(PST)],[CC + CI; flipud(CC - CI)],[1 1 1]*.8), hold on
plot(PST,CC), hold on
plot([0 0],[min(CC - CI) max(CC + CI)],'-.'), hold on
plot([D D],[min(CC - CI) max(CC + CI)]),      hold on
title('Covariance','FontSize',16)
xlabel('lag (ms)')
ylabel('covariance')
axis square

% Coherence function
%-------------------------------------------------------------------------- 
subplot(2,2,2); hold off
CI  = 1.6449*std(COH')';
CC  = mean(COH')';
fill([HZ; flipud(HZ)],[CC + CI; flipud(CC - CI)],[1 1 1]*.8), hold on
plot(HZ,CC), hold on
title('Coherence','FontSize',16)
xlabel('a.u.')
ylabel('covariance')
axis square
 
% Phase delay function
%--------------------------------------------------------------------------
subplot(2,2,3); hold off
CI  = 1000*1.6449*std(FSD')';
CC  = 1000*mean(FSD')';
fill([HZ; flipud(HZ)],[CC + CI; flipud(CC - CI)],[1 1 1]*.8), hold on
plot(HZ,CC,[min(HZ) max(HZ)],[D D],'-.'), hold on
title('Phase delay','FontSize',16)
xlabel('delay (ms)')
ylabel('covariance')
axis square
 
% Conditional density over delay
%--------------------------------------------------------------------------
x  = linspace(1/128,32,128);
qp = spm_Npdf(log(x/16),DCM.Ep.D(i1,i2),Vp.D(i1,i2));
 
subplot(2,2,4),hold off
plot(x,qp,[D D],[0 max(qp)],'-.')
title('Conduction delay','FontSize',16)
xlabel('delay (ms)')
ylabel('conditional density')
axis square

