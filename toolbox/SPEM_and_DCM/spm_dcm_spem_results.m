function DCM = spm_dcm_spem_results(DCM)
% Display (DCM results of) of smooth pursuit eye movements
% FORMAT DCM = spm_dcm_spem_results(DCM)
%
% DCM
%     name: name string
%       xY: data   [1x1 struct]
%       xU: design [1x1 struct]
%       pE: prior expectation
%       pC: prior covariance
%
% and (if inverted)
%
%       Y{i}   - predicted responses
%       DEM{i} - ADEM inversion structure
%       Ep     - posterior expectation
%       Cp     - posterior covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_spem_results.m 6014 2014-05-23 15:00:35Z guillaume $
 
% model specification
%--------------------------------------------------------------------------
 
% preliminaries
%--------------------------------------------------------------------------
N    = length(DCM.xY.y{1});
if ~isfield(DCM.xY,'u')
    DCM.xY.u{1} = cos(2*pi/N*(0:N - 1)');
end
if ~isfield(DCM.xY,'dt')
    DCM.xY.dt = 8;
end
if ~isfield(DCM.xY,'occ')
    DCM.xY.occ = @(x) x - x + 1;
end
 
% occlusion times
%--------------------------------------------------------------------------
u    = DCM.xY.u{1};
o    = find(diff(DCM.xY.occ(u))) + 1;
t    = (1:N)*DCM.xY.dt;
y    = DCM.xY.y;
col  = {'r','b','g','c','m'};
 
% empirical responses
%==========================================================================
subplot(3,2,1)
plot(t,u,'-.k'), hold on
str   = {'target';};
for i = 1:length(y)
    plot(t,y{i} + u,col{i}), hold on
    str{i + 1} = sprintf('condition %i',i);
end
for i = 1:length(o)
    plot(t(o(i))*[1 1],[-1 1],':')
end
xlabel('Time (sec)')
ylabel('Displacement')
title('Empirical SPEM','FontSize',16)
hold off
spm_axis square
spm_axis tight
 
subplot(3,2,2)
plot(t,t - t,'-.k'), hold on
for i = 1:length(y)
    plot(t,y{i},col{i})
end
for i = 1:length(o)
    plot(t(o(i))*[1 1],[-1 1]/4,':')
end
xlabel('Time (sec)')
ylabel('Lag')
title('Empirical lag','FontSize',16)
hold off 
spm_axis square
spm_axis tight
legend(str)
 
% Predicted responses
%==========================================================================
if ~isfield(DCM,'Y')
    return
else
    y    = DCM.Y;
    u    = DCM.DEM{1}.pU.v{2}(1,:)';
    u    = u/max(u);
end
 
% predicted responses
%--------------------------------------------------------------------------
subplot(3,2,3)
plot(t,u,'-.k'), hold on
for i = 1:length(y)
    plot(t,y{i} + u,col{i}), hold on
end
for i = 1:length(o)
    plot(t(o(i))*[1 1],[-1 1],':')
end
xlabel('Time (sec)')
ylabel('Displacement')
title('Predicted SPEM','FontSize',16)
hold off
spm_axis square
spm_axis tight
 
subplot(3,2,4)
plot(t,t - t,'-.k'), hold on
for i = 1:length(y)
    plot(t,y{i},col{i})
end
for i = 1:length(o)
    plot(t(o(i))*[1 1],[-1 1]/4,':')
end
xlabel('Time (sec)')
ylabel('Lag')
title('Predicted lag','FontSize',16)
hold off 
spm_axis square
spm_axis tight
 
% Parameter estimates
%==========================================================================
if ~isfield(DCM,'Ep')
    return
else
    Ep = DCM.Ep;
    Cp = spm_unvec(diag(DCM.Cp),Ep);
end
 
subplot(3,2,5)
spm_plot_ci(Ep.A,Cp.A)
xlabel('Parameter')
ylabel('Mean and 90% CI')
title('Parameter estimates','FontSize',16)
spm_axis square
 
subplot(3,2,6)
spm_plot_ci(Ep.B,Cp.B)
xlabel('Parameter')
ylabel('Mean and 90% CI')
title(sprintf('and (%i) Experimental effects',numel(Ep.B)),'FontSize',16)
spm_axis square
