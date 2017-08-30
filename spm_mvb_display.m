function spm_mvb_display(MVB)
% Model display for MVB
% FORMAT spm_mvb_display(MVB)
% MVB  - multivariate Bayes structure, select one if not provided
%__________________________________________________________________________
% Copyright (C) 2007-2017 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mvb_display.m 7162 2017-08-30 11:47:07Z guillaume $
 

if ~nargin
    [filename, sts] = spm_select(1,'^MVB.*\.mat','Select MVB to display');
    if ~sts, return; end
    load(filename);
end
 
%-Get figure
%--------------------------------------------------------------------------
Fmvb  = spm_figure('GetWin','MVB');
spm_clf(Fmvb);
spm_figure('Focus',Fmvb);
 
%-Display specified model
%==========================================================================
M     = MVB.M;
XYZ   = MVB.XYZ;
VOX   = MVB.VOX;
X0    = MVB.X0;
X     = MVB.X;
Y     = MVB.Y;
 
%-Model comparison and selection, relative to the null model
%--------------------------------------------------------------------------
F     = M.F(2:end) - M.F(1);
P     = exp(F);
P     = P./(P + 1);
 
subplot(3,2,1)
bar(F), hold on
plot([0 length(F) + 1], [3 3],'r')
plot([0 length(F) + 1],-[3 3],'r')
plot([0 length(F) + 1], [5 5],'r:')
plot([0 length(F) + 1],-[5 5],'r:')
hold off
xlabel('partitions')
axis square
grid on
title({'log-evidence';sprintf('maximum p = %.2f%s',100*max(P),'%')})
 
subplot(3,2,2)
hist(M.qE,64)
xlabel('voxel-weight')
ylabel('frequency')
axis square
grid on
title({'distribution of weights'})
 
 
%-Posterior probabilities
%--------------------------------------------------------------------------
P   = 1 - spm_Ncdf(0,abs(M.qE),M.qC);
str{1,1} = 'Posterior probabilities at maxima  ';
str{2,1} = '________________________________';
str{3,1} = 'p(|w| > 0)    location (x,y,z)  weight (w)';
str{4,1} = '________________________________';
while length(str) < 16  && any(P)
    
    [p,i]  = max(P);
    str{end + 1,1} = sprintf('p = %.3f   %2.1f,%2.1f,%2.1fmm   q = %.4f;',...
                 p,XYZ(1,i),XYZ(2,i),XYZ(3,i),M.qE(i));
    P      = P.*(((XYZ(1,:) - XYZ(1,i)).^2 + ...
                  (XYZ(2,:) - XYZ(2,i)).^2 + ...
                  (XYZ(3,:) - XYZ(3,i)).^2) > 4^2)';
 
end
str{end + 1,1} = '________________________________';
str{end + 1,1} = sprintf('%i voxels; %i scans',length(P),length(X));
subplot(3,2,4)
text(0,1/2,str,'FontSize',10)
axis off
 
%-Maximum intensity projection
%--------------------------------------------------------------------------
subplot(3,2,3)
i  = P > .5;
spm_mip(P(i),XYZ(1:3,i),VOX)
axis image
title(['PPM: ' MVB.name ' (' MVB.contrast ')'])
 
%-Residual forming matrix
%--------------------------------------------------------------------------
R    = speye(size(X0,1)) - X0*pinv(X0);
Ns   = 1:size(X0,1);
 
%-Predictions and target (adjusted)
%--------------------------------------------------------------------------
X    = R*X;
P    = R*Y*M.qE;
SNR  = var(P)/var(X - P);
str  = sprintf('SNR (variance) %.2f',SNR);
 
subplot(3,2,5)
plot(Ns,X,'k',Ns,P,'r')
xlabel('scans')
ylabel('adjusted response')
axis square, grid on
try
    title({MVB.name [' (prior: ' M.priors ')']})
catch
    title(MVB.name)
end
legend({'target','prediction'})
 
subplot(3,2,6)
plot(X,P,'r.')
xlabel('contrast')
ylabel('prediction')
title({'observed and predicted contrast',str})
axis square, grid on

%-Print
%==========================================================================
spm_print('',Fmvb);
