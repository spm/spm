function [CVA] = spm_cva_results(CVA)
% displays the results of a CVA analysis
% FORMAT [CVA] = spm_cva_results(CVA)
%
% xSPM   - structure containing specific SPM details
% SPM    - structure containing generic analysis details
% hReg   - Handle of results section XYZ registry (see spm_results_ui.m)
%
% CVA.contrast =  contrast name
% CVA.name     =  CVA name
% CVA.c        =  contrast weights
% CVA.X        =  contrast subspace
% CVA.Y        =  whitened and adjusted data
% CVA.X0       =  null space of contrast
% 
% CVA.XYZ      =  locations of voxels (mm)
% CVA.xyz      =  seed voxel location (mm)
% 
% CVA.V        =  canonical vectors  (data)
% CVA.v        =  canonical variates (data)
% CVA.W        =  canonical vectors  (design)
% CVA.w        =  canonical variates (design)
% CVA.C        =  canonical contrast (design)
% 
% CVA.chi      =  Chi-squared statistics testing D >= i
% CVA.df       =  d.f.
% CVA.p        =  p-values
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_cva_results.m 3717 2010-02-08 16:44:42Z guillaume $
 

% get CVA if necessary
%--------------------------------------------------------------------------
if nargin<1
    [file,sts] = spm_select(1,'mat',...
        'Select CVA to display',[],[],'^CVA.*\.mat$');
    if ~sts, CVA = []; return; end
    file = load(file);
    CVA  = file.CVA;
end

% show results
%==========================================================================
spm_figure('GetWin','MVB');

% unpack
%--------------------------------------------------------------------------
VOX      = CVA.VOX;
XYZ      = CVA.XYZ;
 
% maximum intensity projection (first canonical image)
%--------------------------------------------------------------------------
subplot(2,2,1)
spm_mip(CVA.V(:,1).*(CVA.V(:,1) > 0),XYZ(1:3,:),diag(VOX))
axis image
title({'(Principal) canonical image',[CVA.name ':' CVA.contrast]})
 
% inference and canonical variates
%--------------------------------------------------------------------------
Xstr{1} = 'Dimensionality';
Xstr{2} = ['Chi-squared: ' sprintf('%6.1f ',  CVA.chi)];
Xstr{3} = ['           df: ' sprintf('%6.0f ',CVA.df) ];

subplot(2,2,2)
bar(log(CVA.p)); hold on
plot([0 (length(CVA.p) + 1)],log(0.05)*[1 1],'r:','LineWidth',4), hold off
xlabel(Xstr)
ylabel('log p-value')
axis square
title({'Test of dimensionality';sprintf('minimum p = %.2e',min(CVA.p))})
 
subplot(2,2,3)
plot(CVA.w,CVA.v,'.')
xlabel('prediction')
ylabel('response')
axis square
title('Canonical variates')


% canonical contrast
%--------------------------------------------------------------------------
i     = find(CVA.p < 0.05);
str   = 'Significant canonical contrasts';
if isempty(i)
    i   = 1;
    str = 'first canonical contrast';
end
subplot(2,2,4)
bar(CVA.C(:,i))
xlabel('Parameter')
axis square
title(str)
