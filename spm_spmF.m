function spm_spmF
% Display SPM{F} and design matrix
% FORMAT spm_spmF
%___________________________________________________________________________
%
% spm_spmF prompts for details of a SPM{F} that is then
% displayed in the Graphics window with the design matrix,
% the degrees of freedom and other parameters of the analysis 
%
% The SPM{F} is thresholded at a specified (uncorrected) threshold
%
%__________________________________________________________________________
% %W% %E%


% get SPM{F} and threshold
%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','SPM{F}')

tmp   = spm_get(1,'.mat','select SPMF.mat for analysis','SPMF');
U     = spm_input('threshold e.g. 2.8 or 0.001',1);

%---------------------------------------------------------------------------
set(2,'Pointer','watch')

CWD   = strrep(tmp,'/SPMF.mat',''); % Get directory name
K     = [];
names = [];
load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMF'])

if U < 1; U = spm_invFcdf(1 - U,Fdf); end

Q    = SPMF > U;


% display SPM{F}
%---------------------------------------------------------------------------
figure(3); spm_clf; subplot(2,1,1)
spm_mip(SPMF(Q),XYZ(:,Q),V(1:6))
xlabel(sprintf('SPM{F}: p < %0.3f {uncorrected}',1 - spm_Fcdf(U,Fdf)))
set(get(gca,'Xlabel'),'Visible','on')

% display design matrix and other details of the analysis
%---------------------------------------------------------------------------
axes('Position',[0.1 0.1 0.26 0.36]);
image((spm_DesMtxSca([K H C B G],names) + 1)*32)
title('Design matrix','Fontweight','Bold')
xlabel 'effect'; ylabel 'scan'

axes('Position',[0.46 0.1 0.42 0.3]); axis off
text(0,1.05,'Results directory:')
text(0,.95,CWD,'Fontweight','Bold')
text(0,.8,sprintf('Search volume %0.0f voxels',S))
text(0,.7,sprintf('Resolution {FWHM} %-4.1f %-4.1f %-4.1f mm',FWHM))
text(0,.6,sprintf('Image size  %-4.0f %-4.0f %-4.0f voxels',V(1:3)))
text(0,.5,sprintf('Voxel size  %-4.1f %-4.1f %-4.1f mm',V(4:6)))
text(0,.35,'Degrees of freedom','Fontsize',10)
text(0,.25,sprintf('Effects    %0.0f ',Fdf(1)),'Fontsize',10)
text(0,.20,sprintf('Confounds  %0.0f ',rank([K B G])),'Fontsize',10)
text(0,.15,sprintf('Residuals  %0.0f ',Fdf(2)),'Fontsize',10)
line([0 1],[1.15 1.15],'LineWidth',3);
line([0 1],[.86 .86],'LineWidth',3);
line([0 1],[0.3 0.3],'LineWidth',1);
line([0 1],[.06 .06],'LineWidth',1);

%---------------------------------------------------------------------------
set(2,'Pointer','Arrow','Name',''); figure(2); clf











