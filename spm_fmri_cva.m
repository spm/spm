function spm_fmri_cva
% User interface for spm_fmri_mancova.m - canonical variates analysis
% FORMAT spm_fmri_cva
%___________________________________________________________________________
%
% spm_fmri_cva prompts for the selection of CVA.mat and then
% displays canonical images and variates.  These canonical images
% are the spatial modes that reflect condition effects relative to error.
% The canonical variates represent the expression of these spatial modes over
% scans.  Due to the nature of the model used in the ManCova the canonical
% variates can be decomposed into a set of condition-dependent transients
% that are an estimate of the hemodynamic response to each condition (these
% are called canonical transients).
%
% see spm_fmri_mancova
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','Canonical Images, Variates & Transients')
tmp    = spm_get(1,'.mat','select CVA.mat you wish to analyse','CVA');
CWD    = strrep(tmp,'/CVA.mat','');
cd(CWD(CWD ~= ' '))
load CVA

e = spm_input(sprintf('Canonical image? %0.0f - %0.0f',1,size(CU,2)),1);


%---------------------------------------------------------------------------
set(2,'Pointer','Watch')

% display canonical spatial modes
%---------------------------------------------------------------------------
figure(3); spm_clf

if V(3) == 1						% 2-dimensional data
	subplot(2,1,1)
	spm_mip(CU(:,e),XYZ,V(1:6));
	title(sprintf('Canonical image %0.0f',e));
else
	axes('Position',[0.05 0.5 0.45 0.4])
	spm_mip(CU(:,e),XYZ,V(1:6));
	title(sprintf('Canonical image %0.0f {+ve}',e));

	axes('Position',[0.50 0.5 0.45 0.4])
	spm_mip(-CU(:,e),XYZ,V(1:6));
	title(sprintf('Canonical image %0.0f {-ve}',e));
end

axes('Position',[0 0.5 1 1]); axis off
text(0.1,0.02,'Canonical image analysis:','Fontsize',16,'FontWeight','Bold')
text(0.1,0,sprintf('P{dimension >= %0.0i} = %0.4f',e,pV(e)),'FontWeight','Bold')
text(0.52,0.02,CWD)


% canonical variate
%---------------------------------------------------------------------------
subplot(2,2,3)
d     = 1:size(CV,1);
plot(d,CV(:,e),d,X*E(:,e),'.');
set(gca,'XLim',([0 length(CV)] + 1))
axis square
xlabel('scan')
ylabel('relative expression')
title('Canonical Variate')



% canonical transients
%---------------------------------------------------------------------------
CF    = [];
d     = size(W,2);
for i = 1:(size(H,2)/d)
	CF = [CF W*BETA(([1:d] + d*(i - 1)),:)*E(:,e)]; end

subplot(2,2,4)
plot(CF)
set(gca,'Xlim',[1 size(CF,1)])
axis square
xlabel('scan')
ylabel('relative expression')
title('Canonical Transients')
grid on


%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','','Pointer','Arrow'); figure(3)
