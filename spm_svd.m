
% eigenimage display
% FORMAT spm_svd
% requires the following to be in working memory
% u    - eigenimages or spatial modes
% v    - component scores
% s    - eigenvalues
% e    - selected mode
% XYZ  - voxel locations in mm
%____________________________________________________________________________
%
% spm_svd displays the positive and negative parts of an eigenimage
% as maximum intensity projections.  The eigenvalue spectrum and component
% scores are also displayed.  spm_svd is essentially a script invoked by
% spm_svd_ui.m
%
%__________________________________________________________________________
% %W% Karl Friston %E%


%---------------------------------------------------------------------------
global CWD

% find and clear Interactive window
%---------------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph=spm_figure('Create','Graphics'); end
figure(Fgraph), spm_clf(Fgraph)


% eigenimages
%---------------------------------------------------------------------------
if V(3) == 1						% 2-dimensional data
	subplot(2,1,1)
	spm_mip(u(:,e),XYZ,V(1:6));
	title(sprintf('eigenimage %0.0f',e));
else
	axes('Position',[0.05 0.5 0.45 0.4])
	spm_mip(u(:,e),XYZ,V(1:6));
	title(sprintf('eigenimage %0.0f {+ve}',e));

	axes('Position',[0.50 0.5 0.45 0.4])
	spm_mip(-u(:,e),XYZ,V(1:6));
	title(sprintf('eigenimage %0.0f {-ve}',e));
end

axes('Position',[0 0.5 1 1]); axis off
text(0.1,0,'Eigenimage analysis:','Fontsize',16,'FontWeight','Bold')
text(0.42,0,CWD)


% eigenvalue spectrum
%---------------------------------------------------------------------------
subplot(2,2,3)
d     = zeros(size(s)); d(e) = s(e);
[x y] = bar(s);
fill(x,y,[1 1 1]*0.9); hold on
[x y] = bar(d);
fill(x,y,[1 0 0]); hold off
axis square; set(gca,'Xlim',[0 size(u,2)])
xlabel eigenimage
ylabel eigenvalue
title(sprintf('%0.1f percent of variance',s(e)/length(s)*100))

% component scores
%---------------------------------------------------------------------------
subplot(2,2,4)
[x y] = bar(v(:,e));
fill(x,y,[1 1 1]*0.9);
axis square; set(gca,'Xlim',[0 (size(u,2) + 1)])
xlabel 'observation'
ylabel 'component score'


