
% displays the results of a canonical variates analysis
% FORMAT spm_cva
%___________________________________________________________________________
%
% spm_cva prompts for the selection of CVA.mat and displays the 
% selected canonical image in the results window along with the
% appropriate canonical vector
%
%__________________________________________________________________________
% %W% %E%

%---------------------------------------------------------------------------
figure(2); clf; set(2,'Name','Canonical Images')
tmp = spm_get(1,'.mat','select MAN[C]OVA you wish to analyse','CVA');
CWD = strrep(tmp,'/CVA.mat',''); % Get directory name
i = spm_input(sprintf('Canonical image? %0.0f - %0.0f',1,size(CV,2)),1);


% display canonical images
%---------------------------------------------------------------------------
load([CWD,'/XYZ'])
load([CWD,'/CVA'])

figure(3); spm_clf
if V(3,1) == 1					% 2-dimensional data

	subplot(2,1,1)
	spm_mip(CV(:,i),XYZ,V(1:6));
	title(sprintf('canonical image %0.0f',i));
else
	axes('Position',[0.05 0.5 0.45 0.4])
	spm_mip(CV(:,i),XYZ,V(1:6));
	title(sprintf('canonical image %0.0f {+ve}',i));

	axes('Position',[0.55 0.5 0.45 0.4])
	spm_mip(-CV(:,i),XYZ,V(1:6));
	title(sprintf('canonical image %0.0f {-ve}',i));
end

% component scores
%---------------------------------------------------------------------------
subplot(2,2,3);
[x y] = bar(S(:,i));
fill(x,y,[1 1 1]*0.9);
xlabel 'condition'
ylabel 'component score'
title(sprintf('canonical value = %0.0f ',v(i)))
axis square
set(gca,'Xlim',[0 (size(S,1) + 1)]);


% textual information
%---------------------------------------------------------------------------
subplot(2,2,4); axis off; axis square
y  = 1;
x  = -0.2;
text(x,y,'Multivariate analysis','FontSize',16,'Fontweight','Bold');
y = y - 0.2;
text(x,y,pwd,'Fontweight','Bold');
y = y - 0.2;
text(x,y,sprintf('Chi-squared = %0.2f',chi));
y = y - 0.2;
text(x,y,sprintf('p value = %0.6f',pV),'FontSize',16);

