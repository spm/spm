
% graphical display of adjusted data using paramter estimates
% FORMAT spm_spmF_plot
%____________________________________________________________________________
%
% spm_plot is a CallBack script that uses variables in working memory to
% produce plots of adjusted activity at the location that is nearest to
% the coordinate in vector L = [x y z] {mm}
%
%__________________________________________________________________________
% %W% %E%


% find nearest voxel [in a Euclidean sense] in the point list of locations XYZ
%----------------------------------------------------------------------------
[d i] = min(sum(([(XYZ(1,:) - L(1));(XYZ(2,:) - L(2));(XYZ(3,:) - L(3))]).^2));
L     = XYZ(:,i);


% reset the pointer and posotion strings created by spm_sections_ui.m
%----------------------------------------------------------------------------
if V(3) == 1
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(X1,'Position',[L(1)  L(2) 1]);
else
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(h3,'String',sprintf('%0.0f',L(3)));
	set(X1,'Position',[(124 + L(2))  (248 + L(1)) 1]);
	set(X2,'Position',[(124 + L(2))  (112 - L(3)) 1]);
	set(X3,'Position',[(276 + L(1))  (112 - L(3)) 1]);
end

% get voxel specific data y, Z score and P(Zmax > Z) from the SPM{Z} 
%----------------------------------------------------------------------------
y      = XA(:,i);
y      = y - mean(y);
p      = [H C]*BETA([1:size([H C],2)],i);

% if hold is 'on' compare current and previous regressions
%----------------------------------------------------------------------------
subplot(2,1,2);
if strcmp(get(gca,'NextPlot'),'add') & exist('y0')
	DG  = [  [H C];  [H C]];
	DH  = [0*[H C]; -[H C]];
 	dfh = rank(DH);
	dfr = (size(DG,1) - rank([DH DG]));
	R   = [y; y0] - [DH DG]*pinv([DH DG])*[y; y0];
	N   = [y; y0] - [   DG]*pinv([   DG])*[y; y0];
 	F   = (N'*N  - R'*R)/dfh./(R'*R/dfr);

	TITLE = sprintf('difference F = %0.2f df; %0.0f,%0.0f p = %0.3f',...
	F,dfh,dfr,(1 - spm_Fcdf(F,[dfh dfr])));
else
	TITLE = pwd;
end
y0     = y;


%----------------------------------------------------------------------------
for i = 0:size(B,2);		% for each block

	if ~i
		if size(B,2) > 0;
			d   = (B(:,1) == -1);
		else
			d   = ones(size([H C B G]),1);
		end
	else
		d   = (B(:,i) == 1);
	end
	x   = 1:sum(d);

	plot(x,p(d),'b'); hold on
	line(x,y(d),'LineStyle','.','MarkerSize',16)

end
axis square; hold off;
xlabel('observation')
ylabel('adjusted activity {mean corrected}')
title(TITLE,'FontSize',16,'FontWeight','Bold')

