
% graphical display of adjusted data
% FORMAT spm_plot
%____________________________________________________________________________
%
% spm_plot is a CallBack script that uses variables in working memory to
% produce plots of adjusted activity at the location that is nearest to
% the coordinate in vector L = [x y z] {mm}
% These mean activities are the average estimates over subjects.  If these
% estimates derive from an activation study (the effects are factors
% described by the H partition) the data are plotted as a bar chart,  If
% they derive from covariates (C) then they are plotted as (i) functions of
% of the [linear compond specified by the contrast of the] covariate
% and (ii) scan number.
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

% condition effects or covariate ?
%----------------------------------------------------------------------------
d      = min(find(CONTRAST(con(1),:)));

% delete previous axis
%----------------------------------------------------------------------------
subplot(2,1,2); delete(gca)
subplot(2,1,2);
DIR    = CWD(max([1 max(find(CWD == '/'))]):length(CWD));

if d > (size(H,2) + size(K,2));	% covariate

	if size(B,2)		% more than one subject
		y     = [];
		for j = 1:size(B,2);
			y      = [y XA((B(:,j) > 0),i)];
		end
	else
		y      = XA(:,i);
	end
	v     = y;
	v(:)  = [K H C B G]*(CONTRAST(con(1),:)');
	subplot(2,2,3)
	plot(v,y,'.'); axis square
	line(v,y,'LineStyle','.','MarkerSize',12)
	xlabel(['compound of covariates'])
	ylabel('adjusted activity')
	d     = (max(v(:)) - min(v(:)))/4;
	set(gca,'XLim',[(min(v(:)) - d) (max(v(:)) + d)]);

	title(DIR,'FontSize',16,'FontWeight','Bold')

	subplot(2,2,4)
	v     = y;
	x     = 1:size(y,1);
	D     = [K H C B G];
	v(:)  = D*(pinv(D)*XA(:,i));
	plot(x,v,'b',x,y,':'); axis square
	line(x,y,'LineStyle','.','MarkerSize',12)
	xlabel('scan')
	ylabel('response and estimate')


else			% condition effects

	% organize adjusted data and get mean
	%-------------------------------------------------------------------
	y      = XA(:,i);
	x      = [];
	Y      = [];
	for j  = 1:length(y); x(j) = find(H(j,:));    end
	for j  = 1:size(H,2); Y(j) = mean(y(H(:,j))); end


	% plot
	%-------------------------------------------------------------------
	[u v]  = bar(Y);
	fill(u,v,[1 1 1]*.9)
	line(x,y,'LineStyle','.','MarkerSize',16)
	xlabel('level')
	ylabel('adjusted response')
	axis([0 (size(H,2) + 1) min(y)-1 max(y)+1]); axis square
	title(DIR,'FontSize',16,'FontWeight','Bold')
end
