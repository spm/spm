
% graphical display of adjusted data
% FORMAT spm_plot
%____________________________________________________________________________
%
% spm_graph is a CallBack script that uses variables in working memory to
% produce plots of adjusted activity at the significant (p < 0.05
% uncorrected according the the F statistic following the AnCova) that is
% nearest to the point selected.  These mean activities are the average
% estimates over subjects.  If these estimates derive from an activation
% study (the effects are factors described by the H partition) the data
% are plotted as a bar chart,  If they derive from covariates (C) then
% they are plotted as against [a compound of] the covariates or scan
% number as specified by you.  The adjusted data for the selected voxel
% are displayed in the command window for potential use outside SPM.  The
% vatiables x and y contain the ordinates and activities respectively in
% the order in which the scans were entered (i.e. the same as the design
% matrix, The variable Y represents the fitted responses [e.g. means]
% accross subjects).
%
%__________________________________________________________________________
% @(#)spm_plot.m	1.5 96/06/10


% find nearest voxel [in a Euclidean sense] in the point list of locations XYZ
%----------------------------------------------------------------------------
[d i] = min(sum(([(XYZ(1,:) - L(1));(XYZ(2,:) - L(2));(XYZ(3,:) - L(3))]).^2));
L     = XYZ(:,i);


% reset the pointer and posotion strings created by spm_sections_ui.m
%----------------------------------------------------------------------------
if V(3) == 1
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(hXstr,'String',sprintf('x = %0.0f',L(1)));
	set(hYstr,'String',sprintf('y = %0.0f',L(2)));
	set(X1,'Position',[L(1)  L(2) 1]);
else
	set(h1,'String',sprintf('%0.0f',L(1)));
	set(h2,'String',sprintf('%0.0f',L(2)));
	set(h3,'String',sprintf('%0.0f',L(3)));
	set(hXstr,'String',sprintf('x = %0.0f',L(1)));
	set(hYstr,'String',sprintf('y = %0.0f',L(2)));
	set(hZstr,'String',sprintf('z = %0.0f',L(3)));
	set(X1,'Position',[(P1 + L(2))  (P2 + L(1)) 1]);
	set(X2,'Position',[(P1 + L(2))  (P3 - L(3)) 1]);
	set(X3,'Position',[(P4 + L(1))  (P3 - L(3)) 1]);

end



% plot
%----------------------------------------------------------------------------
x    = [];
y    = [];
PLOT = 1;
if size(H,2)
	PLOT = spm_input('plot in terms of',1,'b','conditions|specify',[0 1]);
end


% delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
Finter = spm_figure('FindWin','Interactive');
figure(Fgraph)
subplot(2,1,2); delete(gca), spm_figure('DeletePageControls')
subplot(2,1,2); axis off

if PLOT

	% more than one subject with balanced design
	%--------------------------------------------------------------------
	if ~any(diff(sum(B > 0))) & size(B,2)		
		y     = [];
		for j = 1:size(B,2);
			y = [y XA((B(:,j) > 0),i)];
		end
	else
		y      = XA(:,i);
	end

	% get ordinate
	%--------------------------------------------------------------------
	str   = 'plot as a function of';
	if size(C,2)
		j    = spm_input(str,1,'b','cov|scan|specify',[2 1 3]);
	else
		j    = spm_input(str,1,'b','scan|specify',[1 3]);
	end

	if     j == 1
		x    = [1:size(y,1)]'*ones(1,size(y,2));
		Xlab = 'scan';

	elseif j == 2
		x    = y;
		str  = sprintf('which covariate 1 - %i',size(C,2));
		q    = spm_input(str,2,'e',1);
		x(:) = C(:,q);
		Xlab = sprintf('covariate %i',q);

	elseif j == 3
	    while any(size(x) - size(y))
		str  = sprintf('enter {%i x %i} ordinate',size(y,1),size(y,2));
		x    = spm_input(str,2);
		if any(size(x) == 1)
			x    = x(:)*ones(1,size(y,2));
		end
		if any(size(x) - size(y))
			x = x';
		end
	    end
	   str  = sprintf('ordinate name');
	   Xlab = spm_input(str,3,'s');

	end

	spm_clf(Finter);
	figure(Fgraph)
	subplot(2,1,2)
	Y     = y
	D     = [K H C B G];
	Y(:)  = D*(pinv(D)*XA(:,i));
	for j = 1:size(x,2)
		[p q] = sort(-x(:,j));
		plot(x(q,j),Y(q,j),'b',x(q,j),y(q,j),':'); hold on;
	end
	hold off; axis square
	line(x,y,'LineStyle','.','MarkerSize',12)
	xlabel(Xlab)
	ylabel('response and estimate')


% bar chart of condition effects
%---------------------------------------------------------------------------
else			
	% organize adjusted data and get mean
	%-------------------------------------------------------------------
	y      = XA(:,i)
	x      = [];
	Y      = [];
	for j  = 1:length(y); x(j) = find(H(j,:));    end
	for j  = 1:size(H,2); Y(j) = mean(y(H(:,j))); end


	% plot
	%-------------------------------------------------------------------
	[u v]  = bar(Y);
	fill(u,v,[1 1 1]*.9)
	line(x,y,'LineStyle','.','MarkerSize',12)
	xlabel('level')
	ylabel('adjusted response')
	axis([0 (size(H,2) + 1) (min(y) - 1) (max(y) + 1)])
	axis square
end

%----------------------------------------------------------------------------
if SPMZ

	Pu    = spm_P(1,W,t(i),0,S);		% voxel-level p value
	Pz    = 1 - spm_Ncdf(t(i));		% uncorrected p value (Z)
	TITLE = sprintf('Z = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

elseif SPMF

	Pu    = spm_pF(S,W,df,t(i));		% voxel-level p value
	Pz    = 1 - spm_Fcdf(t(i),df);		% uncorrected p value (F)
	TITLE = sprintf('F = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

end
title(TITLE,'FontSize',16)
spm_clf(Finter)
