
% graphical display of adjusted data
% FORMAT spm_graph
%_______________________________________________________________________
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
%_______________________________________________________________________
% %W% Karl Friston %E%


% Find nearest voxel [Euclidean distance] in the point list of locations XYZ
%-----------------------------------------------------------------------
[d i] = min( sum(([	(XYZ(1,:) - L(1));...
			(XYZ(2,:) - L(2));...
			(XYZ(3,:) - L(3))	]).^2));


% Reset the pointer and position strings created by spm_results_ui
%-----------------------------------------------------------------------
L     = XYZ(:,i);
spm_mip_ui('SetCoords',L);


% Get adjusted data y and fitted effects Y
%-----------------------------------------------------------------------
y     = spm_readXA(QQ(i));			% adjusted data
BETA  = pinv([H C 0*B 0*G])*y;			% parameter estimate
Y     = [H C 0*B 0*G]*BETA;			% fitted data
R     = y - Y;					% residuals
RES   = sum(R.^2);				% SSQ of residuals
SE    = sqrt(RES*diag(BCOV));			% standard error of estimates
HC    = [H C];
MSize = 8;
COL   = ['r' 'b' 'g' 'c' 'y'];

% Inference (for title)
%-----------------------------------------------------------------------
if SPMZ

	Pu    = spm_P(1,W,t(i),0,S);		% voxel-level p value
	Pz    = 1 - spm_Ncdf(t(i));		% uncorrected p value (Z)
	STR   = sprintf('Z = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

%-----------------------------------------------------------------------
elseif SPMF

	Pu    = spm_pF(S,W,df,t(i));		% voxel-level p value
	Pz    = 1 - spm_Fcdf(t(i),df);		% uncorrected p value (F)
	STR   = sprintf('F = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

end


% Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
Fgraph = spm_figure('FindWin','Graphics');
Finter = spm_figure('FindWin','Interactive');
figure(Fgraph)
subplot(2,1,2); spm_figure('DeletePageControls')
if ~strcmp(get(gca,'NextPlot'),'add'); 
	delete(gca); subplot(2,1,2); axis off;
end


% modeling evoked responses
%----------------------------------------------------------------------
if exist('ERI')
    if spm_input('plot event-related responses',1,'b','yes|no',[1 0]);
	j     = 1;
	if size(ERI,2) > 1
		j = spm_input('which events','!+1','e','1 2');
	end

	Cplot = str2mat(...
			'Fitted response',...
			'Fitted response +/- standard error of response.',...
			'Fitted response +/- standard error of onset.',...
			'Fitted response and adjusted data');
	str   = 'plot in terms of';
	Cp    = spm_input(str,'!+1','m',Cplot,[1:size(Cplot,1)]);
	TITLE = deblank(Cplot(Cp,:));


	% cycle over selected events
	%--------------------------------------------------------------
	figure(Fgraph)
	hold on
	dx    = 0.1;
	x     = -4:dx:32;
	if length(x) ~= size(DER,1)
		DER = [zeros(4/dx,size(DER,2)), DER];
		x   = [0:(size(DER,1) - 1)]*dx - 4;
	end
	K     = spm_sptop(SIGMA*RT/dx,length(x));
	KDER  = K*DER;
	for i = j
		Y      = KDER*BETA(ERI(:,i));
		se     = sqrt(diag(KDER*BCOV(ERI(:,i),ERI(:,i))*KDER')*RES);
		pst    = PST(:,i);
		d      = round(pst + 4)/dx;
		q      = find( (d > 0) & (d <= size(KDER,1)) );
		y      = KDER(d(q),:)*BETA(ERI(:,i)) + R(q);
		d      = min(find(abs(Y) > max(abs(Y))/3));
		T      = x(d);
		dYdt   = gradient(Y')'/dx;
		seT    = se(d)./dYdt(d);

		% plot
		%------------------------------------------------------
		if Cp == 1
			plot(x,Y,COL(i))

		elseif Cp == 2
			plot(x,Y,x,Y + se,'-.',x,Y - se,'-.','Color',COL(i))

		elseif Cp == 3
			plot(x,Y,'Color',COL(i))
			line(([0-seT seT] + T),[Y(d) Y(d)],...
				'LineWidth',6,'Color',COL(i))

		elseif Cp == 4
			plot(x,Y,pst,y,'.',...
				'MarkerSize',MSize,'Color',COL(i))

		end
	end

	hold off; axis on
	set(gca,'XLim',[min(x) max(x)])
	xlabel('peri-stimulus time {secs}')
	ylabel(STR,'FontSize',8)
	title(TITLE,'FontSize',16)

	spm_graph_ui(gca)
	return

    end
end


% find out what to plot
%----------------------------------------------------------------------
Cplot = str2mat(...
		'Parameter estimates',...
		'Parameter estimates +/- standard error',...
		'Fitted response',...
		'Adjusted data',...
		'Fitted response and adjusted data');
str   = 'plot in terms of';
Cp    = spm_input(str,1,'m',Cplot,[1:size(Cplot,1)]);
TITLE = deblank(Cplot(Cp,:));

% plot parameter estimates
%----------------------------------------------------------------------
if Cp == 1 | Cp == 2
	i = 1:size(HC,2);
	if length(i) > 1
		if spm_input('which effects','!+1','b','all|specify',[0 1]);
			i = spm_input('which effects {columns}','!+1','e',1);
		end
	end
	BETA = BETA(i);
	SE   = SE(i);

	% bar chart
	%--------------------------------------------------------------
	figure(Fgraph)
	[p q]  = bar(BETA);
	fill(p,q,[1 1 1]*.9)
	if Cp == 2
	  for j = 1:length(BETA)
	    line([j j],([SE(j) 0 - SE(j)] + BETA(j)),'LineWidth',6,'Color','r')
	  end
	end
	xlabel('effect')
	ylabel(STR,'FontSize',8)
	title(TITLE,'FontSize',16)

	spm_graph_ui(gca)
	return



% All fitted effects or selected effects
%-----------------------------------------------------------------------
elseif Cp >= 3

	% bar chart of condition means if H exists
	%---------------------------------------------------------------
	if size(H,2)
	
		% organize adjusted data and get mean
		%-------------------------------------------------------
		i      = any(H')';
		y      = y(i);
		h      = H(i,:);
		x      = [];
		Y      = [];
		for j  = 1:length(y)
			x(j) = find(h(j,:) > 0);
		end
		for j  = 1:size(h,2)
			Y(j) = mean(y(h(:,j)));
		end

		% bar chart
		%-------------------------------------------------------
		figure(Fgraph)
		if Cp == 3

			[p q]  = bar(Y);
			fill(p,q,[1 1 1]*.9)

		elseif Cp == 4

			plot(x,y,'LineStyle','.','MarkerSize',MSize)

		elseif Cp == 5

			[p q]  = bar(Y);
			fill(p,q,[1 1 1]*.9)
			line(x,y,'LineStyle','.','MarkerSize',MSize)

		end

		xlabel('effect')
		ylabel(STR,'FontSize',8)
		title(TITLE,'FontSize',16)

		spm_graph_ui(gca)
		return

	end


	% get ordinates
	%---------------------------------------------------------------
	Cplot = str2mat(...
			'A covariate',...
			'scan or time.',...
			'user specified ordinate');
	str   = 'plot as a function of';
	Cx    = spm_input(str,'!+1','m',Cplot,[1:size(Cplot,1)]);

	if Cx == 1

		str  = sprintf('covariate 1 - %i',size(C,2));
		q    = spm_input(str,'!+1','e',1);
		x    = C(:,q);
		Xlab = sprintf('covariate %i',q);

	elseif Cx == 2

		if length(RT)
			x    = RT*[1:size(Y,1)]';
			Xlab = 'time {seconds}';
		else
			x    = [1:size(Y,1)]';
			Xlab = 'scan';

		end

	elseif Cx == 3

		x    = [];
		str  = sprintf('enter {1 x %i} ordinate',size(Y,1));
		while length(x) ~= length(Y)
			x = spm_input(str,'!+1');
			x = x(:);
		end
		Xlab = 'ordinate';

	end

	% plot
	%---------------------------------------------------------------
	figure(Fgraph)
	[p q] = sort(x);

	if Cp == 3

		plot(x(q),Y(q))

	elseif Cp == 4

		plot(x(q),y(q))

	elseif Cp == 5

		if all(diff(x(q)))
			plot(x(q),Y(q),'b',x(q),y(q),':')
		else
			plot(x(q),Y(q),'b')
		end
		axis square
		line(x(q),y(q),'LineStyle','.','MarkerSize',MSize)

	end
	xlabel(Xlab)
	ylabel(STR,'FontSize',8)
	title(TITLE,'FontSize',16)

	spm_graph_ui(gca)
	return

end


