
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

if ~exist('PST'); PST = []; end

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
y     = spm_readXA(QQ(i));
Y     = [K H C B G]*(pinv([K H C B G])*y);
h     = H;
c     = C;


% Inference (for title)
%-----------------------------------------------------------------------
if SPMZ

	Pu    = spm_P(1,W,t(i),0,S);		% voxel-level p value
	Pz    = 1 - spm_Ncdf(t(i));		% uncorrected p value (Z)
	TITLE = sprintf('Z = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

%-----------------------------------------------------------------------
elseif SPMF

	Pu    = spm_pF(S,W,df,t(i));		% voxel-level p value
	Pz    = 1 - spm_Fcdf(t(i),df);		% uncorrected p value (F)
	TITLE = sprintf('F = %0.2f, p = %0.3f (%.3f corrected)',t(i),Pz,Pu);

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

% More than one subject or session
%----------------------------------------------------------------------
b     = size(B,2);
if  b > 1

	BLOCK = spm_input('which sessions/subjects',...
		1,'b','all|mean|specify',[0 1 2]);

	% mean over block
	%--------------------------------------------------------------
	if BLOCK == 1
		d     = 0;
		D     = 0;
		if size(C,2)
			c = 0;
		end
		if size(H,2)
			h = 0;
		end
		for j = 1:size(B,2);
			d = d + y((B(:,j) > 0))/b;
			D = D + Y((B(:,j) > 0))/b;
			if size(C,2)
				c = c + C((B(:,j) > 0),:)/b;
			end
			if size(H,2)
				h = h + H((B(:,j) > 0),:)/b;
			end
		end
		y     = d;
		Y     = D;

	% subset of block
	%--------------------------------------------------------------
	elseif BLOCK == 2

		j     = spm_input('sessions/subjects eg 1:3','!+1');
		j     = B(:,j);
		j     = any(j' > 0);
		y     = y(j(:));
		Y     = Y(j(:));
		if size(C,2)
			c = C(j(:),:);
		end
		if size(H,2)
			h = H(j(:),:);
		end
	end
end

% Bar chart of condition effects
%-----------------------------------------------------------------------
BAR    = 0;
if size(H,2)
	BAR = spm_input('condition bar chart',1,'b','yes|no',[1 0]);
end


% All fitted effects or selected effects
%-----------------------------------------------------------------------
if spm_input('which effects or events',1,'b','all|specify',[0 1])

	if size(PST,2)
		str = sprintf('which epoch or event (1-%0.0f)',size(PST,2));
		e   = spm_input(str,1,'e',1);

		% columns of design matrix (C)
		%-------------------------------------------------------
		d   = size(c,2)/size(PST,2);
		j   = [1:d] + (e - 1)*d;
	else
		% get effects
		%-------------------------------------------------------
		str = sprintf('which effects (1-%0.0f)',size([h c],2));
		j   = spm_input(str,1,'e',1);
	end

	% re-fit effects
	%---------------------------------------------------------------
	g       = [h c];
	c       = g(:,j);
	g(:,j)  = [];
	ey      = mean(y);
	y       = y - ey;
	BETA    = pinv([c g])*y;
	y       = y - [0*c g]*BETA + ey;
	Y       = [c 0*g]*BETA + ey;

	% re-fit effects
	%---------------------------------------------------------------
	if BAR
		h = h(:,j);
	end
end


% epoch/event-related responses
%----------------------------------------------------------------------
if size(PST,2) & ~BAR
	ER   = spm_input('plot event-related responses',1,'b','yes|no',[1 0]);
end

% get ordinates
%----------------------------------------------------------------------
MSize = 12;
if ER

	str   = sprintf('plot in relation to event (1-%0.0f)',size(PST,2));
	e     = spm_input(str,1,'e',1);

	% get ordinates from PST
	%--------------------------------------------------------------
	x     = PST(:,e);
	Xlab  = 'peri-stimulus time {secs}';
	MSize = 4;

elseif ~BAR

	% get ordinate
	%--------------------------------------------------------------
	COV  = 0;
	if size(C,2)
		str  = 'plot as a function of';
		COV  = spm_input(str,'!+1','b','covariate|scan',[1 0]);
	end
	if COV
		str  = sprintf('covariate 1 - %i',size(C,2));
		q    = spm_input(str,'!+1','e',1);
		x    = C(:,q);
		Xlab = sprintf('covariate %i',q);

	else
		x    = [1:size(y,1)]';
		Xlab = 'scan';
	end

end

% Plot
%----------------------------------------------------------------------
spm_clf(Finter);
figure(Fgraph)
if BAR

	% organize adjusted data and get mean
	%---------------------------------------------------------------
	i      = any(h')';
	y      = y(i);
	h      = h(i,:);
	x      = [];
	Y      = [];
	for j  = 1:length(y)
		x(j) = find(h(j,:));
	end
	for j  = 1:size(h,2)
		Y(j) = mean(y(h(:,j)));
	end

	% bar chart
	%---------------------------------------------------------------
	[p q]  = bar(Y);
	fill(p,q,[1 1 1]*.9)
	line(x,y,'LineStyle','.','MarkerSize',12)
	xlabel('level')
	ylabel('adjusted response')
	axis([0 (size(h,2) + 1) (min(y) - 1) (max(y) + 1)])
	axis square

else

	% plot
	%---------------------------------------------------------------
	[p q] = sort(-x);
	if all(diff(x(q)))
		plot(x(q),Y(q),'b',x(q),y(q),':')
	else
		plot(x(q),Y(q),'b')
	end
	axis square
	line(x,y,'LineStyle','.','MarkerSize',MSize)
	xlabel(Xlab)
	ylabel('adjusted and fitted response')

	uicontrol(3,'Style','Pushbutton','String','hold on','Callback',...
	['hold on'],  'Position',...
	[40 60 60 20],'Interruptible','yes');
	uicontrol(3,'Style','Pushbutton','String','hold off','Callback',...
	['hold off'], 'Position',...
	[40 40 60 20],'Interruptible','yes');

end

%----------------------------------------------------------------------
title(TITLE,'FontSize',16)
spm_clf(Finter)
