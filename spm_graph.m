function [Y,y,beta,SE] = spm_graph(SPM,VOL,xX,xSDM,hReg)
% Graphical display of adjusted data
% FORMAT [Y y beta SE] = spm_graph(SPM,VOL,xX,xSDM,hReg)
%
% SPM    - structure containing SPM, distribution & filtering detals
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .c     - contrast(s) - in cell array
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .XYZmm - location of voxels {mm}
% .QQ    - indices of volxes in Y.mad file
%
%
% VOL    - structure containing details of volume analysed
%        - required fields are:
% .R     - search Volume {resels}
% .iM    - mm -> voxels matrix
% 
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xSDM   - structure containing contents of SPM.mat file
%        - required fields are:
% .Vbeta
% .VResMS
%          ( see spm_spm.m for contents...
%
% hReg   - handle of MIP register
%
% Y      - fitted   data for the selected voxel
% y      - adjusted data for the selected voxel
% beta   - parameter estimates
% SE     - standard error of parameter estimates
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_graph is a CallBack script that uses the strcutures above to
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


%-Get Graphics figure handle
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('Clear',Fgraph,2);


%-Load ER.mat (event-related) file if it exists
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% str = fullfile(SPM.swd,'ER.mat');
% if exist(str,'file'); load(str); end


%-Find nearest voxel [Euclidean distance] in point list & update GUI
%-----------------------------------------------------------------------
if ~length(SPM.XYZmm)
	msgbox('No voxels survive masking & threshold(s)!',...
		sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'help','modal')
	Y=[]; y=[]; beta=[]; SE=[];
	return
end

[xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),SPM.XYZmm);
spm_XYZreg('SetCoords',xyz,hReg);
rcp     = VOL.iM(1:3,:)*[xyz;1];



%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%     Therefore, parameters & ResMS should be read from the image files,
%     rather than recomputing them on the basis of the Y.mad data.
%-----------------------------------------------------------------------
Yfname = fullfile(SPM.swd,'Y.mad');
if exist(Yfname)~=2
	msgbox(['Graphing unavailable: ',...
		'No timecourse data saved with this analysis'],...
		['SPM: ',mfilename],'error','modal')
	return
elseif SPM.QQ(i)==0
	msgbox(['Graphing unavailable: ',...
		'No timecourse data saved for this voxel'],...
		['SPM: ',mfilename],'error','modal')
	return
else
	y = spm_extract([SPM.swd,'/Y.mad'],SPM.QQ(i));
end


%-Get parameter estimates, ResMS, (compute) fitted data & residuals
%-NB: Data in Y.mad is raw, must (re)apply temporal smoothing in K
%     Fitted data & residuals are for temporally smoothed model.
%-----------------------------------------------------------------------
%-Parameter estimates: beta = xX.pKX * xX.K*y;
beta  = ones(length(xSDM.Vbeta),1);
for i = 1:length(beta)
	beta(i) = spm_sample_vol(xSDM.Vbeta(i),rcp(1),rcp(2),rcp(3),0);
end

Y     = xX.xKXs.X * beta;		%-Fitted data (KYhat)
R     = xX.K*y - Y;			%-Residuals (Kres)

%-Residual mean square: ResMS = sum(R.^2)/xX.trRV;
ResMS = spm_sample_vol(xSDM.VResMS,rcp(1),rcp(2),rcp(3),0);

SE    = sqrt(ResMS*diag(xX.Bcov));	%-S.E. of parameter estimates

COL   = ['r','b','g','c','y','m','r','b','g','c','y','m'];



% inference (for xlabel)
%-----------------------------------------------------------------------
Z      = SPM.Z(i);
Pz     = spm_P(1,0,Z,SPM.df,SPM.STAT,1,    SPM.n);
Pu     = spm_P(1,0,Z,SPM.df,SPM.STAT,VOL.R,SPM.n);
STR    = [SPM.STAT sprintf(' = %0.2f, p = %0.3f (%.3f corrected.)',Z,Pz,Pu)];


% find out what to plot
%----------------------------------------------------------------------
Cplot = {'parameter estimates','responses','event-related responses'};
Cp    = spm_input('Plot',-1,'m',Cplot);
TITLE = Cplot{Cp};


% plot parameter estimates
%----------------------------------------------------------------------
if     Cp == 1

	% specify [contrasts] of parameter estimate to bar
	%--------------------------------------------------------------
	Cplot = {	'all parameters',...
			'parameters specified by contrast',...
			'contrast of parameters'};
	%-Don't offer option 3 for conjunctions or F-contrasts!
	if ( size(SPM.c,2)>1 | size(SPM.c{1},2)>1 ), Cplot(3)=[]; end
	Cp    = spm_input('Estimates to plot',-1,'m',Cplot);
	TITLE = Cplot{Cp};
	XLAB  = 'effect';


	if     Cp == 1
		
		%-beta & SE already OK

	elseif Cp == 2

		tmp  = find(any(cat(2,SPM.c{:}),2));
		beta = beta(tmp);
		SE   = SE(tmp);

	elseif Cp == 3

		%-** Should really read c'b from file (when that's implemented)
		beta = SPM.c{1}'*beta;
		SE   = sqrt(ResMS*SPM.c{1}'*xX.Bcov*SPM.c{1});
		XLAB = 'contrast';

	end


	% bar chart
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	h     = bar(beta);
	set(h,'FaceColor',[1 1 1]*.8)
	for j = 1:length(beta)
		line([j j],([SE(j) 0 - SE(j)] + beta(j)),...
			    'LineWidth',3,'Color','r')
	end



% all fitted effects or selected effects
%-----------------------------------------------------------------------
elseif Cp == 2

	% fitted data
	%---------------------------------------------------------------
	Cplot = {	'all effects',...
			'subspace spanned by contrast',...
			'subspace of the contrast',...
			'specified effects'};
	Cx    = spm_input('Fit',-1,'m',Cplot);
	TITLE = [TITLE,': ',Cplot{Cx}];


	if     Cx == 1
		
		Y    = xX.xKXs.X * beta;

	elseif Cx == 2

		i    = find(any(SPM.c{1},2));
		Y    = xX.xKXs.X(:,i)*beta(i);

	elseif Cx == 3

		X    = xX.xKXs.X*SPM.c{1};
		Y    = xX.xKXs.X*(xX.pKX*y);

	elseif Cx == 4

		str  = sprintf('Which columns or effects [1-%d]?',length(beta));
		i    = spm_input(str,'!+1','n');
		Y    = xX.xKXs.X(:,i)*beta(i);
	end


	% adjusted data
	%---------------------------------------------------------------
	y     = Y + R;


	% get ordinates
	%---------------------------------------------------------------
	Cplot = {	'a column of design matrix',...
			'contrast of design matrix',...
			'scan or time',...
			'a user specified ordinate'};
	%-Don't offer option 2 for conjunctions or F-contrasts!
	i = 1:length(Cplot);
	if (size(SPM.c,2)>1 | size(SPM.c{1},2)>1), i(2)=[]; end

	Cx    = spm_input('plot against',-1,'m',Cplot(i),i);

	if     Cx == 1

		str  = sprintf('Which column or effect [1-%d]?',length(beta));
		i    = spm_input(str,'!+1','n1','',1);
		x    = xX.xKXs.X(:,i);
		XLAB = sprintf('%s (explanatory variable %d)',xX.Xnames{i},i);

	elseif Cx == 2

		x    = xX.xKXs.X*SPM.c{1};
		XLAB = 'Contrast of explanatory variables';

	elseif Cx == 3

		if isfield(xX,'RT') & ~isempty(xX.RT)
			x    = xX.RT*[1:size(Y,1)]';
			XLAB = 'time {seconds}';
		else
			x    = [1:size(Y,1)]';
			XLAB = 'scan number';
		end

	elseif Cx == 4

		x    = spm_input('enter ordinate','!+1','e','',size(Y,1));
		XLAB = 'ordinate';

	end

	% plot
	%---------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	[p q] = sort(x);
	if all(diff(x(q)))
		plot(x(q),y(q),':b'); hold on
		plot(x(q),y(q),'.b','MarkerSize',8); hold on
		plot(x(q),Y(q),'r' ); hold off

	else
		plot(x(q),y(q),'.b', 'MarkerSize',8); hold on
		plot(x(q),Y(q),'.r','MarkerSize',16); hold off

	end
	set(gca,'XLim',[-1 1] + get(gca,'XLim'))



% modeling evoked responses
%----------------------------------------------------------------------
elseif Cp == 3

	j     = 1;
	if size(ERI,2) > 1
		j    = spm_input('which events',-1,'n','1 2');
	end
	Cplot = {	'fitted response',...
			'fitted response and PSTH',...
			'fitted response +/- standard error of response',...
			'fitted response +/- standard error of onset',...
			'fitted response and adjusted data'};
	Cp    = spm_input('plot in terms of','!+1','m',Cplot);
	TITLE = deblank(Cplot(Cp,:));


	% cycle over selected events
	%--------------------------------------------------------------
	dx    = 0.1;
	x     = [0:(size(DER,1) - 1)]'*dx - 4;

	% reconstruct response without smoothing
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	XLAB  = 'peri-stimulus time {secs}';
	hold on
	u     = 1;
	for i = j
		Y      = DER*beta(ERI(:,i));
		se     = sqrt(diag(DER*xX.Bcov(ERI(:,i),ERI(:,i))*DER')*ResMS);
		pst    = PST(:,i);
		bin    = round((pst + 4)/dx);
		q      = find( (bin >= 1) & (bin <= size(DER,1)) & pst);
		bin    = bin(q);
		pst    = pst(q);
		y      = DER(bin,:)*beta(ERI(:,i)) + R(q);
		v      = min(find(abs(Y) > max(abs(Y))/2));
		T      = x(v);
		dYdt   = gradient(Y')'/dx;
		seT    = se(v)./dYdt(v);

		% PSTH
		%------------------------------------------------------
		dBIN   = 2/dx;
		BIN    = 1/dx:dBIN:32/dx;
		PSTH   = zeros(length(BIN) - 1,1);
		SEM    = zeros(length(BIN) - 1,1);
		for k  = 1:(length(BIN) - 1)
			q = find(bin > BIN(k) & bin <= BIN(k + 1));
			n = length(q);
			if n
				PSTH(k) = mean(y(q));
				SEM(k)  = std(y(q))/sqrt(n);
			end
		end
		BIN    = (BIN(1:k) + dBIN/2)*dx - 4;
		

		% plot
		%------------------------------------------------------
		if Cp == 1
			plot(x,Y,COL(u))

		elseif Cp == 2
			errorbar(BIN,PSTH,SEM,[':' COL(u)])
			plot(BIN,PSTH,['.' COL(u)],'MarkerSize',16), hold on
			plot(BIN,PSTH,COL(u),'LineWidth',2)
			plot(x,Y,['-.' COL(u)])
			TITLE = 'Peristimulus histogram (2s bins with sem)';

		elseif Cp == 3
			plot(x,Y,COL(u),x,Y + se,...
				['-.' COL(u)],x,Y - se,['-.' COL(u)])

		elseif Cp == 4
			plot(x,Y,COL(u))
			line(([-seT seT] + T),[Y(v) Y(v)],'LineWidth',6)

		elseif Cp == 5
			plot(x,Y,COL(u),pst,y,['.' COL(u)],...
				'MarkerSize',8,'LineWidth',2)

		end
		XLAB = str2mat(XLAB,[sprintf('Trial type %d - ',i) COL(u)]);
		u    = u + 1;
	end

	hold off; axis on
	set(gca,'XLim',[min(x) max(x)])

end


%-Label and call Plot UI
%----------------------------------------------------------------------
axis square
XLAB      = str2mat(XLAB,STR);
YLAB      = 'effect size';
xlabel(XLAB, 'FontSize',10)
ylabel(YLAB,'FontSize',10)
title(TITLE,'FontSize',16)

spm_results_ui('PlotUi',gca)
