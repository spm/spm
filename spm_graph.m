function [Y,y,beta,SE] = spm_graph(SPM,VOL,xX,xCon,xSDM,hReg)
% Graphical display of adjusted data
% FORMAT [Y y beta SE] = spm_graph(SPM,VOL,xX,xCon,xSDM,hReg)
%
% SPM    - structure containing SPM, distribution & filtering detals
%        - required fields are:
% .swd   - SPM working directory - directory containing current SPM.mat
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .Ic    - indicies of contrasts (in xCon)
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
% xCon   - Contrast definitions structure
%        - required fields are:
% .c     - contrast vector/matrix
%          ( see spm_SpUtil.m for details of contrast structure... )
%
% xSDM   - structure containing contents of SPM.mat file
%        - required fields are:
% .Vbeta
% .VResMS
%          ( see spm_spm.m for contents... )
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
% spm_graph is a Callback script that uses the structures above to:  (i)
% send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
% 
% a) Contrasts of parameter estimates (e.g.  activations) and their
% standard error.
% 
% b) Fitted and adjusted responses that can be plotted against time,
% scan, or an indicator variable in the design matrix.
% 
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses to an event or epoch that would have been seen in the
% absence of other effects.
% 
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g.  in
% constructing explanatory variables such as in a psychophysiological
% interaction). To remove or correct for specific effects, specify an
% appropriate F contrast and simply plot the fitted (and adjusted)
% responses after selecting that F contrast.  The vectors Y and y in the
% workspace will now be corrected for the effects in the reduced design
% matrix (X0) specified in the contrast manager with the column indices
% (iX0) of the confounds in this adjustment.
% 
% Plotting data:
% All data and graphics use filtered data and residuals.    In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrast previously defined can be plotted.
% This ensures that the parameters plotted are meaningful even when there
% is collinearity among the design matrix subpartitions.
% 
%
%_______________________________________________________________________
% %W% Karl Friston %E%


%-Get Graphics figure handle
%-----------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');


%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('Clear',Fgraph,2);


%-Find nearest voxel [Euclidean distance] in point list & update GUI
%-----------------------------------------------------------------------
if ~length(SPM.XYZmm)
	msgbox('No voxels survive masking & threshold(s)!',...
		sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'help','modal')
	Y = []; y = []; beta = []; SE = [];
	return
end

[xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),SPM.XYZmm);
spm_XYZreg('SetCoords',xyz,hReg);


%-Extract required data from results files
%=======================================================================
cwd  = pwd;				%-Note current working directory
cd(SPM.swd)				%-Temporarily move to results dir

%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%     Therefore, parameters & ResMS should be read from the image files,
%     rather than recomputing them on the basis of the Y.mad data.
%-----------------------------------------------------------------------
if exist(fullfile('.','Y.mad')) ~= 2
	sf_noYwarn([])
	y    = [];

elseif SPM.QQ(i) == 0
	sf_noYwarn
	if spm_input('move to nearest voxel',-1,'y/n',[1 0])
		Q       = find(SPM.QQ);
		[xyz,i] = spm_XYZreg('NearestXYZ',...
			  spm_XYZreg('GetCoords',hReg),SPM.XYZmm(:,Q));
		i       = Q(i);
		y       = spm_extract('Y.mad',SPM.QQ(i));

	else
		y       = [];
	end

else
	y    = spm_extract('Y.mad',SPM.QQ(i));
end

% reset pointer and compute voxel indices
%-----------------------------------------------------------------------
spm_XYZreg('SetCoords',xyz,hReg);
rcp     = VOL.iM(1:3,:)*[xyz;1];

% inference (for xlabel)
%-----------------------------------------------------------------------
Z      = SPM.Z(i);
Pz     = spm_P(1,0,Z,SPM.df,SPM.STAT,1,    SPM.n);
Pu     = spm_P(1,0,Z,SPM.df,SPM.STAT,VOL.R,SPM.n);
STR    = [SPM.STAT sprintf(' = %0.2f, p = %0.3f (%.3f corrected)',Z,Pz,Pu)];



%-Get parameter estimates, ResMS, (compute) fitted data & residuals
%=======================================================================
%-NB: Data in Y.mad is raw, must (re)apply temporal smoothing in K
%     Fitted data & residuals are for temporally smoothed model

%-Parameter estimates: beta = xX.pKX*xX.K*y;
%-----------------------------------------------------------------------
beta  = ones(length(xSDM.Vbeta),1);
for i = 1:length(beta)
	beta(i) = spm_sample_vol(xSDM.Vbeta(i),rcp(1),rcp(2),rcp(3),0);
end

%-Compute residuals
%-----------------------------------------------------------------------
if isempty(y)

	% make R = NaN so it will not be plotted
	%---------------------------------------------------------------
	R   = NaN*ones(size(xX.xKXs.X,1),1);

else
	% residuals
	%---------------------------------------------------------------
	R   = spm_filter('apply',xX.K, y) - xX.xKXs.X*beta;

end

%-Residual mean square: ResMS = sum(R.^2)/xX.trRV;
%-----------------------------------------------------------------------
ResMS = spm_sample_vol(xSDM.VResMS,rcp(1),rcp(2),rcp(3),0);
SE    = sqrt(ResMS*diag(xX.Bcov));
COL   = ['r','b','g','c','y','m','r','b','g','c','y','m'];



%-Return to previous directory
%-----------------------------------------------------------------------
cd(cwd)					%-Go back to original working dir.


%-Plot
%=======================================================================

% find out what to plot
%----------------------------------------------------------------------
Cplot = {	'Contrast of parameter estimates',...
	 	'Fitted and adjusted responses',...
		'event/epoch-related responses'};
if ~length(y), Cplot{2} = 'Fitted responses'; end
if ~isfield(xSDM,'Sess'), Cplot = Cplot(1:2); end
Cp    = spm_input('Plot',-1,'m',Cplot);

% select contrast to plot and compute fitted and adjusted data
%----------------------------------------------------------------------
if Cp < 3

	% determine current contrasts
	%---------------------------------------------------------------
	for i = 1:length(xCon)
		Icstr{i} = xCon(i).name;
	end
	Ic    = spm_input('Which contrast?','!+1','m',Icstr);
	TITLE = {Cplot{Cp} xCon(Ic).name};


	% fitted (corrected) data (Y = X1o*beta)
	%---------------------------------------------------------------
	Y     = spm_FcUtil('Yc',xCon(Ic),xX.xKXs,beta);

	% adjusted data
	%---------------------------------------------------------------
	y     = Y + R;

end



% plot parameter estimates
%----------------------------------------------------------------------
if     Cp == 1

	% comute contrast of parameter estimates and standard error
	%--------------------------------------------------------------
	c     = xCon(Ic).c;
	cbeta = c'*beta;
	SE    = sqrt(ResMS*diag(c'*xX.Bcov*c));

	% bar chart
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	h     = bar(cbeta);
	set(h,'FaceColor',[1 1 1]*.8)
	for j = 1:length(cbeta)
		line([j j],([SE(j) 0 - SE(j)] + cbeta(j)),...
			    'LineWidth',3,'Color','r')
	end
	set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
	XLAB  = {'effect' STR};
	YLAB  = 'effect size';


% all fitted effects or selected effects
%-----------------------------------------------------------------------
elseif Cp == 2

	% get ordinates
	%---------------------------------------------------------------
	Cplot = {	'an explanatory variable',...
			'scan or time',...
			'a user specified ordinate'};

	Cx    = spm_input('plot against','!+1','m',Cplot);

	if     Cx == 1

		str  = 'Which column or effect?';
		x    = xX.xKXs.X(:,spm_input(str,'!+1','m',xX.Xnames));
		XLAB = xX.Xnames{i};

	elseif Cx == 2

		if isfield(xX,'RT') & ~isempty(xX.RT)
			x    = xX.RT*[1:size(Y,1)]';
			XLAB = 'time {seconds}';
		else
			x    = [1:size(Y,1)]';
			XLAB = 'scan number';
		end

	elseif Cx == 3

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
		plot(x(q),y(q),'.b','MarkerSize',8); hold on
		plot(x(q),Y(q),'.r','MarkerSize',16); hold off
		set(gca,'XLim',[-1 1] + get(gca,'XLim'))

	end
	YLAB  = 'response';
	XLAB  = {XLAB STR};



% modeling evoked responses
%----------------------------------------------------------------------
elseif Cp == 3

	
	% get session and trials
	%--------------------------------------------------------------
	ss    = length(xSDM.Sess);
	if ss > 1
		str   = sprintf('which sessions (1 to %d)',ss);
		ss    = spm_input(str,-1,'n','1');
	end
	tr    = length(xSDM.Sess{ss(1)}.name);
	if tr > 1
		str   = sprintf('which trials or conditions (1 to %d)',tr);
		tr    = spm_input(str,-1,'n','1');
	end
	Cplot = {	'fitted response',...
			'fitted response and PSTH',...
			'fitted response +/- standard error of response',...
			'fitted response +/- standard error of onset',...
			'fitted response and adjusted data',...
			'parametric plot'};
	if isempty(y), Cplot = Cplot([1 3 4]); end
	Cp      = spm_input('plot in terms of',-2,'m',Cplot);
	TITLE   = Cplot{Cp};
	YLAB    = 'effect size';
	XLAB{1} = 'peri-stimulus time {secs}';


	% cycle over selected events
	%--------------------------------------------------------------
	dx      = xX.dt;

	% reconstruct response with filtering
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	hold on
	XLim  = 0;
	u     = 1;
	for s = ss
	    for t = tr

		% trial-specific parameters
		%------------------------------------------------------
		i      = xSDM.Sess{s}.row(:);
		j      = xSDM.Sess{s}.col(xSDM.Sess{s}.ind{t});
		Q      = xSDM.Sess{s}.para{t};

		% basis functions, filter and parameters
		%------------------------------------------------------
		B      = beta(j);
		X      = xSDM.Sess{s}.bf{t};
		q      = 1:size(X,1);
		x      = q*dx;
		K{1}   = struct('HChoice',	xX.K{s}.HChoice,...
				'HParam',	xX.K{s}.HParam,...
				'LChoice',	xX.K{s}.LChoice,...
				'LParam',	xX.K{s}.LParam,...
				'row',		q,...
				'RT',		dx);

		% fitted responses, adjusted data and standard error
		%------------------------------------------------------

		KX     = spm_filter('apply',K,X);
		Y      = KX*B;
		se     = sqrt(diag(X*xX.Bcov(j,j)*X')*ResMS);
		pst    = xSDM.Sess{s}.pst{t};
		bin    = round(pst/dx);
		q      = find( (bin >= 0) & (bin <= size(X,1)));
		y      = zeros(size(i));
		y(q)   = Y(bin(q));
		y      = y + R(i);

		% onset
		%------------------------------------------------------
		v      = min(find(abs(Y) > max(abs(Y))/2));
		T      = x(v);
		dYdt   = gradient(Y')'/dx;
		seT    = se(v)./dYdt(v);

		% PSTH
		%------------------------------------------------------
		INT      = min(pst):2:max(pst);
		PSTH   = [];
		SEM    = [];
		PST    = [];
		for k  = 1:(length(INT) - 1)
			q = find(pst > INT(k) & pst <= INT(k + 1));
			n = length(q);
			if n
				PSTH = [PSTH mean(y(q))];
				SEM  = [SEM std(y(q))/sqrt(n)];
				PST  = [PST mean(pst(q))];
			end
		end
		

	Cplot = {	'fitted response',...
			'fitted response and PSTH',...
			'fitted response +/- standard error of response',...
			'fitted response +/- standard error of onset',...
			'fitted response and adjusted data',...
			'parametric plot'};

		% plot
		%------------------------------------------------------
		switch TITLE

			case 'fitted response'
			%----------------------------------------------
			plot(x,Y,COL(u))

			case 'fitted response and PSTH'
			%----------------------------------------------
			errorbar(PST,PSTH,SEM,[':' COL(u)])
			plot(PST,PSTH,['.' COL(u)],'MarkerSize',16), hold on
			plot(PST,PSTH,COL(u),'LineWidth',2)
			plot(x,Y,['-.' COL(u)])
			TITLE = 'Peristimulus histogram (2s bins with sem)';

			case 'fitted response +/- standard error of response'
			%----------------------------------------------
			plot(x,Y,COL(u))
			plot(x,Y + se,['-.' COL(u)],x,Y - se,['-.' COL(u)])

			case 'fitted response +/- standard error of onset'
			%----------------------------------------------
			plot(x,Y,COL(u))
			line(([-seT seT] + T),[Y(v) Y(v)],'LineWidth',6)

			case 'fitted response and adjusted data'
			%----------------------------------------------
			plot(x,Y,COL(u),pst,y,['.' COL(u)],'MarkerSize',8)

			case 'parametric plot'
			%----------------------------------------------
			hold off
			surf(x',Q',Q*Y')
			YLAB  = 'parameter';
		end

		% xlabel
		%------------------------------------------------------
		str  = [xSDM.Sess{s}.name{t} sprintf(' (Session %d) - ',s)];
		XLAB{end + 1} = [str COL(u)];
		u    = u + 1;
		XLim = max([XLim max(x)]);
	    end
	end
	hold off; axis on
	set(gca,'XLim',[-4 XLim])
	XLAB{end + 1}  = STR;

end


%-Label and call Plot UI
%----------------------------------------------------------------------
axis square
xlabel(XLAB,'FontSize',10)
ylabel(YLAB,'FontSize',10)
title(TITLE,'FontSize',16)

spm_results_ui('PlotUi',gca)





%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function sf_noYwarn(Q)
%=======================================================================
if nargin < 1, Q = 0; end

if isempty(Q)
	str = {'Graphing unavailable: ',...
		'No time-series data saved with this analysis'};
elseif Q(1) == 0
	str = {'Graphing unavailable: ',...
		'No time-series data saved for this voxel'};
else
	return
end
msgbox(str,sprintf('%s%s: %s...',spm('ver'),...
	spm('GetUser',' (%s)'),mfilename),'error','modal')
