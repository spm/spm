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
% .QQ    - indices of voxels in Y.mad file
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
%          (see spm_FcUtil.m for details of contrast structure... )
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
% responses after selecting that F contrast.  The vectors Y (fitted)
% and y (adjusted) in the workspace will now be corrected for the
% effects in the reduced design matrix (X0) specified in the contrast
% manager with the column indices (iX0) of the confounds in this
% adjustment.
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
	spm('alert!','No suprathreshold voxels!',mfilename,0);
	Y = []; y = []; beta = []; SE = [];
	return
end

[xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),SPM.XYZmm);
spm_XYZreg('SetCoords',xyz,hReg);


%-Extract required data from results files
%=======================================================================
%-Get (approximate) raw data y from Y.mad file
%-NB: Data in Y.mad file is compressed, and therefore not fully accurate
%     Therefore, parameters & ResMS should be read from the image files,
%     rather than recomputing them on the basis of the Y.mad data.
%-----------------------------------------------------------------------
if exist(fullfile(SPM.swd,'Y.mad')) ~= 2
	spm('alert"',{'No raw data saved with this analysis:',...
			'Data portions of plots will be unavailable...'},...
		mfilename,0,1);
	y    = [];

elseif SPM.QQ(i) == 0
	switch spm_input({'No raw data saved at this location:',...
		'Data portions of plots unavailable at this location.',...
		' ','Jump to the nearest voxel with saved raw data?'},...
		1,'bd',{'jump','stay','cancel'},[],2,mfilename)
	case 'jump'
		Q       = find(SPM.QQ);
		[xyz,i] = spm_XYZreg('NearestXYZ',xyz,SPM.XYZmm(:,Q));
		i       = Q(i);
		y       = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(i));
	case 'stay'
		y       = [];
	case 'cancel'
		Y = []; y = []; beta = []; SE = [];
		return
	end
else
	y    = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(i));
end

%-Reset pointer, compute voxel indices, compute location string
%-----------------------------------------------------------------------
spm_XYZreg('SetCoords',xyz,hReg);
rcp     = VOL.iM(1:3,:)*[xyz;1];
XYZstr  = sprintf(' at [%g, %g, %g]',xyz);


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
	R   = spm_sp('r',xX.xKXs,spm_filter('apply',xX.K,y));

end

%-Residual mean square: ResMS = sum(R.^2)/xX.trRV;
%-----------------------------------------------------------------------
ResMS = spm_sample_vol(xSDM.VResMS,rcp(1),rcp(2),rcp(3),0);
SE    = sqrt(ResMS*diag(xX.Bcov));
COL   = ['r','b','g','c','y','m','r','b','g','c','y','m'];


%-Plot
%=======================================================================

% find out what to plot
%-----------------------------------------------------------------------
Cplot = {	'Contrast of parameter estimates',...
	 	'Fitted and adjusted responses',...
	 	'Event/epoch-related responses',...
	 	'Plots of parametric responses',...
		'Volterra Kernels'};


% ensure options are appropriate
%-----------------------------------------------------------------------
if ~isfield(xSDM,'Sess')

	Cplot = Cplot(1:2);
else
	Sess  = xSDM.Sess;
end
Cp     = spm_input('Plot',-1,'m',Cplot);
Cplot  = Cplot{Cp};

switch Cplot

% select contrast to plot and compute fitted and adjusted data
%----------------------------------------------------------------------
case {'Contrast of parameter estimates','Fitted and adjusted responses'}

	% determine current contrasts
	%---------------------------------------------------------------
	Ic    = spm_input('Which contrast?','!+1','m',{xCon.name});
	TITLE = {Cplot xCon(Ic).name};

	% fitted (corrected) data (Y = X1o*beta)
	%---------------------------------------------------------------
	Y     = spm_FcUtil('Yc',xCon(Ic),xX.xKXs,beta);

	% adjusted data
	%---------------------------------------------------------------
	y     = Y + R;

end

switch Cplot

% plot parameter estimates
%----------------------------------------------------------------------
case 'Contrast of parameter estimates'

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
	XLAB  = 'effect';
	YLAB  = ['size of effect',XYZstr];


% all fitted effects or selected effects
%-----------------------------------------------------------------------
case 'Fitted and adjusted responses'

	% get ordinates
	%---------------------------------------------------------------
	Xplot = {	'an explanatory variable',...
			'scan or time',...
			'a user specified ordinate'};

	Cx    = spm_input('plot against','!+1','m',Xplot);

	if     Cx == 1

		str  = 'Which column or effect?';
		i    = spm_input(str,'!+1','m',xX.Xnames);
		x    = xX.xKXs.X(:,i);
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
		xlim = get(gca,'XLim');
		xlim = [-1 1]*diff(xlim)/4 + xlim;
		set(gca,'XLim',xlim)

	end
	YLAB  = ['response',XYZstr];


% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Event/epoch-related responses'


	% averge over sessions?
	%--------------------------------------------------------------
	ss    = length(Sess);
	if ss > 1

		% determine if the same basis functions have been used
		%------------------------------------------------------
		for s = 1:ss
		    rep     = length(Sess{s}.bf) == length(Sess{1}.bf);
		    if ~rep, break, end
		    for t   = 1:length(Sess{s}.bf)
			rep = all(size(Sess{s}.bf{t}) == size(Sess{1}.bf{t}));
			if ~rep, break, end
			rep = all(all( Sess{s}.bf{t}  == Sess{1}.bf{t} ));
			if ~rep, break, end
		    end
		    if ~rep, break, end
		end

		% average over sessions?
		%------------------------------------------------------
		if rep
			str   = 'average over sessions?';
			rep   = spm_input(str,'+1','y/n',[1 0]);
		end

		% selected sessions
		%------------------------------------------------------
		if rep
			ss    = 1:ss;
		else
			str   = sprintf('which session (1 to %d)',ss);
			ss    = spm_input(str,'+1','n','1',1,ss);
		end
	end

	% get plot type
	%--------------------------------------------------------------
	Rplot = {	'fitted response',...
			'fitted response and PSTH',...
			'fitted response +/- standard error',...
			'fitted response and adjusted data'};

	if isempty(y), Rplot = Rplot([1 3]); end
	Cr      = spm_input('plot in terms of','+1','m',Rplot);
	TITLE   = Rplot{Cr};
	YLAB    = ['response',XYZstr];
	XLAB{1} = 'peri-stimulus time {secs}';


	% get selected trials
	%--------------------------------------------------------------
	tr    = length(Sess{ss(1)}.pst);
	if tr > 1
		str   = sprintf('which trials or conditions (1 to %d)',tr);
		tr    = spm_input(str,'+1','n');
	end


	% plot
	%--------------------------------------------------------------
	switch TITLE
	case 'fitted response and PSTH'
			str = 'bin size for PSTH {secs}';
			BIN = spm_input(str,'+1','r','2',1);

	otherwise,	BIN = 2;   end

	% reconstruct response with filtering
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	hold on
	dx    = xX.dt;
	XLim  = 0;
	u     = 1;
	for t = tr

		for s = ss

			% basis functions, filter and parameters
			%----------------------------------------------
			j    = 1:size(Sess{s}.sf{t},2):length(Sess{s}.ind{t});
			j    = Sess{s}.col(Sess{s}.ind{t}(j));
			B    = beta(j);
			X    = Sess{s}.bf{t};
			q    = 1:size(X,1);
			x    = (q - 1)*dx;
			K{1} = struct(  'HChoice',	'none',...
			 		'HParam',	[],...
					'LChoice',	xX.K{s}.LChoice,...
					'LParam',	xX.K{s}.LParam,...
					'row',		q,...
					'RT',		dx);

			% fitted responses with standard error
			%----------------------------------------------
			KX       = spm_filter('apply',K,X);
			Y(q,s)   = KX*B;
			se(:,s)  = sqrt(diag(X*xX.Bcov(j,j)*X')*ResMS);
		end

		% average over sessions
		%------------------------------------------------------
		Y     = Y*ones(length(ss))/length(ss);

		% peristimulus times and adjusted data (Y + R)
		%------------------------------------------------------
		pst   = [];
		y     = [];
		for s = ss
			p       = Sess{s}.pst{t}(:);
			bin     = round(p/dx);
			q       = find((bin >= 0) & (bin < size(X,1)));
			pst     = [pst; p];
			p       = R(Sess{s}.row(:));
			p(q)    = p(q) + Y(bin(q) + 1);
			y       = [y; p];
		end


		% PSTH
		%------------------------------------------------------
		INT    = -BIN:BIN:max(pst);
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

		% plot
		%------------------------------------------------------
		switch TITLE

			case 'fitted response'
			%----------------------------------------------
			plot(x,Y,COL(u))

			case 'fitted response and PSTH'
			%----------------------------------------------
			errorbar(PST,PSTH,SEM,[':' COL(u)])
			plot(PST,PSTH,['.' COL(u)],'MarkerSize',16)
			plot(PST,PSTH,COL(u),'LineWidth',2)
			plot(x,Y,['-.' COL(u)])

			case 'fitted response +/- standard error'
			%----------------------------------------------
			plot(x,Y,COL(u))
			plot(x,Y + se,['-.' COL(u)],x,Y - se,['-.' COL(u)])

			case 'fitted response and adjusted data'
			%----------------------------------------------
			plot(x,Y,COL(u),pst,y,['.' COL(u)],'MarkerSize',8)

		end

		% xlabel
		%------------------------------------------------------
		XLAB{end + 1} = [Sess{s}.name{t} ' - ' COL(u)];
		u    = u + 1;
		XLim = max([XLim max(x)]);

	end

	hold off; axis on
	set(gca,'XLim',[-4 XLim])

% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Plots of parametric responses'

	% Get session
	%--------------------------------------------------------------
	s     = length(Sess);
	if  s > 1
		s = spm_input('which session','+1','n1',[],s);
	end

	% Get [parametric] trial
	%--------------------------------------------------------------
	Vname = {};
	j     = [];
	for i = 1:length(Sess{s}.Pv)
		if length(Sess{s}.Pv{i})
			Vname{end + 1} = Sess{s}.name{i};
			j              = [j i];
		end
	end
	t     = j(spm_input('which effect','+1','m',Vname));

	% parameter estimates and fitted response
	%-------------------------------------------------------------
	B     = beta(Sess{s}.col(Sess{s}.ind{t}));
	Q     = Sess{s}.Pv{t};
	SF    = Sess{s}.sf{t};
	SF    = SF(find(SF(:,1)),:);
	X     = Sess{s}.bf{t};
	q     = 1:size(X,1);
	x     = q*xX.dt;
	K{1}  = struct(		'HChoice',	'none',...
				'HParam',	[],...
				'LChoice',	xX.K{s}.LChoice,...
				'LParam',	xX.K{s}.LParam,...
				'row',		q,...
				'RT',		xX.dt);

	KX    = spm_filter('apply',K,X);
	p     = size(SF,2);
	b     = [];
	for i = 1:size(KX,2)
		b = [b SF*B([1:p] + (i - 1)*p)];
	end
	Y     = KX*b';

	% plot
	%-------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	surf(x',Q',Y')
	shading flat
	TITLE = Sess{s}.name{t};
	XLAB  = 'perstimulus time (secs)';
	YLAB  = Sess{s}.Pname{t};
	zlabel(['responses',XYZstr]);


% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Volterra Kernels'


	% Get session
	%--------------------------------------------------------------
	s     = length(Sess);
	if  s > 1
		s = spm_input('which session','+1','n1',[],s);
	end

	% Get [non-parametric] trial
	%--------------------------------------------------------------
	Vname = {};
	j     = [];
	for i = 1:length(Sess{s}.name)
		Vname{end + 1} = Sess{s}.name{i};
	end
	t     = spm_input('which effect','+1','m',Vname);

	% Parameter estimates
	%--------------------------------------------------------------
	B     = beta(Sess{s}.col(Sess{s}.ind{t}));

	% plot
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)

	% second order kernel
	%--------------------------------------------------------------
	if iscell(Sess{s}.bf{t})

		Y     = 0;
		for i = 1:length(Sess{s}.bf{t})
			Y = Y + B(i)*Sess{s}.bf{t}{i};
		end
		p     = ([1:size(Y,2)] - 1)*xX.dt;
		q     = ([1:size(Y,1)] - 1)*xX.dt;
		imagesc(p,q,Y)
		axis xy

		TITLE = {'Second order Volterra Kernel' Sess{s}.name{t}};
		XLAB  = 'perstimulus time (secs)';
		YLAB  = 'perstimulus time (secs)';

	% first  order kernel
	%--------------------------------------------------------------
	else
		j     = 1:size(Sess{s}.sf{t},2):length(Sess{s}.ind{t});
		Y     = Sess{s}.bf{t}*B(j);
		p     = ([1:length(Y)] - 1)*xX.dt;
		plot(p,Y)
		grid on

		TITLE = {'First order Volterra Kernel' Sess{s}.name{t}};
		XLAB  = 'perstimulus time (secs)';
		YLAB  = ['responses',XYZstr];

	end
end


%-Label and call Plot UI
%----------------------------------------------------------------------
axis square
if strcmp(get(get(gca,'Children'),'type'),'image')
	axis image
end
xlabel(XLAB,'FontSize',10)
ylabel(YLAB,'FontSize',10)
title(TITLE,'FontSize',16)

spm_results_ui('PlotUi',gca)
