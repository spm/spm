function [Y,y,beta,Bcov] = spm_graph(SPM,VOL,xX,xCon,xSDM,hReg)
% Graphical display of adjusted data
% FORMAT [Y y beta Bcov] = spm_graph(SPM,VOL,xX,xCon,xSDM,hReg)
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
% beta   - parameter estimates (ML or MAP)
% Bcov   - Covariance of parameter estimates (ML or conditional)
%
% see spm_getSPM for details
%_______________________________________________________________________
%
% spm_graph is a Callback script that uses the structures above to:  (i)
% send adjusted (y) and fitted data (Y), for the selected voxel, to the
% workspace and (ii) provide graphics for:
% 
% a) Contrasts of parameter estimates (e.g. activations) and their
% standard error.
% 
% b) Fitted and adjusted responses that can be plotted against time,
% scan, or an indicator variable in the design matrix.
% 
% c) (fMRI only).  Evoked responses using the basis functions to give
% impulse responses that would have been seen in the absence of other effects.
% 
% Getting adjusted data:
% Ensuring the data are adjusted properly can be important (e.g. in
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
% the value of one.  Only contrasts previously defined can be plotted.
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
	Y = []; y = []; beta = []; Bcov = [];
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
		Y = []; y = []; beta = []; Bcov = [];
		return
	end
else
	y    = spm_extract(fullfile(SPM.swd,'Y.mad'),SPM.QQ(i));
end

%-Reset pointer, compute voxel indices, compute location string
%-----------------------------------------------------------------------
spm_XYZreg('SetCoords',xyz,hReg);
rcp    = VOL.iM(1:3,:)*[xyz;1];
XYZstr = sprintf(' at [%g, %g, %g]',xyz);


%-Get parameter estimates, ResMS, (compute) fitted data & residuals
%=======================================================================

%-Parameter estimates: beta = xX.pKX*xX.K*y;
%-----------------------------------------------------------------------
beta   = ones(length(xSDM.Vbeta),1);
for  i = 1:length(beta)
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
	R   = spm_sp('r',xX.xKXs,spm_filter(xX.K,y));

end

%-Residual mean square: ResMS = sum(R.^2)/xX.trRV; or Cov(b|y)
%-----------------------------------------------------------------------
if SPM.STAT ~= 'P'
	ResMS = spm_sample_vol(xSDM.VResMS,rcp(1),rcp(2),rcp(3),0);
	Bcov  = ResMS*xX.Bcov;

else
	% hyperparameter and Taylor approximation
	%--------------------------------------------------------------
	Bcov  = xSDM.PPM.Cby;
	for j = 1:length(xSDM.PPM.l)

		l    = spm_sample_vol(xSDM.VHp(j),rcp(1),rcp(2),rcp(3),0);
		Bcov = Bcov + xSDM.PPM.dC{j}*(l - xSDM.PPM.l(j));
	end
end
CI    = 1.6449;					% = spm_invNcdf(1 - 0.05);


%-Colour specifications and index;
%-----------------------------------------------------------------------
Col   = [0 0 0; .8 .8 .8; 1 .5 .5];

%-Plot
%=======================================================================

% find out what to plot
%-----------------------------------------------------------------------
Cplot = {	'Contrast estimates and 90% C.I.',...
	 	'Fitted responses',...
	 	'Event-related responses',...
	 	'Parametric responses',...
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

% select contrast if
%----------------------------------------------------------------------
case {'Contrast estimates and 90% C.I.','Fitted responses'}

	% determine which contrast
	%---------------------------------------------------------------
	Ic    = spm_input('Which contrast?','!+1','m',{xCon.name});
	TITLE = {Cplot xCon(Ic).name};
	if SPM.STAT == 'P'
		TITLE = {Cplot xCon(Ic).name '(conditional estimates)'};
	end


% select session and trial if
%----------------------------------------------------------------------
case {'Event-related responses','Parametric responses','Volterra Kernels'}

	% get session
	%--------------------------------------------------------------
	s     = length(Sess);
	if  s > 1
		s = spm_input('which session','+1','n1',1,s);
	end

	% effect names
	%--------------------------------------------------------------
	switch Cplot
	case 'Volterra Kernels'
		u = length(Sess{s}.Fcname);
	otherwise
		u = length(Sess{s}.U);
	end
	Uname = {};
	for i = 1:u
		Uname{i} = Sess{s}.Fcname{i};
	end

	% get effect
	%--------------------------------------------------------------
	str   = sprintf('which effect');
	u     = spm_input(str,'+1','m',Uname);

	% bin size
	%--------------------------------------------------------------
	dt    = Sess{s}.U{1}.dt;

end

switch Cplot

% plot parameter estimates
%----------------------------------------------------------------------
case 'Contrast estimates and 90% C.I.'

	% compute contrast of parameter estimates and 90% C.I.
	%--------------------------------------------------------------
	cbeta = xCon(Ic).c'*beta;
	CI    = CI*sqrt(diag(xCon(Ic).c'*Bcov*xCon(Ic).c));

	% bar chart
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	cla
	hold on

	% estimates
	%--------------------------------------------------------------
	h     = bar(cbeta);
	set(h,'FaceColor',Col(2,:))

	% standard error
	%--------------------------------------------------------------
	for j = 1:length(cbeta)
		line([j j],([CI(j) 0 - CI(j)] + cbeta(j)),...
			    'LineWidth',6,'Color',Col(3,:))
	end

	title(TITLE,'FontSize',12)
	xlabel('contrast')
	ylabel(['contrast estimate',XYZstr])
	set(gca,'XLim',[0.4 (length(cbeta) + 0.6)])
	hold off


% all fitted effects or selected effects
%-----------------------------------------------------------------------
case 'Fitted responses'

	% predicted or adjusted response
	%---------------------------------------------------------------
	str   = 'predicted or adjusted response?';
	if spm_input(str,'!+1','b',{'predicted','adjusted'},[1 0]);

		% fitted (predicted) data (Y = X1*beta)
		%--------------------------------------------------------
		Y = xX.X*xCon(Ic).c*pinv(xCon(Ic).c)*beta;
	else

		% fitted (adjusted) data (Y = X1o*beta)
		%-------------------------------------------------------
		Y = spm_FcUtil('Yc',xCon(Ic),xX.xKXs,beta);

	end

	% adjusted data
	%---------------------------------------------------------------
	y     = Y  + R;

	% get ordinates
	%---------------------------------------------------------------
	Xplot = {	'an explanatory variable',...
			'scan or time',...
			'a user specified ordinate'};
	Cx    = spm_input('plot against','!+1','m',Xplot);

	% an explanatory variable
	%---------------------------------------------------------------
	if     Cx == 1

		str  = 'Which explanatory variable?';
		i    = spm_input(str,'!+1','m',xX.Xnames);
		x    = xX.xKXs.X(:,i);
		XLAB = xX.Xnames{i};

	% scan or time
	%---------------------------------------------------------------
	elseif Cx == 2

		if isfield(xX,'RT') & ~isempty(xX.RT)
			x    = xX.RT*[1:size(Y,1)]';
			XLAB = 'time {seconds}';
		else
			x    = [1:size(Y,1)]';
			XLAB = 'scan number';
		end

	% user specified
	%---------------------------------------------------------------
	elseif Cx == 3

		x    = spm_input('enter ordinate','!+1','e','',size(Y,1));
		XLAB = 'ordinate';

	end

	% plot
	%---------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	cla
	hold on
	[p q] = sort(x);
	if all(diff(x(q)))
		plot(x(q),Y(q),'LineWidth',4,'Color',Col(2,:));
		plot(x(q),y(q),':','Color',Col(1,:));
		plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(3,:)); 

	else
		plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
		plot(x(q),Y(q),'.','MarkerSize',16,'Color',Col(1,:));
		xlim = get(gca,'XLim');
		xlim = [-1 1]*diff(xlim)/4 + xlim;
		set(gca,'XLim',xlim)

	end
	title(TITLE,'FontSize',12)
	xlabel(XLAB)
	ylabel(['response',XYZstr])
	hold off

% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Event-related responses'

	% get plot type
	%--------------------------------------------------------------
	Rplot   = {	'fitted response and PSTH',...
			'fitted response and 90% C.I.',...
			'fitted response and adjusted data'};

	if isempty(y)
		TITLE = Rplot{2};
	else
		TITLE = Rplot{spm_input('plot in terms of','+1','m',Rplot)};
	end

	% plot
	%--------------------------------------------------------------
	switch TITLE
	case 'fitted response and PSTH'
			str = 'bin size for PSTH {secs}';
			BIN = spm_input(str,'+1','r','2',1);

	otherwise,	BIN = 2;   end

	% basis functions and parameters
	%--------------------------------------------------------------
	X     = Sess{s}.bf/dt;
	x     = ([1:size(X,1)] - 1)*dt;
	j     = Sess{s}.col(Sess{s}.Fci{u}(1:size(X,2)));
	B     = beta(j);
	
	% fitted responses with standard error
	%--------------------------------------------------------------
	Y     = X*B;
	CI    = CI*sqrt(diag(X*Bcov(j,j)*X'));

	% peristimulus times and adjusted data (y = Y + R)
	%--------------------------------------------------------------
	p     = Sess{s}.U{u}.pst;
	bin   = round(p/dt);
	q     = find((bin >= 0) & (bin < size(X,1)));
	pst   = p;
	p     = R(Sess{s}.row(:));
	p(q)  = p(q) + Y(bin(q) + 1);
	y     = p;

	% PSTH
	%--------------------------------------------------------------
	INT   = -BIN:BIN:max(pst);
	PSTH  = [];
	SE    = [];
	PST   = [];
	for k = 1:(length(INT) - 1)
		q = find(pst > INT(k) & pst <= INT(k + 1));
		n = length(q);
		if n
			PSTH = [PSTH mean(y(q))];
			SE   = [SE  std(y(q))/sqrt(n)];
			PST  = [PST mean(pst(q))];
		end
	end

	% plot
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	hold on
	switch TITLE

		case 'fitted response and PSTH'
		%------------------------------------------------------
		errorbar(PST,PSTH,SE)
		plot(PST,PSTH,'LineWidth',4,'Color',Col(2,:))
		plot(x,Y,'-.','Color',Col(3,:))

		case 'fitted response and 90% C.I.'
		%------------------------------------------------------
		plot(x,Y,'Color',Col(2,:),'LineWidth',4)
		plot(x,Y + CI,'-.',x,Y - CI,'-.','Color',Col(1,:))

		case 'fitted response and adjusted data'
		%------------------------------------------------------
		plot(x,Y,'Color',Col(2,:),'LineWidth',4)
		plot(pst,y,'.','MarkerSize',4,'Color',Col(1,:))

	end

	% label
	%-------------------------------------------------------------
	[i j] = max(Y);
	text(ceil(1.1*x(j)),i,Sess{s}.Fcname{u},'FontSize',8);
	title(TITLE,'FontSize',12)
	xlabel('peristimulus time {secs}')
	ylabel(['response',XYZstr])
	hold off

% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Parametric responses'

	% basis functions
	%--------------------------------------------------------------
	bf    = Sess{s}.bf;
	pst   = ([1:size(bf,1)] - 1)*dt;

	% parameteric variable
	%--------------------------------------------------------------
	str   = 'which parameter';
	p     = spm_input(str,'+1','m',Sess{s}.U{u}.Pname);
	ons   = find(Sess{s}.U{u}.u(:,1));
	P     = Sess{s}.U{u}.P(1:length(ons),p);

	% parameter estimates and fitted response
	%--------------------------------------------------------------
	B     = beta(Sess{s}.Fci{u});

	% reconstruct trial-specific responses
	%--------------------------------------------------------------
	Y     = zeros(size(bf,1),length(ons));
	uj    = Sess{s}.U{u}.Pi{p};
	for i = 1:length(ons)
		U      = sparse(1,size(Sess{s}.U{u}.u,2));
		U(uj)  = Sess{s}.U{u}.u(ons(i),uj);
		X      = kron(U,bf);
		Y(:,i) = X*B;
	end
	[P j] = sort(P);
	Y     = Y(:,j);

	% plot
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,2,3)
	surf(pst,P,Y')
	shading flat
	title(Sess{s}.U{u}.Uname{1},'FontSize',12)
	xlabel('PST {secs}')
	ylabel(Sess{s}.U{u}.Pname{p})
	zlabel(['responses',XYZstr])
	axis square

	% plot
	%--------------------------------------------------------------
	subplot(2,2,4)
	[j i] = max(mean(Y,2));
	plot(P,Y(i,:),'LineWidth',4,'Color',Col(2,:))
	str   = sprintf('response at %0.1fs',i*dt);
	title(str,'FontSize',12)
	xlabel(Sess{s}.U{u}.Pname{p})
	axis square
	grid on


% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Volterra Kernels'

	% Parameter estimates and basis functions
	%------------------------------------------------------
	bf    = Sess{s}.bf/dt;
	pst   = ([1:size(bf,1)] - 1)*dt;

	% second order kernel
	%--------------------------------------------------------------
	if u > length(Sess{s}.U)

		% Parameter estimates and kernel
		%------------------------------------------------------
		B     = beta(Sess{s}.Fci{u});
		i     = 1;
		Y     = 0;
		for p = 1:size(bf,2)
		for q = 1:size(bf,2)
     			Y = Y + B(i)*bf(:,p)*bf(:,q)';
			i = i + 1;
		end
		end

		% plot
		%------------------------------------------------------
		figure(Fgraph)
		subplot(2,2,3)
		imagesc(pst,pst,Y)
		axis xy
		axis image

		title('2nd order Kernel','FontSize',12);
		xlabel('perstimulus time {secs}')
		ylabel('perstimulus time {secs}')

		subplot(2,2,4)
		plot(pst,Y)
		axis square
		grid on

		title(Sess{s}.Fcname{u},'FontSize',12);
		xlabel('perstimulus time {secs}')


	% first  order kernel
	%--------------------------------------------------------------
	else
		B     = beta(Sess{s}.Fci{u}(1:size(bf,2)));
		Y     = bf*B;

		% plot
		%------------------------------------------------------
		figure(Fgraph)
		subplot(2,1,2)
		plot(pst,Y)
		grid on
		axis square

		title({'1st order Volterra Kernel' Sess{s}.Fcname{u}},...
			'FontSize',12);
		xlabel('perstimulus time {secs}')
		ylabel(['impluse response',XYZstr])
	end

end


%-call Plot UI
%----------------------------------------------------------------------
spm_results_ui('PlotUi',gca)
