function [Y,y,beta,Bcov] = spm_graph(xSPM,SPM,hReg)
% Graphical display of adjusted data
% FORMAT [Y y beta Bcov] = spm_graph(xSPM,SPM,hReg)
%
% SPM    - structure containing SPM, distributional & filtering details
%	   about the excusion set
% SPM    - structure containing generic details about the analysis
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
% All data and graphics use filtered/whitened data and residuals. In PET
% studies the parameter estimates and the fitted data are often the same
% because the explanatory variables are simply indicator variables taking
% the value of one.  Only contrasts previously defined can be plotted.
% This ensures that the parameters plotted are meaningful even when there
% is collinearity among the design matrix subpartitions.
%
% Selecting contrasts used for PPMs will automatically give plots
% based on conditonal estimates.
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
if ~length(xSPM.XYZmm)
	spm('alert!','No suprathreshold voxels!',mfilename,0);

	Y = []; y = []; beta = []; Bcov = [];
	return
end

[xyz,i] = spm_XYZreg('NearestXYZ',spm_XYZreg('GetCoords',hReg),xSPM.XYZmm);
spm_XYZreg('SetCoords',xyz,hReg);
XYZ     = xSPM.XYZ(:,i);		% coordinates


%-Extract filtered and whitened data from files
%=======================================================================
try
	y = spm_get_data(SPM.xY.VY,XYZ);
	y = spm_filter(SPM.xX.K,SPM.xX.W*y);
catch
	try
		% remap files in SPM.xY.P if SPM.xY.VY is no longer valid
		%-------------------------------------------------------
		SPM.xY.VY = spm_vol(SPM.xY.P);
		y = spm_get_data(SPM.xY.VY,XYZ);
		y = spm_filter(SPM.xX.K,SPM.xX.W*y);

	catch
		% data has been moved or renamed
		%-------------------------------------------------------
		y = [];
		spm('alert!',{'Original data have been moved or renamed',...
			      'Recomendation: please update SPM.xY.P'},...
			       mfilename,0);
	end
end
XYZstr = sprintf(' at [%g, %g, %g]',xyz);


%-Compute residuals
%-----------------------------------------------------------------------
if isempty(y)

	% make R = NaN so it will not be plotted
	%---------------------------------------------------------------
	R   = NaN*ones(size(SPM.xX.X,1),1);

else
	% residuals (non-whitened)
	%---------------------------------------------------------------
	R   = spm_sp('r',SPM.xX.xKXs,y);

end

%-Get parameter and hyperparameter estimates
%=======================================================================
if xSPM.STAT ~= 'P'

	%-Parameter estimates:   beta = xX.pKX*xX.K*y;
	%-Residual mean square: ResMS = sum(R.^2)/xX.trRV
	%---------------------------------------------------------------
	beta  = spm_get_data(SPM.Vbeta, XYZ);
	ResMS = spm_get_data(SPM.VResMS,XYZ);
	Bcov  = ResMS*SPM.xX.Bcov;

else
	% or conditional estimates with
	% Cov(b|y) through Taylor approximation
	%---------------------------------------------------------------
	beta  = spm_get_data(SPM.VCbeta, XYZ);
	Bcov  = SPM.PPM.Cby;
	for j = 1:length(SPM.PPM.l)

		l    = spm_get_data(SPM.VHp(j),XYZ);
		Bcov = Bcov + SPM.PPM.dC{j}*(l - SPM.PPM.l(j));
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
try
	Sess  = SPM.Sess;
catch
	Cplot = Cplot(1:2);	
end
Cplot  = Cplot{spm_input('Plot',-1,'m',Cplot)};

switch Cplot

% select contrast if
%----------------------------------------------------------------------
case {'Contrast estimates and 90% C.I.','Fitted responses'}

	% determine which contrast
	%---------------------------------------------------------------
	Ic    = spm_input('Which contrast?','!+1','m',{SPM.xCon.name});
	TITLE = {Cplot SPM.xCon(Ic).name};
	if xSPM.STAT == 'P'
		TITLE = {Cplot SPM.xCon(Ic).name '(conditional estimates)'};
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
		u = length(Sess(s).Fc);
	otherwise
		u = length(Sess(s).U);
	end
	Uname = {};
	for i = 1:u
		Uname{i} = Sess(s).Fc(i).name;
	end

	% get effect
	%--------------------------------------------------------------
	str   = sprintf('which effect');
	u     = spm_input(str,'+1','m',Uname);

	% bin size
	%--------------------------------------------------------------
	dt    = SPM.xBF.dt;

end

switch Cplot

% plot parameter estimates
%----------------------------------------------------------------------
case 'Contrast estimates and 90% C.I.'

	% compute contrast of parameter estimates and 90% C.I.
	%--------------------------------------------------------------
	cbeta = SPM.xCon(Ic).c'*beta;
	CI    = CI*sqrt(diag(SPM.xCon(Ic).c'*Bcov*SPM.xCon(Ic).c));

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
		Y = SPM.xX.X*SPM.xCon(Ic).c*pinv(SPM.xCon(Ic).c)*beta;
	else

		% fitted (corrected)  data (Y = X1o*beta)
		%-------------------------------------------------------
		Y = spm_FcUtil('Yc',SPM.xCon(Ic),SPM.xX.xKXs,beta);

	end

	% adjusted data
	%---------------------------------------------------------------
	y     = Y + R;

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
		i    = spm_input(str,'!+1','m',SPM.xX.name);
		x    = SPM.xX.xKXs.X(:,i);
		XLAB = SPM.xX.name{i};

	% scan or time
	%---------------------------------------------------------------
	elseif Cx == 2

		if isfield(SPM.xY,'RT')
			x    = SPM.xY.RT*[1:size(Y,1)]';
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
		plot(x(q),Y(q),'.','MarkerSize',16,'Color',Col(1,:));
		plot(x(q),y(q),'.','MarkerSize',8, 'Color',Col(2,:));
		xlim = get(gca,'XLim');
		xlim = [-1 1]*diff(xlim)/4 + xlim;
		set(gca,'XLim',xlim)

	end
	title(TITLE,'FontSize',12)
	xlabel(XLAB)
	ylabel(['response',XYZstr])
	legend('fitted','plus error')
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


		% build a simple FIR model subpartition (X); bin size = TR
		%------------------------------------------------------
		BIN         = SPM.xY.RT;
		xBF         = SPM.xBF;
		U           = Sess(s).U(u);
		U.u         = U.u(:,1);
		xBF.name    = 'Finite Impulse Response';
		xBF.order   = round(32/BIN);
		xBF.length  = xBF.order*BIN;
		xBF         = spm_get_bf(xBF);
		BIN         = xBF.length/xBF.order;
		X           = spm_Volterra(U,xBF.bf,1);
		k           = SPM.nscan(s);
		X           = X([0:(k - 1)]*SPM.xBF.T + SPM.xBF.T0 + 32,:);

		% place X in SPM.xX.X
		%------------------------------------------------------
		jX          = Sess(s).row;
		iX          = Sess(s).col(Sess(s).Fc(u).i);
		iX0         = [1:size(SPM.xX.X,2)];
		iX0(iX)     = [];
		X           = [X SPM.xX.X(jX,iX0)];
		X           = spm_filter(SPM.xX.K(s),X);

		% Re-estimate to get PSTH and CI
		%------------------------------------------------------
		j           = xBF.order;
		xX          = spm_sp('Set',X);
		pX          = spm_sp('x-',xX);
		PSTH        = pX*y(jX);
		res         = spm_sp('r',xX,y(jX));
		df          = size(X,1) - size(X,2);
		bcov        = pX*pX'*sum(res.^2)/df;
		PSTH        = PSTH(1:j)/dt;
		PST         = [1:j]*BIN - BIN/2;
		PCI         = CI*sqrt(diag(bcov(1:j,(1:j))))/dt;
	end

	% basis functions and parameters
	%--------------------------------------------------------------
	X     = SPM.xBF.bf/dt;
	x     = ([1:size(X,1)] - 1)*dt;
	j     = Sess(s).col(Sess(s).Fc(u).i(1:size(X,2)));
	B     = beta(j);
	
	% fitted responses with standard error
	%--------------------------------------------------------------
	Y     = X*B;
	CI    = CI*sqrt(diag(X*Bcov(j,j)*X'));

	% peristimulus times and adjusted data (y = Y + R)
	%--------------------------------------------------------------
	pst   = Sess(s).U(u).pst;
	bin   = round(pst/dt);
	q     = find((bin >= 0) & (bin < size(X,1)));
	y     = R(Sess(s).row(:));
	pst   = pst(q);
	y     = y(q) + Y(bin(q) + 1);

	% plot
	%--------------------------------------------------------------
	figure(Fgraph)
	subplot(2,1,2)
	hold on
	switch TITLE

		case 'fitted response and PSTH'
		%------------------------------------------------------
		errorbar(PST,PSTH,PCI)
		plot(PST,PSTH,'LineWidth',4,'Color',Col(2,:))
		plot(x,Y,'-.','Color',Col(3,:))

		case 'fitted response and 90% C.I.'
		%------------------------------------------------------
		plot(x,Y,'Color',Col(2,:),'LineWidth',4)
		plot(x,Y + CI,'-.',x,Y - CI,'-.','Color',Col(1,:))

		case 'fitted response and adjusted data'
		%------------------------------------------------------
		plot(x,Y,'Color',Col(2,:),'LineWidth',4)
		plot(pst,y,'.','Color',Col(3,:))

	end

	% label
	%-------------------------------------------------------------
	[i j] = max(Y);
	text(ceil(1.1*x(j)),i,Sess(s).Fc(u).name,'FontSize',8);
	title(TITLE,'FontSize',12)
	xlabel('peristimulus time {secs}')
	ylabel(['response',XYZstr])
	hold off


% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Parametric responses'


	% return gracefully if no parameters
	%--------------------------------------------------------------
	if ~Sess(s).U(u).P(1).h, break, end

	% basis functions
	%--------------------------------------------------------------
	bf    = SPM.xBF.bf;
	pst   = ([1:size(bf,1)] - 1)*dt;

	% orthogonalised expansion of parameteric variable
	%--------------------------------------------------------------
	str   = 'which parameter';
	p     = spm_input(str,'+1','m',cat(2,Sess(s).U(u).P.name));
	P     = Sess(s).U(u).P(p).P;
	q     = [];
	for i = 0:Sess(s).U(u).P(p).h;
		q = [q spm_en(P).^i];
	end
	q     = spm_orth(q);


	% parameter estimates for this effect
	%--------------------------------------------------------------
	B     = beta(Sess(s).Fc(u).i);

	% reconstruct trial-specific responses
	%--------------------------------------------------------------
	Y     = zeros(size(bf,1),size(q,1));
	uj    = Sess(s).U(u).P(p).i;
	for i = 1:size(P,1)
		U      = sparse(1,uj,q(i,:),1,size(Sess(s).U(u).u,2));
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
	title(Sess(s).U(u).name{1},'FontSize',12)
	xlabel('PST {secs}')
	ylabel(Sess(s).U(u).P(p).name)
	zlabel(['responses',XYZstr])
	axis square

	% plot
	%--------------------------------------------------------------
	subplot(2,2,4)
	[j i] = max(mean(Y,2));
	plot(P,Y(i,:),'LineWidth',4,'Color',Col(2,:))
	str   = sprintf('response at %0.1fs',i*dt);
	title(str,'FontSize',12)
	xlabel(Sess(s).U(u).P(p).name)
	axis square
	grid on


% modeling evoked responses based on Sess
%----------------------------------------------------------------------
case 'Volterra Kernels'

	% Parameter estimates and basis functions
	%------------------------------------------------------
	bf    = SPM.xBF.bf/dt;
	pst   = ([1:size(bf,1)] - 1)*dt;

	% second order kernel
	%--------------------------------------------------------------
	if u > length(Sess(s).U)

		% Parameter estimates and kernel
		%------------------------------------------------------
		B     = beta(Sess(s).Fc(u).i);
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

		title(Sess(s).Fc(u).name,'FontSize',12);
		xlabel('perstimulus time {secs}')


	% first  order kernel
	%--------------------------------------------------------------
	else
		B     = beta(Sess(s).Fc(u).i(1:size(bf,2)));
		Y     = bf*B;

		% plot
		%------------------------------------------------------
		figure(Fgraph)
		subplot(2,1,2)
		plot(pst,Y)
		grid on
		axis square

		title({'1st order Volterra Kernel' Sess(s).Fc(u).name},...
			'FontSize',12);
		xlabel('perstimulus time {secs}')
		ylabel(['impluse response',XYZstr])
	end

end


%-call Plot UI
%----------------------------------------------------------------------
spm_results_ui('PlotUi',gca)
