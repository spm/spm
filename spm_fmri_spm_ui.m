function spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT spm_fmri_spm_ui
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual analysis.
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or explanatory variable.
% (e.g. reference waveform, subject).  The paramters are estimated in a
% least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or CONTRAST
% with the t statistic.  The resulting map of t values constitutes the
% SPM{t}.  The SPM{t} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the SPM{t}
% behaves as a smooth stationary Gaussian field.
%
%     From the user's perspective it is important to specify the design
% matrix and contrasts correctly.  The design matrix is built when you
% specify the number of sessions/subjects and conditions.  The covariates 
% (that constitute the columns of the design matrix) can be thought of as
% reference vectors and can be specified as such.  Alternatively one
% can specify reference vectors in terms of response functions or
% waveforms: Waveforms are specified for each EPOCH of scans that
% constitute a particular condition.  Note that if there are two
% conditions, then two covariates are specified each expressing the
% same waveform[s] every time the condition occurs.  A waveform is
% simply a transient response to the onset of an epoch, that lasts for
% the duration of that epoch.  The form of this response waveform can
% be fixed (e.g. a box-car or half sine wave) or allowed to vary
% between conditions by using two (exponentially modulated sine)
% waveforms.  In the latter case there are two covariates for each
% condition (and the CONTRAST must be specified with this in mind).
% These two response functions correspond to the early and late
% components of a modeled response and differential adaptation of the
% hemodynamics during condition-specific epochs can be tested with the
% appropriate CONTRAST (see the final reference below) 
% 
% If you use a box-car (i.e. square wave) function you can optionally
% include its temporal derivative as an additional covariate of
% interest.  The box-car function is delayed by 6 seconds, however this
% may be two much or too little for some brain regions.  The timing
% covariate (temporal derivative) models a small shift in time that best
% fits the data.  Inferences about whether the observed delay is
% significantly different than 6 seconds can be tested using contrasts in
% the usual way (a positive parameter estimate means a shorter delay).
% 
% Covariates of no interest (called confounds) can also be specfied.  You
% will be prompted for some specific confounds such as low frequency
% artifacts (and whole brain activity).
%
% Epochs can vary in length (and order) within and between subjects or runs.
% If multiple subjects or sessions are specified, then subject or run-specific
% waveforms are used.  This means that main effects of conditions and
% interactions between conditions and subjects (or runs) can be evaluated
% with the appropriate contrast.  If you want to treat all your sessions (or
% subjects) as one then specify just one session/subject.
%
% The way that epochs or successive conditions are specified is now more
% intuitive and flexible.  If there are 3 conditions just type in the
% conditions in the order they were presented i.e. 1 2 3 3 2 1 ....
% Later you will be asked to specify the number of scans for each epoch,
% again as a vector (list of numbers).  If the epochs were all the same
% length, then just type in that length once.
%
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% 2:45-53
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage 2:178-181
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% 2:166-172
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165
%
% Josephs O, Turner R and Friston KJ (1997) Event-related fMRI, Hum. Brain
% Map. 0:00-00
%
%___________________________________________________________________________
% %W% Karl Friston, Jean-Baptiste Poline, Christian Buechel %E%



% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
set(Finter,'Name','fMRI analysis'); 

% get design matrix and data
%===========================================================================
if spm_input('Use a pre-specified design?',1,'b','yes|no',[1 0])

	% load pre-specified design matrix
	%-------------------------------------------------------------------
	load(spm_get(1,'.mat','Select fMRIDesMtx.mat'))

	% get filenames
	%-------------------------------------------------------------------
	nsess  = length(Sess);
	nscan  = zeros(1,nsess);
	P      = [];
	if nsess < 16
		for  i = 1:nsess
			nscan(i) = length(Sess{i}.row);
			str      = sprintf('select scans for session %0.0f',i);
			q        = spm_get(nscan(i),'.img',str);
 			P        = strvcat(P,q);
		end
	else
		for  i = 1:nsess
			nscan(i) = length(Sess{i}.row);
		end
		str    = sprintf('select scans for session %0.0f',i);
		P      = spm_get(sum(nscan),'.img',str);
	end

	% Repeat time
	%-------------------------------------------------------------------
	RT     = X.RT;

else

	% get filenames and other user specified parameters
	%-------------------------------------------------------------------
	nsess  = spm_input(['number of sessions'],1,'e',1);
	nscan  = zeros(1,nsess);
	P      = [];
	for  i = 1:nsess
		str      = sprintf('select scans for session %0.0f',i);
		q        = spm_get(Inf,'.img',str);
 		P        = strvcat(P,q);
		nscan(i) = size(q,1);
	end

	% get Repeat time
	%-------------------------------------------------------------------
	RT     = spm_input('Interscan interval {secs}',2);

	% get design matrix
	%-------------------------------------------------------------------
	[X,Sess] = spm_fMRI_design(nscan,RT);

end


% Threshold for F statistic
%---------------------------------------------------------------------------
UFp    = spm_input('threshold for F, p = ',2,'e',0.001);


% Global normalization
%---------------------------------------------------------------------------
str    = 'remove Global effects';
Global = spm_input(str,3,'scale|none',{'Scaling' 'None'});


% Temporal filtering
%===========================================================================

% High-pass filtering
%---------------------------------------------------------------------------
str	= 'Remove low frequencies?';
cLFmenu = {'64 second cut-off',...
	   'specify',...
	   'no'};
cLF     = spm_input(str,4,'m',cLFmenu);


% default based on peristimulus time
% param = cut-off period (max = 256, min = 32)
%---------------------------------------------------------------------------
param   = 256*ones(1,nsess);
for   i = 1:nsess
	for j = 1:length(Sess{i}.pst)
		param(i) = min([param(i) 2*max(Sess{i}.pst{j})]);
	end
end
param   = round(param);
param(param < 32) = 32;


% specify cut-off (default based on peristimulus time)
%---------------------------------------------------------------------------
switch cLF

	case 1
	%-------------------------------------------------------------------
	param   = 64*ones(1,nsess);

	case 2
	%-------------------------------------------------------------------
	str     = 'Cut off period[s] for each session';
	param   = spm_input(str,4,'e',param);
	if length(param) == 1
		 param = param*ones(1,nsess);
	end

	case 3
	%-------------------------------------------------------------------
	param   = param*Inf;
end

% create filterLF struct array
%---------------------------------------------------------------------------
for i = 1:nsess
	filterLF{i} = struct('Choice',cLFmenu{cLF},'Param',param(i));
end



% Low-pass filtering
%---------------------------------------------------------------------------
str     = 'remove high frequencies?';
cHFmenu = {'smooth with hrf',...
	   'smooth with hrf derivative',...
	   'smooth with Gaussian kernel',...
	   'none'};
cHF     = spm_input(str,5,'m',cHFmenu);
param   = [];

% get Gaussian parameter
%---------------------------------------------------------------------------
switch cHF

	case 3
	%-------------------------------------------------------------------
	param = spm_input('Gaussian FWHM (secs)',5,'e',4);
	param = param/sqrt(8*log(2)); 
end

% create filterHF struct
%---------------------------------------------------------------------------
filterHF = struct('Choice',cHFmenu{cHF},'Param',param);


% intrinsic autocorrelations (Vi)
%---------------------------------------------------------------------------
str     = 'Form of intrinsic autocorrelations?';
cVimenu = {'None',...
	   '1/f'};
cVi     = spm_input(str,6,'m',cVimenu);

switch cVi

	case 1
	%-------------------------------------------------------------------
	param   = [];
	Vi      = speye(sum(nscan));

	case 2
	%-------------------------------------------------------------------
	param   = 0.0178;
	Vi      = [];
	for   i = 1:nsess
		k      = nscan(i);
		[x y]  = size(Vi);
		q      = spm_Vintrinsic(k,RT,'1/f',param);
		Vi(([1:k] + x),([1:k] + y)) = q;
	end
	Vi      = sparse(Vi);

end

% create Vi struct
%---------------------------------------------------------------------------
xVi    = struct('Vi',Vi,'Form',cVimenu{cVi},'Param',param);


% the interactive parts of spm_spm_ui are now finished
%---------------------------------------------------------------------------
set(Finter,'Name','thankyou','Pointer','Watch')



% Contruct convolution matrix
%===========================================================================
K     = [];
for i = 1:nsess
	k      = nscan(i);
	[x y]  = size(K);
	q      = spm_make_filter(k,RT,filterHF,filterLF{i});
	K(([1:k] + x),([1:k] + y)) = q;
end
K     = sparse(K);


% get file identifiers and Global values
%===========================================================================
VY     = spm_vol(P);

if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0])
	error('images do not all have the same dimensions'),           end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end


%-Compute Global variate
%---------------------------------------------------------------------------
GM     = 100;
q      = sum(nscan);
g      = zeros(q,1);
for i  = 1:q, g(i) = spm_global(VY(i)); end

% scale if specified (otherwise session specific grand mean scaling)
%---------------------------------------------------------------------------
gSF    = GM./g;
if strcmp(Global,'None')
	for i = 1:nsess
		j      = Sess{i}.row;
		gSF(j) = GM./mean(g(j));
	end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%---------------------------------------------------------------------------
for i = 1:q, VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); end


%-Masking structure
%---------------------------------------------------------------------------
xM = struct(	'T',	ones(q,1),...
		'TH',	g.*gSF,...
		'I',	0,...
		'VM',	{[]},...
		'xs',	struct('Masking','analysis threshold'));


%-Construct full design matrix (X), parameter names (Xnames),
% and design information structure (xX)
%===========================================================================
Xnames = [X.Xname X.Bname];
xX     = struct(	'X',		[X.xX X.bX],...
			'K',		K,...
			'xVi',		xVi,...
			'RT',		X.RT,...
			'dt',		X.dt,...
			'filterLF',	filterLF,...
			'filterHF',	filterHF,...
			'Xnames',	{Xnames'});



%-Effects designated "of interest" - constuct an F-contrast
%--------------------------------------------------------------------------
Fc = [diff(eye(size(X.xX,2))), zeros(size(X.xX,2) - 1,size(X.bX,2))]';


%-Design description (an nx2 cellstr) - for saving and display
%==========================================================================
sGXcalc  = 'mean voxel value';
sGMsca   = 'session specific';
xsDes    = struct(	'Design',			X.DesN(1),...
			'Basis_functions',		X.DesN(2),...
			'Number_of_sessions',		sprintf('%d',nsess),...
			'Interscan_interval',		sprintf('%0.2f',RT),...
			'High_pass_Filter',		filterHF.Choice,...
			'Low_pass_Filter',		filterLF{1}.Choice,...
			'Intrinsic_correlations',	xVi.Form,...
			'Global_calculation',		sGXcalc,...
			'Grand_mean_scaling',		sGMsca,...
			'Global_normalisation',		Global);
%-global structure
%---------------------------------------------------------------------------
xGX.iGXcalc  = Global{1};
xGX.sGXcalc  = sGXcalc;
xGX.rg       = g;
xGX.sGMsca   = sGMsca;
xGX.GM       = GM;
xGX.gSF      = gSF;


%-Save SPMcfg.mat file
%---------------------------------------------------------------------------
save SPMcfg xsDes VY xX xM xGX Fc Sess

%-Display Design report
%===========================================================================
spm_DesRep('DesMtx',xX,{VY.fname}',xsDes)


%-Analysis Proper
%===========================================================================
spm_clf(Finter);
if spm_input('estimate?',1,'b','yes|no',[1 0])
	spm_spm(VY,xX,xM,Fc,Sess);
end


%-End: Cleanup GUI
%---------------------------------------------------------------------------
spm_clf(Finter); spm('Pointer','Arrow')

