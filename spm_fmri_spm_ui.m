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
% @(#)spm_fmri_spm_ui.m	2.4 Karl Friston, Jean-Baptiste Poline, Christian Buechel 99/01/11



% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');


% get filenames and other user specified parameters
%===========================================================================
set(Finter,'Name','fMRI analysis'); 
nsess  = spm_input(['number of sessions'],1,'e',1);
nscan  = zeros(1,nsess);
Q      = [];
for  i = 1:nsess
	str      = sprintf('select scans for session %0.0f',i);
	P        = spm_get(Inf,'.img',str);
 	Q        = str2mat(Q,P);
	nscan(i) = size(P,1);
end
Q(1,:) = [];
P      = Q;
clear Q;

% Threshold for F statistic
%---------------------------------------------------------------------------
global UFp
UFp    = spm_input('threshold for F, p = ',2,'e',0.001);


% Global normalization
%---------------------------------------------------------------------------
str    = 'remove Global effects';
Global = spm_input(str,3,'scale|none',{'Scaling' 'None'});


% get Repeat time
%---------------------------------------------------------------------------
RT     = spm_input('Interscan interval {secs}',4);	% interscan interval


% get design matrix
%===========================================================================
if spm_input('Design Matrix: ',1,'b','create|load',[1 0])
	[X,Sess] = spm_fMRI_design(nscan,RT);
else
	load(spm_get(1,'.mat','Select fMRIDesMtx.mat'))
	if ~exist('Sess'); Sess = []; end
end



% temporal filtering input: after spm_fMRI_design to get the inter stimuli 
% periods taken from Sess if it exists.
%---------------------------------------------------------------------------

%---------------------------------
str	= 'Cut low frequencies (LF) ?';
cLFmenu   = {'Cut LF with defaults',...
	     'None',...
	     'Specify'};
cLF      = spm_input(str,4,'m',cLFmenu);
cLFstr   = deblank(cLFmenu{cLF});
dft_no_sess = 120;
dft_min_period = 32; 					%-- 2*hrf duration for minimum period 
cLF_period	= [];


if cLF == 1 							%-- Cut LF with defaults

	if size(Sess)
	   for  i = 1:nsess
	     	cLF_period = [cLF_period ...
			  	max(dft_min_period, 2*max(cat(2,Sess{i}.pst{:})))];
	   end
	else 
	   	cLF_period = dft_no_sess*[1:nsess]; 	%- no session, default
	end
	
elseif cLF == 2 						%-- Cut none

	cLF_period = 2*max(nscan)*RT*[1:nsess];

else 										%-- Specify

	cLF_period = [];
	while (prod(size(cLF_period)) ~= nsess) & (prod(size(cLF_period)) ~= 1)
	   cLF_period = spm_input('Cut off period(s)',4,'e',dft_no_sess);
	end
	if prod(size(cLF_period)) == 1, cLF_period = cLF_period*ones(1,nsess); end;

end

filterLF(1:nsess) = struct('Menu',{cLFmenu},'Choice',cLF,'Param',cLF_period);


%---------------------------------
str    = 'Cut high frequencies (HF) ? ';
cHFmenu   = {'Cut HF with default hrf',...
	     'Cut HF with hrf derivative',...
	     'None',...
	     'Specify with gaussian filter',...
	     'Specify with cutting frequency'};
cHF      = spm_input(str,4,'m',cHFmenu);
cHFstr   = deblank(cHFmenu{cHF});
sigma = 0;
param	= [];


switch cHF

	case 4
		sigma = spm_input('Gaussian filter FWHM',4,'e',4,1);
		sigma = sigma/sqrt(8*log(2))/RT; % sigma in scan 
		param = ones(1,nsess)*sigma;
	case 5
		cutF = -1;
		while cutF <= 0
			cutF = spm_input('Cut frequency in Hz',4,'e',1/RT/4,1);
		end
		param = ones(1,nsess)*cutF;
end

filterHF(1:nsess) = struct('Menu',{cHFmenu},'Choice',cHF,'Param',param);


%---------------------------------
% Vi : intrinsic autocorrelation
% Vi_par : parameters of form if assumed
% Vi_Flag : 1 if formed assumed 0 otherwise

str    = 'Assume an intrinsic autocorrelation ? ';
if spm_input(str,4,'yes|no',[1 0])
	Vi_par = [0.0178];

	Vi     = [];
	for i = 1:nsess
	    k      = nscan(i);
	    [x y]  = size(Vi);
	    Vi(([1:k] + x),([1:k] + y)) = spm_Vintrinsic(k,RT,'1/f',Vi_par);
	end
	Vi     = sparse(Vi);
	Vi_Flag = 1;	
else 
	Vi_par = [];
	Vi = speye(sum(nscan));
	Vi_Flag = 0;
end
xVi = struct('Vi',Vi,'flag',Vi_Flag,'Param',Vi_par);


% the interactive parts of spm_spm_ui are now finished
%---------------------------------------------------------------------------
set(Finter,'Name','thankyou','Pointer','Watch')



% Contruct convolution matrix
%===========================================================================


K     = [];
for i = 1:nsess
	k      = nscan(i);
	[x y]  = size(K);
	%- K(([1:k] + x),([1:k] + y)) = spm_sptop(sigma,k);
	K(([1:k] + x),([1:k] + y)) = spm_make_filter(k,RT,filterHF(i),filterLF(i),'norm');
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
		j      = find(X.bX(:,i));
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
Xnames = [X.Xname X.Bname X.Cname];
xX     = struct(	'X',		[X.xX X.bX X.cX],...
			'K',		K,...
			'xVi',		xVi,...
			'RT',		RT,...
			'dt',		X.dt,...
			'sigma',	sigma,...
			'iB',		[1:nsess] + size(X.xX,2),...
			'Xnames',	{Xnames'});



%-Pre-specified contrast for F-test on effects of interest
%===========================================================================
if ~isempty(X.xX)

	%-Effects designated "of interest" - constuct an F-contrast
	%------------------------------------------------------------------
	Fc = [diff(eye(size(X.xX,2))), zeros(size(X.xX,2) - 1,...
	     (size(X.bX,2) + size(X.cX,2)))]';
end


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
str =  {X.DesN{2};
	sprintf('number of sessions: %d',nsess);
	sprintf('interscan interval: %0.2f secs',RT);
	sprintf('temporal smoothing: %0.2f secs',sigma*sqrt(8*log(2))/RT);
	};

sGXcalc  = {'mean voxel value'};
sGMsca   = {'scaling of overall grand mean'};
xsDes    = struct(	'Design',			{X.DesN{1}},...
			'Global_calculation',		{sGXcalc},...
			'Grand_mean_scaling',		{sGMsca},...
			'Global_normalisation',		{Global},...
			'Parameters',			{str});
%-global structure
%---------------------------------------------------------------------------
xGX.iGXcalc  = {Global};
xGX.sGXcalc  = {sGXcalc};
xGX.rg       = g;
xGX.sGMsca   = {sGMsca};
xGX.GM       = GM;
xGX.gSF      = gSF;


%-Save SPMcfg.mat file
%---------------------------------------------------------------------------
save SPMcfg xsDes VY xX xM xGX Fc Sess

%-Display Design report & "Run" button
%===========================================================================
spm_DesRep('DesMtx',xX.X,xX.Xnames,{VY.fname}',xsDes)


%-Analysis Proper
%===========================================================================
if spm_input('Estimate',1,'b',['now|later'],[1 0]);
	spm_spm(VY,xX,xM,Fc,Sess)
end


%-End: Cleanup GUI
%---------------------------------------------------------------------------
spm_clf(Finter); spm('Pointer','Arrow')

