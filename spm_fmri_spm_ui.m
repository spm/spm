function [xX,Sess] = spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT [xX,Sess] = spm_fmri_spm_ui
%
% xX            - structure describing design matrix
% xX.X          - design matrix
% xX.dt         - time bin {secs}
% xX.RT         - Repetition time {secs}
% xX.iH         - vector of H partition (condition effects)      indices,
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects)          indices
% xX.iG         - vector of G partition (nuisance variables)     indices
% xX.Xnames     - cellstr of effect names corresponding to columns
%                 of the design matrix
%
% Sess{s}.BFstr   - basis function description string
% Sess{s}.DSstr   - Design description string
% Sess{s}.row     - scan   indices      for session s
% Sess{s}.col     - effect indices      for session s
% Sess{s}.name{i} - of ith trial type   for session s
% Sess{s}.ind{i}  - column indices      for ith trial type {within session}
% Sess{s}.bf{i}   - basis functions     for ith trial type
% Sess{s}.sf{i}   - stick functions     for ith trial type
% Sess{s}.ons{i}  - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}  - peristimulus times  for ith trial type (secs)
% Sess{s}.para{i} - vector of paramters for ith trial type
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% filtering that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual parameter
% estimation.
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or explanatory variable.
% (e.g. regressor or stimulus function).  The parameters are estimated in
% a least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or contrast
% with the T or F statistic.  The resulting statistical map constitutes 
% an SPM.  The SPM{T}/{F} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the components of
% the SPM (i.e. residual fields) behave as smooth stationary Gaussian fields.
%
% spm_fmri_spm_ui allows you to (i) specify a statistical model in terms
% of a design matrix, (ii) review that design, (iii) associate some data
% with a pre-specified design [or (iv) specify both the data and design]
% and then proceed to estimate the parameters of the model.
% Inferences can be made about the ensuing parameter estimates (at a first
% or fixed-effect level) in the results section, or they can be re-entered
% into a second (random-effect) level analysis by treating the session or 
% subject-specific [contrasts of] parameter estimates as new summary data.
% Inferences at any level obtain by specifying appropriate T or F contrasts
% in the results section to produce SPMs and tables of p values and statistics.
%
% spm_fmri_spm calls spm_fMRI_design which allows you to configure a
% design matrix in terms of events or epochs.  This design matrix can be
% specified before or during data specification.  In some instances
% (e.g. with stochastic designs that have to realized before data
% acquisition) it is necessary to build the design matrix first and then
% select the corresponding data.  In others it may be simpler to specify
% the data and then the design.  Both options are supported.  Once the
% design matrix, data and filtering have been specified spm_fmri_spm_ui
% calls spm_spm to estimate the model parameters that are then saved for
% subsequent analysis.
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, where the latter model prolonged
% and possibly time-varying responses to state-related changes in
% experimental conditions.  Event-related response are modelled in terms
% of responses to instantaneous events.  Mathematically they are both
% modeled by convolving a series of delta (or stick) functions,
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions can be very simple, like a box car,
% or may model voxel-specific forms of evoked responses with a linear
% combination of several basis functions (e.g. a Fourier set).  Basis
% functions can be used to plot estimated responses to single events or
% epochs once the parameters (i.e. basis function coefficients) have
% been estimated.  The importance of basis functions is that they provide
% a graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The
% nice thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be sychronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
% 
% Event-related designs may be stochastic or deterministic.  Stochastic
% designs involve one of a number of trial-types occurring with a
% specified probably at successive intervals in time.  These
% probabilities can be fixed (stationary designs) or time-dependent
% (modulated or non-stationary designs).  The most efficient designs
% obtain when the probabilities of every trial type are equal and this is
% enforced in SPM.  The modulation of non-stationary designs is simply
% sinusoidal with a period of 32 seconds.  A critical aspect of
% stochastic event-related designs is whether to include null events or
% not.  If you wish to estimate the evoke response to a specific event
% type (as opposed to differential responses) then a null event must be
% included (even though it is not modelled explicitly).
% 
% The choice of basis functions depends upon the nature of the inference
% sought.  One important consideration is whether you want to make
% inferences about compounds of parameters (i.e.  contrasts).  This is
% the case if (i) you wish to use a SPM{T} to look separately at
% activations and deactivations or (ii) you with to proceed to a second
% (random-effect) level of analysis.  If this is the case then (for
% event-related studies) use a canonical hemodynamic response function
% (HRF) and derivatives with respect to latency (and dispersion).  Unlike
% other bases, contrasts of these effects have a physical interpretation
% and represent a parsimonious way of characterising event-related
% responses.  Bases such as a Fourier set require the SPM{F} for
% inference and preclude second level analyses.
% 
% See spm_fMRI_design for more details about how designs are specified.
%
% Serial correlations in fast fMRI time-series are dealt with as
% described in spm_spm.  At this stage you need to specific the filtering
% that will be applied to the data (and design matrix).  This filtering
% is important to ensure that bias in estimates of the standard error are
% minimized.  This bias results from a discrepancy between the estimated
% (or assumed) auto-correlation structure of the data and the actual
% intrinsic correlations.  The intrinsic correlations will be estimated
% automatically using an AR(1) model during parameter estimation.  The
% discrepancy between estimated and actual intrinsic (i.e. prior to
% filtering) correlations are greatest at low frequencies.  Therefore
% specification of the high-pass component of the filter is particularly
% important.  High pass filtering is now implemented at the level of the
% filtering matrix K (as opposed to entering as confounds in the design
% matrix).  The default cutoff period is twice the maximum time interval
% between the most frequently occurring event or epoch (i.e the minium of
% all maximum intervals over event or epochs).
%
% N.B.
% Burst Mode is a specialist design for intermittent epochs of acquisitions
% (used for example to allow for intercalated EEG recording).  Each burst
% is treated as a session but consistent within-session effects (e.g. T1
% effects) are modeled in X.bX.  The primary use of this mode is to generate
% parameter estimate images for a second level analysis.
%
%-----------------------------------------------------------------------
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
% %W% Karl Friston, Jean-Baptiste Poline, Christian Buchel %E%
SCCSid  = '%I%';


%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: fMRI analysis',0);
spm_help('!ContextHelp',mfilename)


% get design matrix and/or data
%=======================================================================
MType   = {'specify a model',...
	   'review a specified model',...
	   'estimate a specified model',...
	   'specify and estimate a model'};
str     = 'Would you like to';
MT      = spm_input(str,1,'m',MType);


%-Initialise output arguments in case return early
X    = [];
Sess = [];

switch MT
%-----------------------------------------------------------------------

	case 1
	% specify a design matrix
	%---------------------------------------------------------------
	if sf_abort(1), spm_clf(Finter), return, end
	[xX,Sess] = spm_fMRI_design;
	spm_fMRI_design_show(xX,Sess);
	return

	case 2
	% 'review a specified model'
	%---------------------------------------------------------------
	spm_clf(Finter)
	[xX,Sess] = spm_fMRI_design_show;
	return

	case 3
	% load pre-specified design matrix
	%---------------------------------------------------------------
	if sf_abort([2,3]), spm_clf(Finter), return, end
	load(spm_get(1,'.mat','Select fMRIDesMtx.mat'))


	% get filenames
	%---------------------------------------------------------------
	nsess  = length(xX.iB);
	nscan  = zeros(1,nsess);
	for  i = 1:nsess
		nscan(i) = length(find(xX.X(:,xX.iB(i))));
	end
	P      = [];
	if nsess < 16
		for i = 1:nsess
			str = sprintf('select scans for session %0.0f',i);
			q   = spm_get(nscan(i),'.img',str);
 			P   = strvcat(P,q);
		end
	else
		str   = sprintf('select scans for this study');
		P     = spm_get(sum(nscan),'.img',str);
	end

	% Repeat time
	%---------------------------------------------------------------
	RT     = xX.RT;


	case 4
	% get filenames and design matrix
	%---------------------------------------------------------------
	if sf_abort, spm_clf(Finter), return, end
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
	%---------------------------------------------------------------
	RT        = spm_input('Interscan interval {secs}',2);

	% get design matrix
	%---------------------------------------------------------------
	[xX,Sess] = spm_fMRI_design(nscan,RT);

end

% Assemble other deisgn parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)

% get rows
%-----------------------------------------------------------------------
for i = 1:nsess
	row{i} = find(xX.X(:,xX.iB(i)));
end
BFstr  = Sess{1}.BFstr;
DSstr  = Sess{1}.DSstr;


% Global normalization
%-----------------------------------------------------------------------
str    = 'remove Global effects';
Global = spm_input(str,1,'scale|none',{'Scaling' 'None'});


% Burst mode
%-----------------------------------------------------------------------
if length(nscan) > 16 & ~any(diff(nscan))
	BM    = spm_input('Burst mode','+1','y/n',[1 0]);
else
	BM    = 0;
end

% Model scan effects as oppsed to session effects if burst mode
%-----------------------------------------------------------------------
if BM
	k         = nscan(1);
	l         = length(xX.iC);
	xX.X      = xX.X(:,xX.iC);
	xX.Xnames = xX.Xnames(xX.iC);
	xX.X      = [xX.X kron(ones(nsess,1),eye(k))];
	xX.iB     = [1:k] + l;
	for   i = 1:k
		X.Xnames{l + i} = sprintf('scan: %i ',i);
	end
	DSstr = '[burst-mode]';
end


% Temporal filtering
%=======================================================================

% High-pass filtering
%-----------------------------------------------------------------------
if BM
	cLF = 'none';
else
	cLF = spm_input('High-pass filter?','+1','b','none|specify');

end
switch cLF

	case 'specify'

	% default based on peristimulus time
	% param = cut-off period (max = 512, min = 32)
	%-------------------------------------------------------------------
	HParam = 512*ones(1,nsess);
	for  i = 1:nsess
		for j = 1:length(Sess{i}.pst)
			HParam(i) = min([HParam(i) 2*max(RT + Sess{i}.pst{j})]);
		end
	end
	HParam = ceil(HParam);
	HParam(HParam < 32) = 32;
	str    = 'session cutoff period (secs)';
	HParam = spm_input(str,'+1','e',HParam,[1 nsess]);

	% LF description
	%-------------------------------------------------------------------
	LFstr = sprintf('[min] Cutoff period %d seconds',min(HParam));

	case 'none'
	%-------------------------------------------------------------------
	HParam = cell(1,nsess);
	LFstr  = cLF;

end


% Low-pass filtering
%-----------------------------------------------------------------------
if BM
	cHF = 'none';
else
	cHF = spm_input('Low-pass filter?','+1','none|Gaussian|hrf');

end
switch cHF

	case 'Gaussian'
	%-------------------------------------------------------------------
	LParam  = spm_input('Gaussian FWHM (secs)','+1','r',4);
	HFstr   = sprintf('Gaussian FWHM %0.1f seconds',LParam);
	LParam  = LParam/sqrt(8*log(2));

	case {'hrf', 'none'}
	%-------------------------------------------------------------------
	LParam  = [];
	HFstr   = cHF;

end

% create filter struct and band-pass specification
%-----------------------------------------------------------------------
for i = 1:nsess
	K{i} = struct(	'HChoice',	cLF,...
			'HParam',	HParam(i),...
			'LChoice',	cHF,...
			'LParam',	LParam,...
			'row',		row{i},...
			'RT',		RT);
end


% intrinsic autocorrelations (Vi)
%-----------------------------------------------------------------------
str     = 'Model intrinsic correlations?';
cVimenu = {'none','AR(1)'};
cVi     = spm_input(str,'+1','b',cVimenu);


% the interactive parts of spm_spm_ui are now finished: Cleanup GUI
%-----------------------------------------------------------------------
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% Contruct K and Vi structs
%=======================================================================
K       = spm_filter('set',K);

% create Vi struct
%-----------------------------------------------------------------------
Vi      = speye(sum(nscan));
xVi     = struct('Vi',Vi,'Form',cVi);
for   i = 1:nsess
	xVi.row{i} = row{i};
end


% get file identifiers and Global values
%=======================================================================
fprintf('%-40s: ','Mapping files')                                   %-#
VY     = spm_vol(P);
fprintf('%30s\n','...done')                                          %-#

if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0])
	error('images do not all have the same dimensions'),           end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end


%-Compute Global variate
%-----------------------------------------------------------------------
GM     = 100;
q      = sum(nscan);
g      = zeros(q,1);
fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
for i  = 1:q
	fprintf('%s%30s',sprintf('\b')*ones(1,30),sprintf('%3d/%-3d',i,q)) %-#
	g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


% scale if specified (otherwise session specific grand mean scaling)
%-----------------------------------------------------------------------
gSF     = GM./g;
if strcmp(Global,'None')
	for i = 1:nsess
		j      = row{i};
		gSF(j) = GM./mean(g(j));
	end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for  i = 1:q, VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); end


%-Masking structure
%-----------------------------------------------------------------------
xM     = struct('T',	ones(q,1),...
		'TH',	g.*gSF,...
		'I',	0,...
		'VM',	{[]},...
		'xs',	struct('Masking','analysis threshold'));


%-Complete design matrix (xX)
%=======================================================================
xX.K   = K;
xX.xVi = xVi;



%-Effects designated "of interest" - constuct an F-contrast
%-----------------------------------------------------------------------
F_iX0 = xX.iB;


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
for i    = 1:length(Sess), ntr(i) = length(Sess{i}.name); end
sGXcalc  = 'mean voxel value';
sGMsca   = 'session specific';
xsDes    = struct(	'Design',			DSstr,...
			'Basis_functions',		BFstr,...
			'Number_of_sessions',		sprintf('%d',nsess),...
			'Conditions_per_session',	sprintf('%-3d',ntr),...
			'Interscan_interval',		sprintf('%0.2f',RT),...
			'High_pass_Filter',		LFstr,...
			'Low_pass_Filter',		HFstr,...
			'Intrinsic_correlations',	xVi.Form,...
			'Global_calculation',		sGXcalc,...
			'Grand_mean_scaling',		sGMsca,...
			'Global_normalisation',		Global);
%-global structure
%-----------------------------------------------------------------------
xGX.iGXcalc  = Global{1};
xGX.sGXcalc  = sGXcalc;
xGX.rg       = g;
xGX.sGMsca   = sGMsca;
xGX.GM       = GM;
xGX.gSF      = gSF;


%-Save SPMcfg.mat file
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPMstats configuration')                   %-#
save SPMcfg xsDes VY xX xM xGX F_iX0 Sess
fprintf('%30s\n','...SPMcfg.mat saved')                              %-#


%-Display Design report
%=======================================================================
spm_DesRep('DesMtx',xX,{VY.fname}',xsDes)


%-Analysis Proper
%=======================================================================
spm_clf(Finter);
if spm_input('estimate?',1,'b','now|later',[1,0],1)
	spm('Pointer','Watch')
	spm('FigName','Stats: estimating...',Finter,CmdLine);
	spm_spm(VY,xX,xM,F_iX0,Sess);
	spm('Pointer','Arrow')
else
	spm('FigName','Stats: configured',Finter,CmdLine);
	spm('Pointer','Arrow')
	% spm_DesRep(****
end


%-End: Cleanup GUI
%-----------------------------------------------------------------------
fprintf('\n\n')



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function abort = sf_abort(i)
%=======================================================================
if nargin<1, i=[1:3]; end
tmp    = zeros(1,3);
tmp(i) = 1;
tmp = tmp & [	exist(fullfile('.','fMRIDesMtx.mat'),'file')==2 ,...
		exist(fullfile('.','SPMcfg.mat'),    'file')==2 ,...
		exist(fullfile('.','SPM.mat'),       'file')==2 ];
if any(tmp)
	str = {	'        SPM fMRI design matrix config. (fMRIDesMtx.mat)',...
		'        SPMstats configuration (SPMcfg.mat)',...
		'        SPMstats results files (inc. SPM.mat)'};
	str = {	'Current directory contains existing SPMstats files:',...
		str{tmp},['(pwd = ',pwd,')'],' ',...
		'Continuing will overwrite existing files!'};
	abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename);
	if abort, fprintf('%-40s: %30s\n\n',...
		'Abort...   (existing SPMstats files)',spm('time')), end
else
	abort = 0;
end
