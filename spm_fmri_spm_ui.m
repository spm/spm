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
% Sess{s}      -  see spm_fMRI_design
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
% design matrix in terms of events or epochs.  
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, The only distinction is the duration
% of the underlying input or stimulus function. Mathematically they are both
% modeled by convolving a series of delta (stick) or box functions (u),
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions model the hemodynamic convolution,
% applied by the brain, to the inputs.  This convolution can be first-order
% or a generalized convolution modeled to second order (if you specify the 
% Volterra option). [The same inputs are used by the hemodynamic model or
% or dynamic causal models which model the convolution explicitly in terms of 
% hidden state variables (see spm_hdm_ui and spm_dcm_ui).]
% Basis functions can be used to plot estimated responses to single events 
% once the parameters (i.e. basis function coefficients) have
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
% obtain when the probabilities of every trial type are equal.
% A critical issue in stochastic designs is whether to include null events
% If you wish to estimate the evoke response to a specific event
% type (as opposed to differential responses) then a null event must be
% included (even if it is not modelled explicitly).
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
% inference.
% 
% See spm_fMRI_design for more details about how designs are specified.
%
% Serial correlations in fast fMRI time-series are dealt with as
% described in spm_spm.  At this stage you need to specify the filtering
% that will be applied to the data (and design matrix) to give a
% generalized least squares (GLS) estimate of the parameters required.
% This filtering is important to ensure that the GLS estimate is
% efficient and that the error variance is estimated in an unbiased way.
% 
% The serial correlations will be estimated with a ReML (restricted
% maximum likelihood) algorithm using an AR(1) plus white noise model
% during parameter estimation.  This estimate assumes the same
% correlation structure for each voxel, within each session.  The ReML
% estimates are then used to correct for non-sphericity during inference
% by adjusting the statistics and degrees of freedom appropriately.  The
% discrepancy between estimated and actual intrinsic (i.e. prior to
% filtering) correlations are greatest at low frequencies.  Therefore
% specification of the high-pass filter is particularly important.
% Furthermore, high-pass filtering whitens the data and renders parameter
% estimation more efficient (by Gauss-Markov theorum).  Low-pass
% filtering (i.e. smoothing) has been removed as an option because the
% ReML estimator of serial correlations resolves the robustenss issues to
% a degree. High-pass filtering is implemented at the level of the
% filtering matrix K (as opposed to entering as confounds in the design
% matrix).  The default cutoff period is 128 seconds.  Use 'explore design'
% to ensure this cuttof is not removing too much experimental variance.
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
%_______________________________________________________________________
% %W% Karl Friston, Jean-Baptiste Poline, Christian Buchel %E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'model',index_of_Analysis};
%_______________________________________________________________________

SCCSid  = '%I%';

global BCH; %- used as a flag to know if we are in batch mode or not.

%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);
spm_help('!ContextHelp',mfilename)


% get design matrix and/or data
%=======================================================================
MType = {  'specify a model',...
	   'review a specified model',...
	   'estimate a specified model',...
	   'specify and estimate a model'};
MT    = spm_input('What would you like to do?',1,'m',MType,...
                  'batch',{},'types');

%-Initialise output arguments
%-----------------------------------------------------------------------
xX   = [];
Sess = {};

switch MT
%-----------------------------------------------------------------------

	case 1
	% specify a design matrix
	%---------------------------------------------------------------
	if sf_abort, spm_clf(Finter), return, end
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
	if isempty(BCH)
	   load(spm_get(1,'fMRIDesMtx.mat','Select SPM_fMRIDesMtx.mat'));
	else
	   load('SPM_fMRIDesMtx.mat');
	end


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
		   str    = sprintf('select scans for session %0.0f',i);
		   if isempty(BCH)
			q = spm_get(nscan(i),'.img',str);
		   else
			q = sf_bch_get_q(i);
		   end %- 
 		   P      = strvcat(P,q);
		end
	else
		str       = sprintf('select scans for this study');
		if isempty(BCH)
			P = spm_get(sum(nscan),'.img',str);
		else
		   for  i = 1:nsess
			q = sf_bch_get_q(i);
			P = strvcat(P,q);
		   end
		end
	end


	% Repeat time
	%---------------------------------------------------------------
	RT     = xX.RT;

	case 4
	% get filenames and design matrix
	%---------------------------------------------------------------
	if sf_abort, spm_clf(Finter), return, end
	spm_input('Scans & sessions...',1,'d',mfilename,'batch')
	nsess  = spm_input(['number of sessions'],'+1','e',1,...
                'batch',{},'nsess');
	nscan  = zeros(1,nsess);
	P      = [];
	for  i = 1:nsess
		str  = sprintf('select scans for session %0.0f',i);
		if isempty(BCH)
		   q = spm_get(Inf,'.img',str);
		else
		   q = sf_bch_get_q(i);
		end
 		P        = strvcat(P,q);
		nscan(i) = size(q,1);
	end

	% get Repeat time
	%---------------------------------------------------------------
	RT  = spm_input('Interscan interval {secs}','+1','batch',{},'RT');

	% get design matrix
	%---------------------------------------------------------------
	[xX,Sess] = spm_fMRI_design(nscan,RT);

end


% Assemble remaining design parameters
%=======================================================================
spm_help('!ContextHelp',mfilename)
spm_input('Global intensity normalisation...',1,'d',mfilename,'batch')


% Global normalization
%-----------------------------------------------------------------------
str    = 'remove Global effects';
Global = spm_input(str,'+1','scale|none',{'Scaling' 'None'},...
		    'batch',{},'global_effects');

% Temporal filtering
%=======================================================================
spm_input('Temporal autocorrelation options','+1','d',mfilename,'batch')

% High-pass filtering
%-----------------------------------------------------------------------
cF    = spm_input('High-pass filter?','+1','b','none|specify',...
			'batch',{},'HF_fil');
switch cF

	case 'specify'

	% default 128 seconds
	%---------------------------------------------------------------
	HParam = 128*ones(1,nsess);
	str    = 'session cutoff period (secs)';
	HParam = spm_input(str,'+1','e',HParam,[1 nsess],...
	                  'batch',{},'HF_cut');

	% Filter description
	%---------------------------------------------------------------
	Fstr   = sprintf('[min] Cutoff period %d seconds',min(HParam));

	case 'none'
	%---------------------------------------------------------------
	HParam = Inf*ones(1,nsess);;
	Fstr   = cF;
end


% create and set filter struct
%-----------------------------------------------------------------------
for i = 1:nsess
	K{i} = struct(	'HParam',	HParam(i),...
			'row',		Sess{i}.row,...
			'RT',		xX.RT);
end
K       = spm_filter(K);



% intrinsic autocorrelations (Vi)
%-----------------------------------------------------------------------
str     = 'Correct for serial correlations?';
cVi     = {'none','AR(1) + w'};
cVi     = spm_input(str,'+1','b',cVi,'batch',{},'int_corr');


%-Estimate now or later?
%-----------------------------------------------------------------------
bEstNow = spm_input('estimate?','_','b','now|later',[1,0],1,...
		'batch',{},'now_later');

%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm_clf(Finter);
spm('FigName','Configuring, please wait...',Finter,CmdLine);
spm('Pointer','Watch');


% Contruct Vi structure for non-sphericity ReML estimation
%=======================================================================

% create Vi struct
%-----------------------------------------------------------------------
switch cVi

	case 'AR(1) + w'
	%---------------------------------------------------------------
	Q  = spm_Ce(nscan,2);

	case 'none'
	%---------------------------------------------------------------
	Q  = speye(sum(nscan));

end
xVi.Vi   = Q;
xVi.form = cVi;


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
	fprintf('%s%30s',sprintf('\b')*ones(1,30),sprintf('%4d/%-4d',i,q)) %-#
	g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#


% scale if specified (otherwise session specific grand mean scaling)
%-----------------------------------------------------------------------
gSF     = GM./g;
if strcmp(Global,'None')
	for i = 1:nsess
		j      = Sess{i}.row;
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


%-Effects designated "of interest" - constuct F-contrast structure array
%-----------------------------------------------------------------------
if length(xX.iC)
	F_iX0  = struct(	'iX0',		xX.iB,...
				'name',		'effects of interest');
else
	F_iX0  = [];
end

%-Trial-specifc effects specified by Sess
%-----------------------------------------------------------------------
for s = 1:nsess
	str   = sprintf('Session %d: ',s);
	for i = 1:length(Sess{s}.Fci)
		q                  = 1:size(xX.X,2);
		FciX               = Sess{s}.col(Sess{s}.Fci{i});
		q(FciX)            = [];
		F_iX0(end + 1).iX0 = q;
		F_iX0(end).name    = [str Sess{s}.Fcname{i}];
	end
end


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
for i   = 1:length(Sess), ntr(i) = length(Sess{i}.U); end
sGXcalc = 'mean voxel value';
sGMsca  = 'session specific';
xsDes   = struct(	'Basis_functions',	Sess{1}.Bfname,...
			'Number_of_sessions',	sprintf('%d',nsess),...
			'Trials_per_session',	sprintf('%-3d',ntr),...
			'Interscan_interval',	sprintf('%0.2f',xX.RT),...
			'High_pass_Filter',	Fstr,...
			'Serial_correlations',	xVi.form,...
			'Global_calculation',	sGXcalc,...
			'Grand_mean_scaling',	sGMsca,...
			'Global_normalisation',	Global);
%-global structure
%-----------------------------------------------------------------------
xGX.iGXcalc  = Global{:};
xGX.sGXcalc  = sGXcalc;
xGX.rg       = g;
xGX.sGMsca   = sGMsca;
xGX.GM       = GM;
xGX.gSF      = gSF;


%-Save SPMcfg.mat file
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPMstats configuration')                   %-#
save('SPMcfg','SPMid','xsDes','VY','xX','xM','xGX','F_iX0','Sess');
fprintf('%30s\n','...SPMcfg.mat saved')                              %-#

 
%-Display Design report
%=======================================================================
fprintf('%-40s: ','Design reporting')                                %-#
spm_DesRep('DesMtx',xX,reshape({VY.fname},size(VY)),xsDes)
fprintf('%30s\n','...done')                                          %-#


%-Analysis Proper
%=======================================================================
spm_clf(Finter);
spm('FigName','fMRI stats models',Finter,CmdLine);
if bEstNow
	spm('Pointer','Watch')
	spm('FigName','Stats: estimating...',Finter,CmdLine);
	spm_spm(VY,xX,xM,F_iX0,Sess,xsDes);
	spm('Pointer','Arrow')
else
	spm_clf(Finter)
	spm('FigName','Stats: configured',Finter,CmdLine);
	spm('Pointer','Arrow')
	spm_DesRep('DesRepUI',struct(	'xX',		xX,...
					'VY',		VY,...
					'xM',		xM,...
					'F_iX0',	F_iX0,...
					'Sess',		{Sess},...
					'xsDes',	xsDes,...
					'swd',		pwd,...
					'SPMid',	SPMid,...
					'cfg',		'SPMcfg'));
end


%-End: Cleanup GUI
%-----------------------------------------------------------------------
fprintf('\n\n')



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function abort = sf_abort(i)
%=======================================================================
if nargin < 1, i = [1:3]; end
tmp    = zeros(1,3);
tmp(i) = 1;
tmp = tmp & [	exist(fullfile('.','SPM_fMRIDesMtx.mat'),'file')==2 ,...
		exist(fullfile('.','SPMcfg.mat'),        'file')==2 ,...
		exist(fullfile('.','SPM.mat'),           'file')==2 ];
if any(tmp)
	str = {	'    SPM fMRI design matrix definition (SPM_fMRIDesMtx.mat)',...
		'    SPMstats configuration            (SPMcfg.mat)',...
		'    SPMstats results files            (inc. SPM.mat)'};
	str = {	'Current directory contains existing SPMstats files:',...
		str{tmp},['(pwd = ',pwd,')'],' ',...
		'Continuing will overwrite existing files!'};

	abort = spm_input(str,1,'bd','stop|continue',[1,0],1,mfilename,...
                         'batch',{},'stop_writing');
	if abort, fprintf('%-40s: %30s\n\n',...
		'Abort...   (existing SPMstats files)',spm('time')), end
else
	abort = 0;
end

function q = sf_bch_get_q(i)
%=======================================================================
% This is to deal with a specific case where the sampling   
% isn't regular. Only implemented in bch mode.
% 
q 	     = spm_input('batch',{},'files',i);
files 	     = q;
t_sampl      = spm_input('batch',{},'time_sampl',i);
remain       = spm_input('batch',{},'remain',i);
q(remain,:)  = files;

%- fills the gap with the first images ...
%-----------------------------------------------------------------------
%- there should be enough first images!
q(t_sampl,:) = files(1:length(t_sampl),:);
