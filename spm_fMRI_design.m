function [X,Sess] = spm_fMRI_design(nscan,RT)
% Assembles a design matrix for fMRI studies
% FORMAT [X,Sess] = spm_fMRI_design(nscan,RT)
%
% nscan   - n vector {nscan(n) = number of scans in session n}
% RT      - intercans interval {seconds}
%
% X.dt    - time bin {secs}
% X.RT    - Repetition time {secs}
% X.xX    - regressors
% X.bX    - session effects
% X.Xname - names of subpartiton columns {1xn}
% X.Bname - names of subpartiton columns {1xn}
%
% Sess{s}.BFstr    - basis function description string
% Sess{s}.DSstr    - Design description string
% Sess{s}.rep      - session replication flag
% Sess{s}.row      - scan   indices      for session s
% Sess{s}.col      - effect indices      for session s
% Sess{s}.name{i}  - of ith trial type   for session s
% Sess{s}.ind{i}   - column indices      for ith trial type {within session}
% Sess{s}.bf{i}    - basis functions     for ith trial type
% Sess{s}.sf{i}    - stick functions     for ith trial type
% Sess{s}.ons{i}   - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}   - peristimulus times  for ith trial type (secs)
% Sess{s}.Pv{i}    - vector of paramters for ith trial type
% Sess{s}.Pname{i} - name   of paramters for ith trial type
%
% saves fMRIDesMtx.mat (X Sess)
%___________________________________________________________________________
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, where the latter model prolonged
% and possibly time-varying responses to state-related changes in
% experimental conditions.  Event-related response are modelled in terms
% of responses to instantaneous events.  Mathematically they are both
% modelled by convolving a series of delta (or stick) functions,
% indicating the onset of an event or epoch with a set of basis
% functions.  These basis functions can be very simple, like a box car,
% or may model voxel-specific forms of evoked responses with a linear
% combination of several basis functions (e.g.  a Fourier set).  Basis
% functions can be used to plot estimated responses to single events or
% epochs once the parameters (i.e.  basis function coefficients) have
% been estimated.  The importance of basis functions is that they provide
% a graceful transition between simple fixed response models (like the
% box-car) and finite impulse response (FIR) models, where there is one
% basis function for each scan following an event or epoch onset.  The
% nice thing about basis functions, compared to FIR models, is that data
% sampling and stimulus presentation does not have to be sychronized
% thereby allowing a uniform and unbiased sampling of peri-stimulus time.
%
% spm_fMRI_design allows you to combine both event- and epoch-related
% responses in the same model.  You are asked to specify the number
% of condition (epoch) or trial (event) types.  Epoch and event-related 
% responses are modeled in exactly the same way by first specifying their
% onsets [in terms of stimulus onset asynchronies (SOAs) or explicit onset 
% times] and then convolving with appropriate basis functions (short ones
% for event-related models and longer ones for epoch-related respones).
% Enter 0 to skip these if you only want to use regressors you have designed 
% outside spm_fMRI_design).
%
% Interactions or response modulations can enter at two levels.  Firstly
% the stick function itself can be modulated by some parametric variate
% (this can be time or some trial-specific variate like reaction time)
% modeling the interaction between the trial and the variate or, secondly
% interactions among the trials themselves can be modeled using a Volterra
% series formulation that accommodates interactions over time (and therefore
% within and between trial types).  The first sort of interaction is
% specified by extra (modulated) stick functions in Sess{s}.sf{i}.  If
% a polynomial expansion of the specified variate is requested there will
% be more than one additional column.  The corresponding name of the
% explanatory variables in X.Xname is Sn(s) trial i(p)[q] for the qth
% order expansion of the variate convolved with the pth basis function
% of the ith trial in the sth session.  If no parametric variate is
% specified the name is simply Sn(s) trial i(p).  Interactions among
% and within trials enter as new trial types but do not have .pst or .ons
% fields.  These interactions can be characterized later, in results, in
% terms of the corresponding second order Volterra Kernels.
%
% The design matrix is assembled on a much finer time scale (X.dt) than the 
% TR and is then subsampled at the acquisition times.
%
% Sess{s}.ons{i} contains stimulus onset times in seconds relative to the
% timing of the first scan and are provided to contruct stimuli when using 
% stochastic designs
%
%
% Notes on spm_get_ons, spm_get_bf and spm_Volterra are included below
% for convenience.
%_______________________________________________________________________
%
% spm_get_ons contructs a cell of sparse delta functions specifying the
% onset of events or epochs (or both). These are convolved with a basis set
% at a later stage to give regressors that enter into the design matrix.
% Interactions of evoked responses with some parameter (time or a specified 
% variate Pv) enter at this stage as additional columns in sf with each delta
% function multiplied by the [expansion of the] trial-specific parameter.
% If parametric modulation is modeled, P contains the original variate and
% Pname is its name.  Otherwise P{i} = [] and Pname{i} = '';
%-----------------------------------------------------------------------
%
%___________________________________________________________________________
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%---------------------------------------------------------------------------
%___________________________________________________________________________
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in SF by the basis functions in BF to create
% a design matrix X.  For second order expansions new entries appear
% in IND, BF and name that correspond to the interaction among the
% orginal causes (if the events are sufficiently close in time).
% The basis functions for these are two dimensional and are used to
% assemble the second order kernel in spm_graph.m.  Second order effects
% are computed for only the first column of SF.
%---------------------------------------------------------------------------
% @(#)spm_fMRI_design.m	2.9 Karl Friston 99/04/20

% construct Design matrix {X} - cycle over sessions
%===========================================================================

% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
STOC   = 0;

% global parameters
%-----------------------------------------------------------------------
global fMRI_T; 
global fMRI_T0; 
if isempty(fMRI_T),  fMRI_T  = 16; end;
if isempty(fMRI_T0), fMRI_T0 = 16; end;

% get nscan and RT if no arguments
%---------------------------------------------------------------------------
if nargin < 2
	RT     = spm_input('Interscan interval {secs}','+1','r',[],1);
end
if nargin < 1
	nscan  = spm_input(['scans per session e.g. 64 64 64'],'+1');
	STOC   = 1;
end
nsess  = length(nscan);
dt     = RT/fMRI_T;					% time bin {secs}


% separate specifications for non-relicated sessions
%--------------------------------------------------------------------------
rep = 0;
if nsess > 1 & ~any(nscan - nscan(1))
	rep = spm_input(['are sessions replicated exactly'],2,'yes|no',[1 0]);
end
xX    = [];
bX    = [];
Xname = {};
Bname = {};
Sess  = {};
for s = 1:nsess

	% specify design for this session?
	%-------------------------------------------------------------------
	set(Finter,'Name',sprintf('Session %0.0f',s));
        if   (s == 1) | ~rep
	X       = [];
	B       = [];
	D       = [];
	PST     = {};
	ONS     = {};
	IND     = {};

	% Event/epoch related responses			
	%===================================================================
	k       = nscan(s);

	% event/epoch onsets {ons} and window lengths {W}
	%-------------------------------------------------------------------
	[SF,name,Pv,Pname,DSstr] = spm_get_ons(k,fMRI_T,dt,STOC);
	ntrial     = size(SF,2);


	% get basis functions for responses
	%-------------------------------------------------------------------
	[BF BFstr] = spm_get_bf(name,fMRI_T,dt);


	% peri-stimulus {PST} and onset {ONS} times in seconds
	%-------------------------------------------------------------------
	for   i = 1:ntrial
		on     = find(SF{i}(:,1))*dt;
		pst    = [1:k]*RT - on(1);			
		for  j = 1:length(on)
			u      = [1:k]*RT - on(j);
			v      = find(u >= -1);
			pst(v) = u(v);
		end
		ONS{i} = on;
		PST{i} = pst;
	end


	% convolve with basis functions
	%-------------------------------------------------------------------
	str     = 'interactions among trials (Volterra)';
	if spm_input(str,'+1','y/n',[1 0])

		[X Xn IND BF name] = spm_Volterra(SF,BF,name,2);
	else
		[X Xn IND]         = spm_Volterra(SF,BF,name,1);
	end


	% Resample design matrix {X} at acquisition times
	%-------------------------------------------------------------------
	X     = X([0:k-1]*fMRI_T + fMRI_T0,:);


	% get user specified regressors
	%===================================================================
	c     = spm_input('user specified regressors',1,'w1',0);
	while size(D,2) < c
		str   = sprintf('regressor %i',size(D,2) + 1);
		D     = [D spm_input(str,2,'e',[],[k Inf])];
	end
	if      c & length(DSstr)
		DSstr = [DSstr '& user specified covariates '];
	elseif  c
		DSstr = 'User specified covariates ';
	end

	% append regressors and names
	%-------------------------------------------------------------------
	for i = 1:size(D,2)
		X           = [X D(:,i)];
		str         = sprintf('regressor: %i',i);
		Xn{end + 1} = spm_input(['name of ' str],2,'s',str);
	end


	% Confounds: Session effects 
	%===================================================================
	B      = ones(k,1);
	Bn{1}  = sprintf('constant');

	% end specification
	%-------------------------------------------------------------------
	end

	% Session structure
	%-------------------------------------------------------------------
	Sess{s}.BFstr = BFstr;
	Sess{s}.DSstr = DSstr;
	Sess{s}.rep   = rep;
	Sess{s}.row   = size(xX,1) + [1:k];
	Sess{s}.col   = size(xX,2) + [1:size(X,2)];
	Sess{s}.name  = name;
	Sess{s}.ind   = IND;
	Sess{s}.bf    = BF;
	Sess{s}.sf    = SF;
	Sess{s}.pst   = PST;
	Sess{s}.ons   = ONS;
	Sess{s}.Pv    = Pv;
	Sess{s}.Pname = Pname;

	% Append names
	%-------------------------------------------------------------------
	q     = length(Xname);
	for i = 1:length(Xn) 
		Xname{q + i}  = [sprintf('Sn(%i) ',s) Xn{i}];
	end
	q     = length(Bname);
	for i = 1:length(Bn) 
		Bname{q + i}  = [sprintf('Sn(%i) ',s) Bn{i}];
	end

	% append into xX and bX
	%===================================================================
	[x y]   = size(xX);
	[i j]   = size(X);
	xX(([1:i] + x),([1:j] + y)) = spm_detrend(X);
	[x y]   = size(bX);
	[i j]   = size(B);
	bX(([1:i] + x),([1:j] + y)) = B;
end

% finished
%---------------------------------------------------------------------------
X.dt    = dt;
X.RT    = RT;
X.xX    = xX;
X.bX    = bX;
X.Xname = Xname;
X.Bname = Bname;

%-End: Cleanup GUI
%---------------------------------------------------------------------------
save fMRIDesMtx X Sess
spm_clf(Finter);
