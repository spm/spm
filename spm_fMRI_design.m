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
% Sess{s}.BFstr   - basis function description string
% Sess{s}.DSstr   - Design description string
% Sess{s}.row     - scan   indices      for session s
% Sess{s}.col     - effect indices      for session s
% Sess{s}.name{i} - of ith trial type   for session s
% Sess{s}.ind{i}  - column indices      for ith trial type {within session}
% Sess{s}.bf{i}   - basis functions     for ith trial type
% Sess{s}.ons{i}  - stimuli onset times for ith trial type (secs)
% Sess{s}.pst{i}  - peristimulus times  for ith trial type (secs)
% Sess{s}.para{i} - vector of paramters for ith trial type
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
%
% see spm_fmri_spm_ui for details
%
% NB:
%
% spm_fMRI_design allows you to combine both event- and epoch-related
% responses in the same model.  You are asked to specify the number
% of event and epoch types. Enter 0 to skip either (or both if you only
% want to use regressors you have designed outside spm_fMRI_design).
% 
% The design matrix is assembled on a much finer time scale (X.dt) than the 
% TR and is then subsampled at the acquisition times.
%
% Sess{s}.ons{i} contains stimulus onset times in seconds relative to the
% timing of the first scan and are provided to contruct stimuli.
%
% Parametric variations in evoked responses are modelled by modulating
% event or epoch specific delta functions with a spcified vector of
% parameters.  This could be trial-specific reaction times or another
% perfomance or stimulus attribute.  Time effects are modelled at this
% level.
%
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
	RT     = spm_input('Interscan interval {secs}','+1');
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
	qx      = 1;
	Xn      = {};
	Bn      = {};
	PST     = {};
	ONS     = {};
	BF      = {};
	IND     = {};

	% Event/epoch related responses			
	%===================================================================
	k       = nscan(s);

	% event/epoch onsets {ons} and window lengths {W}
	%-------------------------------------------------------------------
	[ons,W,name,para,DSstr] = spm_get_ons(k,fMRI_T,dt,STOC);
	

	% get basis functions for responses
	%-------------------------------------------------------------------
	[BF BFstr] = spm_get_bf(W,dt);


	% peri-stimulus {PST} and onset {ONS} times in seconds
	%-------------------------------------------------------------------
	ntrial  = size(ons,2);
	for   i = 1:ntrial

		on    = find(ons(:,i))*dt;
		pst   = [1:k]*RT - on(1);			
		for j = 1:length(on)
			u      = [1:k]*RT - on(j);
			v      = find(u >= -1);
			pst(v) = u(v);
		end
		ONS{i} = on;
		PST{i} = pst;
	end


	% convolve with basis functions and resample at acquisition times to 
	% create design matrix {X}
	%-------------------------------------------------------------------
	for i = 1:ntrial
		IND{i} = [];
		for  j = 1:size(BF{i},2)
			x      = conv(full(ons(:,i)),BF{i}(:,j));
			X      = [X x([0:k-1]*fMRI_T + fMRI_T0,:)];
			%-- X      = [X x([1:k]*T,:)];
			Xn{qx} = [name{i} sprintf(': bf: %i',j)];
			IND{i} = [IND{i} qx];
			qx     = qx + 1;
		end
	end


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
		X      = [X D(:,i)];
		str    = sprintf('regressor: %i',i);
		Xn{qx} = spm_input(['name of ' str],2,'s',str);
		qx     = qx + 1;
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
	Sess{s}.row   = size(xX,1) + [1:k];
	Sess{s}.col   = size(xX,2) + [1:size(X,2)];
	Sess{s}.name  = name;
	Sess{s}.ind   = IND;
	Sess{s}.bf    = BF;
	Sess{s}.pst   = PST;
	Sess{s}.ons   = ONS;
	Sess{s}.para  = para;

	% Append names
	%-------------------------------------------------------------------
	q     = length(Xname);
	for i = 1:length(Xn) 
		Xname{q + i}  = [sprintf('Sn: %i ',s) Xn{i}];
	end
	q     = length(Bname);
	for i = 1:length(Bn) 
		Bname{q + i}  = [sprintf('Sn: %i ',s) Bn{i}];
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
