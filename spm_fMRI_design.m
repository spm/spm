function [xX,Sess] = spm_fMRI_design(nscan,RT)
% Assembles a design matrix for fMRI studies
% FORMAT [xX,Sess] = spm_fMRI_design(nscan,RT)
%
% nscan         - n vector {nscan(n) = number of scans in session n}
% RT            - interscan interval {seconds}
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
% Sess{s}.U		- cell of (n) trial-specific structures;
% Sess{s}.bf		- basis functions
% Sess{s}.row		- scan   indices for session s
% Sess{s}.col		- effect indices for session s
% Sess{s}.Fci		- indices X1 for F-contrast {within session s}
% Sess{s}.Fcname	- trial-specific F-contrast names;
% Sess{s}.Bfname	- basis function description string;
% Sess{s}.rep 		- session replication flag       {0 or 1};
% Sess{s}.V		- order of [Volterra] covolution {1 or 2};
%
% Sess{s}.U             - {1 x n} cell of (n) trial-specific structures
%
%       U{i}.Uname - {1 x j} cell of names for each input or cause
% 	U{i}.u     - (k x j) inputs or stimulus function matrix
% 	U{i}.dt    -         time bin (seconds)
% 	U{i}.ons   - (q x 1) onsets for q trials (time bins + 32)
% 	U{i}.off   - (q x 1) offset (time bins + 32)
% 	U{i}.pst   - (k x 1) peristimulus times (seconds)
% 	U{i}.Pname - {1 x p} cell of parameter name
% 	U{i}.P     - (q x p) parameter matrix
% 	U{i}.Pi    - {1 x p} cell of sub-indices of U{i}.u
%
% saves SPM_fMRIDesMtx.mat (xX Sess)
%_______________________________________________________________________
%
% spm_fMRI_design allows you to build design matrices with separable
% session-specific partitions.  Each partition may be the same (in which
% case it is only necessary to specify it once) or different.  Responses
% can be either event- or epoch related, where the latter model prolonged
% and possibly time-varying responses to state-related changes in
% experimental conditions.  Event-related response are modelled in terms
% of responses to instantaneous events.  Mathematically they are both
% modelled by convolving a series of delta (stick) or box-car functions,
% encoding the input or stimulus function. with a set of basis
% functions.  
%
% spm_fMRI_design allows you to combine both event- and epoch-related
% responses in the same model.  You are asked to specify the number
% of trial (event or epoch) types.  Epoch and event-related 
% responses are modeled in exactly the same way by first specifying their
% onsets [in terms of onset times] and then their durations.  Events are 
% specified with a duration or 0.  If you enter a single number for the
% durations it will be assumed that all trials conform to this duration.
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
%
%                           ----------------
%
% spm_get_ons contructs a cell of sparse input functions U{i}.u specifying 
% occurrence events or epochs (or both). These are convolved with a basis set
% at a later stage to give regressors that enter into the design matrix.
% Interactions of evoked responses with some parameter (time or a specified 
% variate P) enter at this stage as additional columns in U{i}.u with each
% trial multiplied by the [expansion of the] trial-specific parameter.
% If parametric modulation is modeled, P contains the original variate and
% Pname is its name.  Otherwise P{i} = [] and Pname{i} = '';
%
%                           ----------------
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
%
%                           ----------------
%
% For first order expansions spm_Volterra simply convolves the causes
% (e.g. stick functions) in U{i}.u by the basis functions in Sess{s}.bf
% to create design matrix X.  For second order expansions new entries appear
% in the deisgn matrix that correspond to the hemodynamic interaction among the
% orginal causes (if the events are sufficiently close in time).
% The basis functions for these are two dimensional and are used to
% assemble the second order kernel in spm_graph.m.  Second order effects
% are computed for only the first column of U{i}.u.
%
%_______________________________________________________________________
% %W%   Karl Friston %E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% This routine uses the BCH gobal variable for spm_input : 
%    BCH.bch_mat 
%    BCH.index0  = {'model',index_of_Analysis};
%_______________________________________________________________________

SCCSid  = '%I%';

%-GUI setup
%-----------------------------------------------------------------------
SPMid = spm('SFnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','fMRI stats model setup',0);
spm_help('!ContextHelp',mfilename)


% construct Design matrix {X} - cycle over sessions
%=======================================================================

% global parameters
%-----------------------------------------------------------------------
global fMRI_T fMRI_T0 UNITS
if isempty(fMRI_T),  fMRI_T  = 16; end;
if isempty(fMRI_T0), fMRI_T0 = 1;  end;

% get nscan and RT if no arguments
%-----------------------------------------------------------------------
if nargin < 2
	spm_input('Basic parameters...',1,'d',mfilename,'batch')
	RT = spm_input('Interscan interval {secs}','+1','r',[],1,...
                       'batch',{},'RT');
end
if nargin < 1
        nscan = spm_input(['scans per session e.g. 64 64 64'],'+1',...
                          'batch',{},'nscans');
end
nsess = length(nscan);


% time units, dt = time bin {secs}
%-----------------------------------------------------------------------
dt    = RT/fMRI_T;					
str   = 'specify design in';
UNITS = spm_input(str,'+1','scans|secs','batch',{},'units');

% separate specifications for non-relicated sessions
%-----------------------------------------------------------------------
rep = 0;
if nsess > 1 & ~any(diff(nscan))
	str = 'are sessions exact replications';
	rep = spm_input(str,'+1','yes|no',[1 0],'batch',{},'replicated');
end
Xx    = [];
Xb    = [];
Xname = {};
Bname = {};
Sess  = {};

% Get basis functions
%-----------------------------------------------------------------------
[bf Bfn] = spm_get_bf(dt);

% 1st or 2nd order Volterra expansion?
%-----------------------------------------------------------------------
str   = 'model interactions (Volterra)';
V     = spm_input(str,'+1','y/n',[2 1],'batch','volterra');


% get session specific design parameters
%=======================================================================
for s = 1:nsess

	% set prompt string
	%---------------------------------------------------------------
	if rep
		Fstr = 'All sessions';
	else
		Fstr = sprintf('Session %d/%d',s,nsess);
	end

	k     = nscan(s);
	if (s == 1) | ~rep


		% create convolved stimulus functions
		%=======================================================

		% Inputs, neuronal casues or stimulus functions U
		%-------------------------------------------------------
		U              = spm_get_ons(k,fMRI_T,dt,s);


		% Convolve stimulus functions with basis functions
		%-------------------------------------------------------
		[X,Xn,Fci,Fcn] = spm_Volterra(U,bf,V);


		% Resample regressors at acquisition times (32 bin offset)
		%-------------------------------------------------------
		if size(X,2)
			X      = X([0:(k - 1)]*fMRI_T + fMRI_T0 + 32,:);
		end


		% get user specified regressors
		%=======================================================
		spm_input('Other regressors',1,'d',Fstr,'batch')
		D       = [];
		c       = spm_input('user specified regressors','+1','w1',0,...
    	                          'batch',{},'regressors_nb',s);
                while size(D,2) < c
		      str = sprintf('regressor %i',size(D,2) + 1);
		      D   = [D spm_input(str,2,'e',[],[k Inf],...
                                    'batch',{'regressors',s},'values')];
		end

		% append regressors and names
		%-------------------------------------------------------
		for i = 1:size(D,2)
			X           = [X D(:,i)];
			str         = sprintf('regressor %i',i);
			Xn{end + 1} = spm_input(['name of ' str],...
			'+0','s',str,'batch',{'regressors',s},'names',i);
		end


		% Confounds: Session effects 
		%=======================================================
		B      = ones(k,1);
		Bn{1}  = sprintf('constant');

	end

	% Session structure
	%---------------------------------------------------------------
	Sess{s}.U      = U;
	Sess{s}.bf     = bf;
	Sess{s}.row    = size(Xx,1) + [1:k];;
	Sess{s}.col    = size(Xx,2) + [1:size(X,2)];;
	Sess{s}.Fci    = Fci;
	Sess{s}.Fcname = Fcn;
	Sess{s}.Bfname = Bfn;
	Sess{s}.rep    = rep;
	Sess{s}.V      = V;

	% Append names
	%---------------------------------------------------------------
	for i = 1:length(Xn) 
		Xname{end + 1} = [sprintf('Sn(%i) ',s) Xn{i}];
	end
	for i = 1:length(Bn) 
		Bname{end + 1} = [sprintf('Sn(%i) ',s) Bn{i}];
	end

	% append into Xx and Xb
	%===============================================================
	Xx    = blkdiag(Xx,X);
	Xb    = blkdiag(Xb,B);

end %- for s = 1:nsess

% finished
%-----------------------------------------------------------------------
xX     = struct(	'X',		[Xx Xb],...
			'RT',		RT,...
			'iH',		[],...
			'iC',		[1:size(Xx,2)],...
			'iB',		[1:size(Xb,2)] + size(Xx,2),...
			'iG',		[],...
			'Xnames',	{[Xname Bname]});


%-End: Save SPM_fMRIDesMtx.mat
%-----------------------------------------------------------------------
fprintf('\t%-32s: ','Saving fMRI design')                            %-#
save('SPM_fMRIDesMtx','SPMid','xX','Sess');
fprintf('%30s\n\n','...SPM_fMRIDesMtx.mat saved')                    %-#
spm_input('!DeleteInputObj')

