function [sf,Cname,Pv,Pname,DSstr] = spm_get_ons(k,T,dt,STOC,Fstr,v,Cname,s)
% returns onset times for events
% FORMAT [sf,Cname,Pv,Pname,DSstr] = spm_get_ons(k,T,dt,STOC,Fstr,v,Cname,s)
%
% k     - number of scans
% T     - time bins per scan
% dt    - time bin length (secs)
% STOC  - flag to enable stochastic designs [0 or 1]
% Fstr  - Prompt string (usually indicates session)
% v     - number of conditions or trials 		: can be empty
% Cname - {1 x v}   cell of names for each condition 	: can be empty
% s	- session number (used for batch system)
%
% sf    - {1 x n}   cell of stick function matrices
% Cname - {1 x n}   cell of names for each condition
% Pv    - {1 x n}   cell of parametric vectors
% Pname - {1 x n}   cell of names for each parameter
% DSstr - Design string
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
%
% Notes on responding to questions:
%
% 'number of conditions or trials':  The number of conditions, trials,
%        events or epochs in the design.  Generally the baseline condition
%        (epoch-related) or null event (event-related) should not be included
%        e.g. for a simple ABABAB.. design enter 1
% 
% STOCHASTIC DESIGNS
%
% 'stochastic design': If you want a random design select yes.  The ensuing
%        design matrix and onset times in Sess are then used in 
%        subsequent analysis of the data and stimulus design respectively.
%
%       'include a null event': for stochastic designs a null event should
%                be included if you want to estimate responses common to
%                all trial types
%
%       'SOA (scans)': Stimulus onset asynchrony for the sucessive occurrence
%                of trials.  This is the time (in scans) between the onset
%                of sucessive stimuli or trials (usually a fraction of a scan)
%
%       'relative frequency [trial 1,..n null]':  Enter a vector with a
%                relative frequency of presentation for each trial type
%                (and the null event if included).  The null event is last.
%                The most efficient designs are given when all the frequencies
%                are equal.
%
%       'stationary|modulated': If the occurence probabilities are
%                the same for all scans then choose 'stationary'.  Modulated
%                designs are more efficient but entail 'runs' of the
%                same trial type.
%
% NON STOCHASTIC DESIGNS
%
% 'Fixed|Variable':  If the event of epoch starts with a fixed
%        SOA choose 'Fixed'. If the SOA changes within any trial type
%        choose variable.
%
%        'vector of onsets (scans) for trial n':  If the SOA are variable
%                you have to enter a vector of onet times for each event or
%                epoch.  Time is specified in terms of scans, where the
%                start of the session begins at 0.
% 
%        'SOA (scans)' and 'first trial (scans)':  If the SOA is fixed you
%                only have to specify what it is and when the first condition 
%                starts. 
%
% 'parametric modulation':  This allows you to model time of other effects
%         on eveoked responses in terms of an interaction with the specified
%         variate.
%
% SLCIE TIMIING
%
% With longs TRs you may want to shift the regressors so that they are
% aligned to a particular slice.  This is effected by resetting the
% values of fMRI_T and fMRI_T0 in som_defaults.  fMRI_T is the number of
% time-bins per scan used when building regressors.  Onsets are defined
% in temporal units of scans starting at 0.  fMRI_T0 is the first
% time-bin at which the regressors are resampled to coincide with data
% acquisition.  If fMRI_T0 = 1 then the regressors will be appropriate
% for the first slice.  If you want to temporally realign the regressors
% so that they match responses in the middle slice then make fMRI_T0 =
% fMRI_T/2 (assuming there is a negligible gap between volume
% acquisitions. Default values are fMRI_T = 16 and fMRI_T0 = 1.
%
%
%_______________________________________________________________________
% @(#)spm_get_ons.m	2.19 Karl Friston 99/08/13

global batch_mat;
global iA;

%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin < 5, Fstr = ''; end
spm_input(Fstr,1,'d','batch',batch_mat)

% initialize variables
%-----------------------------------------------------------------------
sf     = {};
Pv     = {};
Pname  = {};
DSstr  = '';

% get stick functions {ons} and names
%=======================================================================

% get trials
%-----------------------------------------------------------------------
if isempty(v)
	v     = spm_input('number of conditions or trials',2,'w1',...
                    'batch',batch_mat,{'model',iA,'conditions',s},'number');
end
if isempty(Cname)
	Cname = {};
	for i = 1:v
		% get names
		%---------------------------------------------------------------
		str         = sprintf('name for condition/trial %d ?',i);
		Cname{i}    = spm_input(str,3,'s',sprintf('trial %d',i),...
                      'batch',batch_mat,{'model',iA,'conditions',s},'names',i);
	end
end


% event/epoch-related responses
%-----------------------------------------------------------------------
if v    

	% stochastic designs
	%---------------------------------------------------------------
	spm_input('Trial specification...',1,'d',Fstr,'batch',batch_mat)
	if STOC
		 STOC = spm_input('stochastic design','+1','y/n',[1 0],...
                       'batch',batch_mat,{'model',iA,'stochastics',s},'specify');
	end
	if STOC

		% minimum SOA
		%-------------------------------------------------------
		ne      = spm_input('include a null event','+1','y/n',[1 0],....
                         'batch',batch_mat,{'model',iA,'stochastics',s},'null_event');
		soa     = spm_input('SOA (scans)','+1','r',2,...
                         'batch',batch_mat,{'model',iA,'stochastics',s},'soa')*T;
		on      = fix(1:soa:(k*T));
		ns      = length(on);
		DSstr   = [DSstr sprintf('Stochastic: %.2fsec SOA ',soa*dt)];

		% occurence probabilities - stationary
		%-------------------------------------------------------
		if ne
		    str = sprintf('relative frequency [trial 1,..%d null]',v);
		else
		    str = sprintf('relative frequency [trial 1,..%d]',v);
		end
		P       = ones(1,(v + ne));
		P       = spm_input(str,'+1','r',P,[1 (v + ne)],...
                          'batch',batch_mat,...
			  {'model',iA,'stochastics',s},'relative_frequency');
		str     = 'occurence probability';
		if spm_input(str,'+1','stationary|modulated',[1 0], ...
			'batch', batch_mat, ...
			{'model',iA,'stochastics',s},'stationary_or_modulated')		
			DSstr = [DSstr '(stationary) '];
			P     = P(:)*ones(1,ns);
 
		% occurence probabilities - modulated (32 sec period)
		%-------------------------------------------------------
		else
			DSstr = [DSstr '(modulated) '];
			p     = ones((v + ne),ns);
			dc    = 32/dt;
			for i = 1:(v + ne);
				q      = sin(2*pi*(on/dc + (i - 1)/(v + ne)));
				p(i,:) = 1 + q;
			end
			P     = diag(P)*p;
		end

		% assign trials
		%-------------------------------------------------------
		P     = [zeros(1,ns); cumsum(P)];
		P     = P*diag(1./max(P));
		q     = zeros(size(on));
		Q     = rand(size(on));
		for i = 1:(v + ne);
			j       = find(Q >= P(i,:) & Q < P(i + 1,:));
			q(j)    = i;
		end

		% create stick functions
		%-------------------------------------------------------
		ons   = sparse(on,q,1,k*T,v + ne);

		% stick function array (and delete null event)
		%-------------------------------------------------------
		for  i = 1:v
			sf{i}   = full(ons(:,i));
		end

	% non-stochastic designs
	%---------------------------------------------------------------
	else

	    % get onsets
	    %-----------------------------------------------------------
	    Sstr  = spm_input('SOA',2,'Fixed|Variable',...
                 'batch',batch_mat,{'model',iA,'conditions',s},'fix_var_SOA');
	    DSstr = [DSstr  Sstr ' SOA '];
	    i     = 0;
	    while i < v

		% get onsets
		%-------------------------------------------------------
		switch Sstr
			%- In batch mode, Sstr is always 'Variable'
			case 'Fixed'
			%-----------------------------------------------
			str   = ['SOA (scans) for ' Cname{i + 1}];
			soa   = spm_input(str,3,'r');
			on    = spm_input('time to first trial (scans)',4,'r',0);
			on    = on:soa:k;

			case 'Variable'
			%-----------------------------------------------
			str   = ['vector of onsets (scans) for ' Cname{i + 1}];
			on    = spm_input(str,3,'batch',batch_mat,...
                     {'model',iA,'conditions',s},'onsets',i+1);
		end

		if iscell(on)

			% create stick functions
			%-----------------------------------------------
			for j = 1:length(on)
				i     = i + 1;
	    			ons   = sparse(k*T,1);
				ons(round(on{j}*T + 1)) = 1;
				sf{i} = ons(1:(k*T));
			end
		else
			% create stick functions
			%-----------------------------------------------
			i     = i + 1;
	    		ons   = sparse(k*T,1);
			ons(round(on*T + 1)) = 1;
			sf{i} = ons(1:(k*T));
		end
	    end
	end


	% get parameters, contruct interactions and append
	%================================================================
	spm_input('Parametric specification...','+1','d',Fstr,'batch',batch_mat)

	% paramteric representation of causes - defaults for main effects
	%----------------------------------------------------------------
	for i = 1:v
		Pv{i}     = [];
		Pname{i} = '';
	end

	% get parameter type
	%----------------------------------------------------------------
	Ptype = {'none',...
		 'time',...
		 'other'};
	Ptype = spm_input('parametric modulation','+1','b',Ptype,...
                'batch',batch_mat,{'model',iA,'parametrics',s},'type');
	switch Ptype

		case 'none'
		%--------------------------------------------------------
		return

		case 'other'
		%--------------------------------------------------------
		Pstr   = spm_input('name of parameter','+1','s',...
                        'batch',batch_mat,{'model',iA,'parametrics',s},'name');

		case 'time'
		%--------------------------------------------------------
		Pstr   = Ptype;
	end

	% get parameters of expansion
	%----------------------------------------------------------------
	Etype = {'linear',...
		 'exponen',...
		 'polynom'};
	Etype = spm_input('expansion','+1','b',Etype,'batch',batch_mat,...
	                    {'model',iA,'parametrics',s},'exp_type');
	DSstr = [DSstr  '[ x ' Pstr ' (' Etype ')] '];
	switch Etype

		case 'exponen'
		%--------------------------------------------------------
		if strcmp(Ptype,'time')
			h = round(k*T*dt/4);
	           	h = spm_input('time constant {secs}','+1','r',h,...
                      		'batch',batch_mat,...
				{'model',iA,'parametrics',s},'time_cst');

		else
		   h = spm_input('decay constant','+1','r',...
                      'batch',batch_mat,{'model',iA,'parametrics',s},'decay_cst');
		end

		case 'polynom'
		%--------------------------------------------------------
		str       = 'order of polynomial expansion';
		h         = spm_input(str,'+1','r',2,'batch',batch_mat,...
                                       {'model',iA,'parametrics',s},'order');

	end


	% cycle over selected trial types
	%----------------------------------------------------------------
	str   = sprintf('which trial[s] 1 to %d',v);
	Ypos  = spm_input('!NextPos','batch',batch_mat);

	for i = spm_input(str,'+1','e',1,'batch',batch_mat,...
                            {'model',iA,'parametrics',s},'trials')
	   	spm_input(Cname{i},Ypos,'d',Fstr,'batch',batch_mat);
		on    = find(sf{i}(:,1));
		ns    = length(on);

		% get parameters
		%-------------------------------------------------------
		switch Ptype

			case 'other'
			%-----------------------------------------------
			str   = ['parameters for ' Cname{i}];
			p     = spm_input(str,'+1','r',[],[ns,1],...
				    'batch',batch_mat,...
                                    {'model',iA,'parametrics',s},'parameters',i);

			case 'time'
			%-----------------------------------------------
			p     = on*dt;

		end

		% expansion
		%--------------------------------------------------------
		switch Etype


			case 'polynom'
			%------------------------------------------------
			u              = spm_detrend(p(:));
			v              = zeros(size(u,1),h + 1);
			q              = sparse(size(sf{i},1),h);
			for j = 0:h
 				v(:,(j + 1)) = (u.^j) - v*(pinv(v)*(u.^j));
			end
			for j = 1:h
				u      = v(:,(j + 1));
				q(:,j) = sparse(on,1,u,size(sf{i},1),1);
			end

			case 'exponen'
			%------------------------------------------------
			q              = exp(-p/h);
			q              = spm_detrend(q(:));
			q              = sparse(on,1,q,size(sf{i},1),1);

			case 'linear'
			%------------------------------------------------
			q              = spm_detrend(p(:));
			q              = sparse(on,1,q,size(sf{i},1),1);


		end

		% append as modulated stick functions
		%--------------------------------------------------------
		sf{i}    = [sf{i} q];
		Pv{i}    = p;
		Pname{i} = Pstr;

	end
end
