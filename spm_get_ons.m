function [U] = spm_get_ons(k,T,dt,s)
% returns input [designed effects] structures
% FORMAT [U] = spm_get_ons(k,T,dt,s)
%
% k  - number of scans
% T  - time bins per scan
% dt - time bin length (secs)
% s  - session number (used by batch system)
%
% U  - {1 x n}   cell of (n) trial-specific structures
% 	U{i}.Uname - cell of names for each input or cause
% 	U{i}.u     - inputs or stimulus function matrix
% 	U{i}.dt    - time bin (seconds)
% 	U{i}.ons   - onsets (time bins + 32)
% 	U{i}.off   - offset (time bins + 32)
% 	U{i}.pst   - peristimulus times (seconds)
% 	U{i}.Pname - cell of parameter name
% 	U{i}.P     - parameter matrix
% 	U{i}.Pi    - cell of sub-indices of u pertaining to columns of P
%_______________________________________________________________________
%
%
% SLICE TIMIING
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
% %W% Karl Friston %E%

% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'model',index_of_Analysis};
%_______________________________________________________________________

global BCH UNITS

%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)


% time units
%-----------------------------------------------------------------------
if ~length(UNITS) 
	UNITS = 'scans';
end
switch UNITS

	case 'scans'
	%----------------------------------------------------------------
	TR = T*dt;

	case 'secs'
	%----------------------------------------------------------------
	TR = 1;
end

%-prompt string
%-----------------------------------------------------------------------
if nargin < 4
	Fstr = ''; 
else
	Fstr = sprintf('Session %d: trial specification in %s',s,UNITS);
end
spm_input(Fstr,1,'d','batch')


% initialize ouput structure
%-----------------------------------------------------------------------
U     = {};

% get stick functions {ons} and names
%=======================================================================

% get trials
%-----------------------------------------------------------------------
v     = spm_input('number of conditions or trials',2,'w1',...
                   'batch',{},'conditions_nb',s);
for i = 1:v

	% initialize
	%---------------------------------------------------------------
	Uname    = {};
	Pname    = {};
	u        = [];
	P        = [];
	Pi       = {};

	% get names
	%---------------------------------------------------------------
	str      = sprintf('name for condition/trial %d ?',i);
	Uname{1} = spm_input(str,3,'s',sprintf('trial %d',i),...
                                       'batch',{'conditions',s},'names',i);

	% get main [trial] effects
	%================================================================

	% onsets
	%---------------------------------------------------------------
	str      = ['vector of onsets - ' Uname{1}];
	ons      = spm_input(str,4,'batch',{'conditions',s},'onsets',i);
	ons      = ons(:);

	% durations
	%---------------------------------------------------------------
	str      = 'duration[s] (events = 0)';
	while 1
		dur = spm_input(str,5,'e','batch',...
                                        {'conditions',s},'durations',i);
		if length(dur) == 1
			dur  = dur*ones(size(ons));
		end
		if length(dur) == length(ons), break, end
		str = sprintf('enter a scalar or [%d] vector',length(ons));

	end
	dur      = dur(:);

	% peri-stimulus times {seconds}
	%---------------------------------------------------------------
	pst    = [1:k]*T*dt - ons(1)*TR;			
	for  j = 1:length(ons)
		w      = [1:k]*T*dt - ons(j)*TR;
		v      = find(w >= -1);
		pst(v) = w(v);
	end



	% add parameters x trial interactions
	%================================================================

	% paramteric representation of causes (u) - 1st = main effects
	%----------------------------------------------------------------
	Ptype    = {'none',...
		    'time',...
		    'other'};
	Ptype    = spm_input('parametric modulation',6,'b',Ptype,...
                          'batch',{},'parametrics_type',s);

	switch Ptype

		case 'none'
		%--------------------------------------------------------
		u        = 1;

		case 'time'
		%--------------------------------------------------------
		Pname{1} = 'time';
		Pi{1}    = [1 2];
		str      = 'time constant for adaptation {secs}';
		p        = ons*TR;
		P        = p;

		% create interaction term
		%--------------------------------------------------------
		h        = round(k*T*dt/4);
		h        = spm_input(str,7,'r',h,...
                      		'batch',{'parametrics',s},'time_cst');
		u        = [p.^0 exp(-p/h)];
		Uname{2} = 'adapation';

		case 'other'
		%--------------------------------------------------------
		str   = ['# parameters (' Uname{1} ')'];
		n     = spm_input(str,7,'n1',1,'batch',{'pnum',s},'pnum');
		u     = ones(length(ons),1);;
		for v = 1:n

			% get names and variates
			%------------------------------------------------
			str      = sprintf('parameter %d name',v);
			Pname{v} = spm_input(str,7,'s',...
				          'batch',{'parametrics',s},'name');

			p     = spm_input(Pname{v},7,'r',[],[length(ons),1],...
				'batch',{'parametrics',s},'parameters');
			p     = p(:);
			P     = [P p];


			% polynomial expansion of u = f(p)
			%------------------------------------------------
			str   = 'order of polynomial expansion';
			h     = spm_input(str,8,'n1',1,...
                                     'batch',{'parametrics',s},'order');

			% sub-indices and inputs
			%------------------------------------------------
			Pi{v} = [1, ([1:h] + size(u,2))];
			for j = 1:h
 				u   = [u p.^j];
				str = sprintf('%sx%s^%d',Uname{1},Pname{v},j);
				Uname{end + 1} = str;
			end

		end

	end

	% orthogonalize inputs
	%---------------------------------------------------------------
	u          = spm_orth(u);

	% and scale so sum(u*dt) = number of events, if event-related 
	%---------------------------------------------------------------
	if ~any(dur)
		u  = u/dt;
	end

	% create stimulus functions (32 bin offset)
	%---------------------------------------------------------------
	sf         = sparse(k*T + 32,size(u,2));
	ons        = round(ons*TR/dt) + 32;			% onsets
	sf(ons,:)  = u;
	off        = round(dur*TR/dt) + ons + 1;		% offsets
	sf(off,:)  = -u;
	sf         = cumsum(sf);				% integrate
	sf         = sf(1:(k*T + 32),:);


	% place in ouputs structure
	%---------------------------------------------------------------
	U{i}.Uname = Uname;			% - input names
	U{i}.u     = sf;			% - stimulus function matrix
 	U{i}.dt    = dt;			% - time bin (seconds)
	U{i}.ons   = ons;			% - onsets (time bins)
	U{i}.off   = off;			% - offset (time bins)
	U{i}.pst   = pst;			% - pst (seconds)
	U{i}.Pname = Pname;			% - parameter name
	U{i}.P     = P;				% - parameter variates
	U{i}.Pi    = Pi;			% - parameter sub-indices

end
