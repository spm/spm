function [ons,W,name,para] = spm_get_ons(ST,k,T)
% returns onset times for events
% FORMAT [ons,W,name,para] = spm_get_ons(ST,k,T)
% ST   - Study type
% 	1  - 'event-related fMRI (single events)'
% 	2  - 'event-related fMRI (epochs of events)'
% 	3  - 'epoch-related fMRI (fixed length epochs)'
% 
% k    - scans per session
% T    - time bins per scan
%
% ons  - [m x n]   stick functions {onsets for n trial types over m time bins}
% W    - [1 x n]   vector of window/epoch lengths {time-bins}
% name - {1 x n}   cell of names for each trial type
% para - {1 x n/2} cell of parameter vectors {when specified}
%_______________________________________________________________________
% %W% Karl Friston %E%

% time bins per scan
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
name   = {};
para   = {};

% get onsets (ons) and names {names}
%-----------------------------------------------------------------------
if     ST == 1 | ST == 2
 
	% cycle over events
	%---------------------------------------------------------------
	v     = spm_input('number of event or trial types',1,'e',1);
	str   = 'are inter-stimulus intervals';
	ISI   = spm_input(str,2,'fixed|variable',[1 0]);
	ons   = zeros(k*T,v);
	for i = 1:v
	
		% get name
		%-------------------------------------------------------
		str         = sprintf('trial %d',i);
		name{i}     = spm_input('event or trial name',2,'s',str);
		set(Finter,'Name',[name{i} ': event ' sprintf('%i',i)])


		% get onsets
		%-------------------------------------------------------
		if ISI
			isi = spm_input('interstimulus interval (scans)',3);
			on  = spm_input('first stimulus (scans)',4,'e',1);
			on  = on:isi:k;
		else
			on  = spm_input('stimulus onset times (scans)',3);
		end
		on    = round(on*T + 1);
		nons  = length(on);


		% get durations if 'epochs of events'
		%-------------------------------------------------------
		d     = ones(1,nons);
		if ST == 2
			d     = [];
			while length(d) ~= nons
				d   = spm_input('durations[s] (scans)',4);
				d   = round(d*T + 1);
				if length(d) == 1
					d    = ones(1,nons)*d;
				end
			end
		end


		% cumulate
		%-------------------------------------------------------
		for j = 1:length(on)
			q        = [0:(d(j) - 1)] + on(j);
			ons(q,i) = ones(d(j),1);
		end

	end

	% basis function lengths (not used subsequently)
	%---------------------------------------------------------------
	W     = zeros(1,v);


elseif ST == 3


	% vector of conditions
	%---------------------------------------------------------------
	str   = 'epoch order eg 1 2 1...  {0 = null}';
	a     = spm_input(str,1);
	v     = max(a);


	% get names
	%-------------------------------------------------------
	for i = 1:v
		str     = sprintf('name of epoch or condition %d',i);
		name{i} = spm_input(str,2,'s');
	end


	% vector of epoch lengths
	%---------------------------------------------------------------
	W     = zeros(1,v);
	while sum(W(a)) ~= k
		str = sprintf('scans per epoch 1 to %d',v);
		W   = spm_input(str,2);
		while length(W) < v, W = [W W]; end
		W   = W(1:v);
	end
	W     = W*T;


	% epoch onsets
	%---------------------------------------------------------------
	lag   = spm_input('start of first epoch {scans}',3,'e',0);
	c     = cumsum(W(a)) - W(a) + lag*T;
	ons   = zeros(k*T,v);
	for i = 1:v
		on       = c(find(a == i));
		j        = round(on + 1);
		ons(j,i) = ones(length(j),1);
	end
end


% clf
%-----------------------------------------------------------------------
spm_clf(Finter)

% get parameters if 'parametric', contruct interactions and pre-pend
%-----------------------------------------------------------------------
if spm_input('model parametric effects',1,'b','no|yes',[0 1])

	% cycle over trial types
	%---------------------------------------------------------------
	oxp   = zeros(k*T,v);
	pname = {};
	for i = 1:v

		% get delta functions
		%-------------------------------------------------------
		on     = find(ons(:,i));
		nons   = length(on);
		p      = [];
		while length(p) ~= nons
			str = [sprintf('[%d]-parameter for ',nons) name{i}];
			p   = spm_input(str,2);
			p   = p(:);
		end
		para{i}     = spm_detrend(p);

		% fill in modulated deltas
		%-------------------------------------------------------
		oxp(on,i)   = para{i};
		pname{i}    = [name{i} ' x p'];
	end

	% append to ons, names and W
	%---------------------------------------------------------------
	ons   = [oxp ons];
	name  = [pname name];
	W     = [W W];
end

% finished
%-----------------------------------------------------------------------
set(Finter,'Name','')

