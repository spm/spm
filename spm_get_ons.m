function [ons,W,name,para,DSstr] = spm_get_ons(k,T,dt)
% returns onset times for events
% FORMAT [ons,W,name,para,DSstr] = spm_get_ons(k,T,dt)
%
% k    - number of scans
% T    - time bins per scan
% dt   - time bin length (secs)
%
% ons     - [m x n]   stick functions {for n trial-types over m = k*T time-bins}
% W       - [1 x n]   vector of window/epoch lengths {time-bins}
% name    - {1 x n}   cell of names for each trial-type
% para    - {1 x n}   cell of parameter vectors
% DSstr   - Design string
%_______________________________________________________________________
% spm_get_ons contructs a matrix of delta functions specifying the onset
% of events or epochs (or both). These are convolved with a basis set at
% a later stage to give regressors that enter into the deisgn matrix.
% The length of the beasis functions is determined by W:
%
% W = 0 - implies an event
% W > 0 - implies an epoch
%-----------------------------------------------------------------------
% %W% Karl Friston %E%


% time bins per scan
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
name   = {};
para   = {};
W      = [];
DSstr  = '';

% get onsets (ons) (and names {names})
%=======================================================================

% get trials
%-----------------------------------------------------------------------
v     = spm_input('event types (without null event)','+1','w1',0);
u     = spm_input('epoch types (without baseline)  ','+1','w1',0);
ons   = sparse(k*T,v + u);
for i = 1:v

	% get names
	%---------------------------------------------------------------
	str         = sprintf('trial %d',i);
	name{i}     = spm_input('event or trial name','+1','s',str);
end
for i = 1:u

	% get names
	%---------------------------------------------------------------
	str         = sprintf('cond %d',i);
	name{i + v} = spm_input('epoch or condition name','+1','s',str);
end

% event-related responses
%-----------------------------------------------------------------------
if v    

	% stochastic designs
	%---------------------------------------------------------------
	set(Finter,'Name','event specification')
	DSstr       = 'event-related ';
	if spm_input('stochastic design',1,'y/n',[1 0])

		% minimum SOA
		%-------------------------------------------------------
		ne      = spm_input('include a null event','+1','y/n',[1 0]);
		soa     = spm_input('SOA (scans)','+1','r',2)*T;
		on      = fix(1:soa:(k*T));
		ns      = length(on);
		DSstr   = [DSstr sprintf('(Stochastic: %.2fsec SOA) ',soa*dt)];

		% occurence probabilities - stationary
		%-------------------------------------------------------
		str     = 'occurence probability';
		if spm_input(str,'+1','stationary|modulated',[1 0])
			DSstr = [DSstr '-stationary '];
			P     = ones((v + ne),ns);
 
		% occurence probabilities - modulated (32 sec period)
		%-------------------------------------------------------
		else
			DSstr = [DSstr '-modulated '];
			P     = ones((v + ne),ns);
			dc    = 32/dt;
			for i = 1:(v + ne);
				p      = sin(2*pi*(on/dc + (i - 1)/(v + 1)));
				P(i,:) = 1 + p;
			end
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

		% delete null event
		%-------------------------------------------------------
		ons   = full(ons(:,1:v));

	% non-stochastic designs
	%---------------------------------------------------------------
	else


	    % get onsets
	    %-------------------------------------------------------
	    Sstr  = spm_input('SOA','+1','fixed|variable');
	    DSstr = [DSstr '(' Sstr ' SOA) '];
	    i     = 0;
	    while i < v

		% set name
		%-------------------------------------------------------
		set(Finter,'Name',name{i + 1})
		switch Sstr

			case 'fixed'
			%-----------------------------------------------
			soa   = spm_input('SOA (scans)','+1','r');
			on    = spm_input('first trial (scans)','+1','r',1);
			on    = on:soa:k;

			case 'variable'
			%-----------------------------------------------
			str   = 'vector [or array] of onsets (scans)'
			on    = spm_input(str,'+1');
		end

		if iscell(on)

			% create stick functions
			%-------------------------------------------------
			for j = 1:length(on)
				i = i + 1
				ons(round(on{j}*T + 1),i) = 1;
			end
		else
			% create stick functions
			%------------------------------------------------
			i = i + 1;
			ons(round(on*T + 1),i) = 1;
		end

	    end
	    v = i;
	end

	% basis function lengths (W == 0 implies event-related)
	%---------------------------------------------------------------
	W     = zeros(1,v);

end

% epoch-related responses
%-----------------------------------------------------------------------
if u

	% set name
	%-------------------------------------------------------
	set(Finter,'Name','epoch specification')
	if length(DSstr)
		DSstr  = [DSstr '& epoch-related '];
	else
		DSstr  = 'epoch-related ';
	end

	% vector of conditions
	%---------------------------------------------------------------
	a     = spm_input('epoch order eg 01020102...',1,'c');
	while (max(a) ~= u) | (min(a) ~= 0)
		str = sprintf('re-enter using 0 to %d',u)
		a   = spm_input(str,2,'c');
	end
	a     = a - min(a) + 1;

	% epoch lengths
	%---------------------------------------------------------------
	w     = zeros(1,u + 1);
	while sum(w(a)) ~= k
		str = sprintf('scans per epoch 0 to %d',u);
		w   = spm_input(str,2,'n');
		while length(w) < (u + 1), w = [w w]; end
		w   = w(1:u + 1);
	end
	w     = w*T;

	% epoch onsets (discounting baseline)
	%---------------------------------------------------------------
	on    = spm_input('start of first epoch {scans}','+1','r',1);
	c     = cumsum(w(a)) - w(a) + on*T;
	for i = 1:u
		ons(round(c(find(a == (i + 1))) + 1),i + v) = 1;
	end

	% basis function lengths
	%---------------------------------------------------------------
	W     = [W w([1:u] + 1)];

end


% get parameters, contruct interactions and append
%=======================================================================
set(Finter,'Name','')
v     = length(W);
if v
    if spm_input('model parametric or time effects',1,'y/n',[1 0])


	% cycle over selected trial types
	%---------------------------------------------------------------
	str   = sprintf('which effect[s] 1 to %d',v);
	for i = spm_input(str,'+1','e',1)

		% basis functions - Type
		%-------------------------------------------------------
		set(Finter,'Name',name{i})
		Ptype = str2mat(...
			'User specified ',...
			'Time (exponential)',...
			'Time (linear)');
		str   = 'Form of parametric modulation';
		Pov   = spm_input(str,3,'m',Ptype,[1:size(Ptype,1)]);
		Pstr  = {'param' 'time' 'time'};
		Pstr  = [name{i} ' x ' Pstr{Pov}];
		on    = find(ons(:,i));
		ns    = length(on);
		p     = [];

		% get parameters
		%-------------------------------------------------------
		if     Pov == 1

			while length(p) ~= ns
				str = 'vector of parameters';
				p   = spm_input(str,'+1','r',[],[1,ns]);
			end

		% exponential adaptation
		%-------------------------------------------------------
		elseif Pov == 2

			tau = round(k*T*dt/4);
			tau = spm_input('time constant {secs}','+1','r',tau);
			p   = exp(-on*dt/tau)';

		% linear adaptation
		%-------------------------------------------------------
		elseif Pov == 3

			p   = on*dt;
		end
		p     = spm_detrend(p(:));


		% append
		%-------------------------------------------------------
		p     = sparse(on,1,p,size(ons,1),1);
		ons   = [ons p];
		name  = [name {Pstr}];
		W     = [W W(i)];

	end
    end
end

% paramteric representation of causes
%-----------------------------------------------------------------------
for i = 1:size(ons,2)
	para{i} = ons(find(ons(:,i)),i);
end

% finished
%-----------------------------------------------------------------------
set(Finter,'Name','')
spm_clf(Finter)
