function [ons,W,name,para] = spm_get_ons(ST,k,T,dt)
% returns onset times for events
% FORMAT [ons,W,name,para] = spm_get_ons(ST,k,T,dt)
%
% ST   - Study type
%    1    - event-related (stochastic: stationary)
%    2    - event-related (stochastic: nodulated)
%    3    - event-related (deterministic: fixed SOA)
%    4    - event-related (deterministic: variable SOA)
%    5    - epoch-related responses
% 
% k    - scans per session
% T    - time bins per scan
% dt   - time bin length (secs)
%
% ons  - [m x n]   stick functions {for n trial-types over m = k*T time-bins}
% W    - [1 x n]   vector of window/epoch lengths {time-bins}
% name - {1 x n}   cell of names for each trial-type
% para - {1 x n}   cell of parameter vectors {when specified}
%_______________________________________________________________________
% %W% Karl Friston %E%

% time bins per scan
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
name   = {};
para   = {};

% get onsets (ons) (and names {names})
%-----------------------------------------------------------------------
if  ST == 1 | ST == 2 | ST == 3 | ST == 4

	% get trials
	%---------------------------------------------------------------
	v     = spm_input('number of event or trial types',1,'e',1);
	ons   = sparse(k*T,v);
	for i = 1:v

		% get name
		%-------------------------------------------------------
		str     = sprintf('trial %d',i);
		name{i} = spm_input('event or trial name',2,'s',str);
	end


	% stochastic designs
	%---------------------------------------------------------------
	if ST == 1 | ST == 2

		% minimum SOA
		%-------------------------------------------------------
		soa     = spm_input('SOA (scans)',1,'e',2)*T;
		on      = fix(1:soa:(k*T));
		ns      = length(on);
		
		% occurence probabilities - stationary
		%-------------------------------------------------------
		if  ST == 1

			P     = ones((v + 1),ns);
 
		% occurence probabilities - modulated (32 sec period)
		%-------------------------------------------------------
		elseif ST == 2

			P     = ones((v + 1),ns);
			dc    = 32/dt;
			for i = 1:(v + 1);
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
		for i = 1:(v + 1);
			j       = find(Q >= P(i,:) & Q < P(i + 1,:));
			q(j)    = i;
		end

		% create stick functions
		%-------------------------------------------------------
		ons   = sparse(on,q,1,k*T,v + 1);

		% delete null event
		%-------------------------------------------------------
		ons   = full(ons(:,1:v));
	end

	% non-stochastic designs
	%---------------------------------------------------------------
	if ST == 3 | ST == 4
	for i = 1:v
	
		% set name
		%-------------------------------------------------------
		set(Finter,'Name',name{i})

		% get onsets
		%-------------------------------------------------------
		if     ST == 3

			soa = spm_input('SOA (scans)',3);
			on  = spm_input('first trial (scans)',4,'e',1);
			on  = on:soa:k;

		elseif ST == 4

			on  = spm_input('trial onset times (scans)',3);
		end


		% create stick functions
		%-------------------------------------------------------
		ons(round(on*T + 1),i) = 1;

	end
	end

	% basis function lengths (not used subsequently)
	%---------------------------------------------------------------
	W     = zeros(1,v);


elseif ST == 5


	% vector of conditions
	%---------------------------------------------------------------
	str   = 'epoch order eg 1 2 1...';
	a     = spm_input(str,1);
	v     = max(a);
	ons   = zeros(k*T,v);


	% epoch names
	%---------------------------------------------------------------
	for i = 1:v
		str     = sprintf('name of epoch %d',i);
		name{i} = spm_input(str,2,'s',sprintf('cond %d',i));
	end

	% epoch lengths
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
	on    = spm_input('start of first epoch {scans}',3,'e',0);
	c     = cumsum(W(a)) - W(a) + on*T;
	for i = 1:v
		ons(round(c(find(a == i)) + 1),i) = 1;
	end

end


% get parameters if 'parametric', contruct interactions and pre-pend
%=======================================================================
if spm_input('model parametric or time effects',1,'b','no|yes',[0 1])


	% cycle over selected trial types
	%---------------------------------------------------------------
	str   = sprintf('which trial[s] 1 to %d',v);
	for i = spm_input(str,2,'e',1)

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
				str = sprintf('[%d]-vector',ns);
				p   = spm_input(str,4);
			end

		% exponential adaptation
		%-------------------------------------------------------
		elseif Pov == 2

			tau = round(k*T*dt/4);
			tau = spm_input('time constant {secs}',4,'e',tau);
			p   = exp(-on*dt/tau)';

		% linear adaptation
		%-------------------------------------------------------
		elseif Pov == 3

			p   = on/max(on);
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



% paramteric representation of causes
%-----------------------------------------------------------------------
for i = 1:size(ons,2)
	para{i} = ons(find(ons(:,i)),i);
end

% finished
%-----------------------------------------------------------------------
set(Finter,'Name','')

% clf
%-----------------------------------------------------------------------
spm_clf(Finter)
