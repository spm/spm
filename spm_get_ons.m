function [sf,Cname,Pv,Pname,DSstr] = spm_get_ons(k,T,dt,STOC)
% returns onset times for events
% FORMAT [sf,Cname,Pv,Pname,DSstr] = spm_get_ons(k,T,dt,STOC)
%
% k     - number of scans
% T     - time bins per scan
% dt    - time bin length (secs)
% STOC  - flag to enable stochastic designs [0 or 1]
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
%-----------------------------------------------------------------------
% %W% Karl Friston %E%


% time bins per scan
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fstr   = get(Finter,'name');
sf     = {};
Cname  = {};
Pv     = {};
Pname  = {};
DSstr  = '';

% get stick functions {ons} and names
%=======================================================================

% get trials
%-----------------------------------------------------------------------
v     = spm_input('number of conditions or trials',1,'w1',0);
for i = 1:v

	% get names
	%---------------------------------------------------------------
	str         = sprintf('trial %d',i);
	Cname{i}    = spm_input('condition or trial name','+1','s',str);
end


% event/epoch-related responses
%-----------------------------------------------------------------------
if v    

	% stochastic designs
	%---------------------------------------------------------------
	set(Finter,'Name','trial specification')
	if STOC
		STOC = spm_input('stochastic design',1,'y/n',[1 0]);
	end
	if STOC

		% minimum SOA
		%-------------------------------------------------------
		ne      = spm_input('include a null event','+1','y/n',[1 0]);
		soa     = spm_input('SOA (scans)','+1','r',2)*T;
		on      = fix(1:soa:(k*T));
		ns      = length(on);
		DSstr   = [DSstr sprintf('Stochastic: %.2fsec SOA ',soa*dt)];

		% occurence probabilities - stationary
		%-------------------------------------------------------
		str     = 'occurence probability';
		if spm_input(str,'+1','stationary|modulated',[1 0])
			DSstr = [DSstr '(stationary) '];
			P     = ones((v + ne),ns);
 
		% occurence probabilities - modulated (32 sec period)
		%-------------------------------------------------------
		else
			DSstr = [DSstr '(modulated) '];
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

		% stick function array (and delete null event)
		%-------------------------------------------------------
		for  i = 1:v
			sf{i}   = full(ons(:,i));
		end

	% non-stochastic designs
	%---------------------------------------------------------------
	else

	    % get onsets
	    %-------------------------------------------------------
	    Sstr  = spm_input('SOA',1,'Fixed|Variable');
	    DSstr = [DSstr  Sstr ' SOA '];
	    i     = 0;
	    while i < v

		% get onsets
		%-------------------------------------------------------
		switch Sstr

			case 'Fixed'
			%-----------------------------------------------
			str   = ['SOA (scans) for ' Cname{i + 1}];
			soa   = spm_input(str,'+1','r');
			on    = spm_input('first trial (scans)','+1','r',1);
			on    = on:soa:k;

			case 'Variable'
			%-----------------------------------------------
			str   = ['vector of onsets (scans) for ' Cname{i + 1}];
			on    = spm_input(str,'+1');
		end

		if iscell(on)

			% create stick functions
			%-------------------------------------------------
			for j = 1:length(on)
				i     = i + 1
	    			ons   = sparse(k*T,1);
				ons(round(on{j}*T + 1)) = 1;
				sf{i} = ons(1:(k*T));
			end
		else
			% create stick functions
			%------------------------------------------------
			i     = i + 1;
	    		ons   = sparse(k*T,1);
			ons(round(on*T + 1)) = 1;
			sf{i} = ons(1:(k*T));
		end
	    end
	end


	% get parameters, contruct interactions and append
	%================================================================
	set(Finter,'Name','Parametric specification')

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
	Ptype = spm_input('parametric modulation ','+1','b',Ptype);
	switch Ptype

		case 'none'
		%--------------------------------------------------------
		set(Finter,'Name',Fstr)
		return

		case 'other'
		%--------------------------------------------------------
		Pstr   = spm_input('name of parameter','+1','s');

		case 'time'
		%--------------------------------------------------------
		Pstr   = Ptype;
	end

	% get parameters of expansion
	%----------------------------------------------------------------
	Etype = {'linear',...
		 'exponen',...
		 'polynom'};
	Etype = spm_input('expansion','+1','b',Etype);
	DSstr = [DSstr  '[ x ' Pstr ' (' Etype ')] '];
	switch Etype

		case 'exponen'
		%--------------------------------------------------------
		if strcmp(Ptype,'time')
			h = round(k*T*dt/4);
			h = spm_input('time constant {secs}','+1','r',h);

		else
			h = spm_input('decay constant','+1','r');
		end

		case 'polynom'
		%--------------------------------------------------------
		str       = 'order of polynomial expansion';
		h         = spm_input(str,'+1','r',2);

	end


	% cycle over selected trial types
	%----------------------------------------------------------------
	str   = sprintf('which trial[s] 1 to %d',v);
	for i = spm_input(str,'+1','e',1)

		set(Finter,'Name',Cname{i})
		on    = find(sf{i});
		ns    = length(on);

		% get parameters
		%-------------------------------------------------------
		switch Ptype

			case 'other'
			%-----------------------------------------------
			str   = ['parameters for ' Cname{i}];
			p     = spm_input(str,'+1','r',[],[ns,1]);

			case 'time'
			%-----------------------------------------------
			p     = on*dt;

		end

		% expansion
		%--------------------------------------------------------
		switch Etype


			case 'polynom'
			%------------------------------------------------
			r              = spm_detrend(p(:));
			q              = [];
			for j = 1:h
				q(:,j) = sparse(on,1,r.^j,size(sf{i},1),1);
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

		% appnd as new trials
		%--------------------------------------------------------
		sf{i}    = [sf{i} q];
		Pv{i}    = p;
		Pname{i} = Pstr;

	end
end


% finished
%------------------------------------------------------------------------
set(Finter,'Name',Fstr)
