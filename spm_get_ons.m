function [ons,W] = spm_get_ons(ER,k,T)
% returns onset times for events
% FORMAT [ons W] = spm_get_ons(ER,k,T)
% ER  - 'event'|'epoch'
% k   - scans per session
% T   - time bines per scan
%
% ons - {n x m} stick functions {onsets for n epochs/events over m time bins}
% W   - {1 x n} vector of window/epoch lengths {time-bins}
%_______________________________________________________________________
% %W% Karl Friston %E%

% time bins per scan
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if     strcmp(ER,'event')
 
	% cycle over events
	%---------------------------------------------------------------
	v     = spm_input('number of event types',1,'e',1);
	str   = 'inter-stimulus intervals';
	ISI   = spm_input(str,'!+1','fixed|variable',[1 0]);
	ons   = zeros(k*T,v);
	for i = 1:v
	
		% get onsets
		%-------------------------------------------------------
		set(Finter,'Name',sprintf('event type %i',i))
		if ISI
			isi = spm_input('interstimulus interval (scans)',2);
			on  = spm_input('first stimulus (scans)','!+0','e',1);
			on  = on:isi:k;
		else
			on  = spm_input('stimulus onset times (scans)','!+0');
		end

		% cumulate
		%-------------------------------------------------------
		j        = round(on*T + 1);
		ons(j,i) = ones(length(j),1);

	end
	W     = ones(1,v);


elseif strcmp(ER,'epoch')


	% vector of conditions
	%---------------------------------------------------------------
	str   = 'epoch order eg 1 2 1...';
	a     = spm_input(str,1);
	v     = max(a);
	
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

% finished
%-----------------------------------------------------------------------
set(Finter,'Name','')

