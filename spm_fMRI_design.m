function [X,Sess] = spm_fMRI_design(nscan,RT,BM)
% Assembles a design matrix for fMRI studies
% FORMAT [X,Sess] = spm_fMRI_design(nscan,RT,BM)
%
% nscan   - n vector {nscan(n) = number of scans in session n}
% RT      - intercans interval {seconds}
% BM      - Burst-mode flag
%
% X.dt    - time bin {secs}
% X.RT    - Repetition time {secs}
% X.DesN  - X.DesN{1} study type, X.DesN{2} basis function description
% X.xX    - regressors
% X.bX    - session effects
% X.Xname - names of subpartiton columns {1xn}
% X.Bname - names of subpartiton columns {1xn}
%
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
% %W% Karl Friston %E%

% construct Design matrix {X} - cycle over sessions
%===========================================================================

% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');

% get nscan and RT if no arguments
%---------------------------------------------------------------------------
if nargin < 3
	BM     = spm_input('Burst mode?',1,'b','no|yes',[0 1]);
end
if nargin < 2
	RT     = spm_input('Interscan interval {secs}',2);
end
if nargin < 1
	nscan  = spm_input(['scans per session e.g. 64 64 64'],1);
end
nsess  = length(nscan);
T      = 16;						% bins per scan
dt     = RT/T;						% time bin {secs}


% separate specifications for non-relicated sessions
%--------------------------------------------------------------------------
rep   = BM;
if nsess > 1 & ~any(nscan - nscan(1)) & ~rep
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
	IND     = [];

	% Event/epoch related responses			
	%===================================================================
	k       = nscan(s);

	% Study type - epoch or event-related fMRI
	%-------------------------------------------------------------------
	SType   = {'event-related (stochastic: stationary)',...
		   'event-related (stochastic: modulated)',...
		   'event-related (deterministic: fixed SOA)',...
		   'event-related (deterministic: variable SOA)',...
		   'epoch-related responses'};
	str     = 'Model';
	ST      = spm_input(str,1,'m',SType);
	STstr   = deblank(SType{ST});


	% event/epoch onsets {ons} and window lengths {W}
	%-------------------------------------------------------------------
	[ons,W,name,para] = spm_get_ons(ST,k,T,dt);


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
			X      = [X x([1:k]*T,:)];
			Xn{qx} = [name{i} sprintf(': bf: %i',j)];
			IND{i} = [IND{i} qx];
			qx     = qx + 1;
		end
	end


	% get user specified regressors
	%===================================================================
	c     = spm_input('user specified regressors',1,'e',0);
	while size(D,2) < c
		str = sprintf('[%d]-variate %i',k,size(D,2) + 1);
		d   = spm_input(str,2);
		if size(d,2) == k, d = d';    end
		if size(d,1) == k, D = [D d]; end	
	end

	% append regressors and names
	%-------------------------------------------------------------------
	for i = 1:size(D,2)
		X      = [X D(:,i)];
		Xn{qx} = sprintf('regressor: %i',i);
		qx     = qx + 1;
	end


	% Confounds: Session effects 
	%===================================================================
	if ~BM
		B      = ones(k,1);
		Bn{1}  = sprintf('constant');
	end

	% end specification
	%-------------------------------------------------------------------
	end

	% Session structure
	%-------------------------------------------------------------------
	Sess{s}.row  = size(xX,1) + [1:k];
	Sess{s}.col  = size(xX,2) + [1:size(X,2)];
	Sess{s}.name = name;
	Sess{s}.ind  = IND;
	Sess{s}.bf   = BF;
	Sess{s}.pst  = PST;
	Sess{s}.ons  = ONS;
	Sess{s}.para = para;

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

% Burst mode - model main effect of scan
%---------------------------------------------------------------------------
if BM
	bX    = kron(ones(nsess,1),eye(k));
	for i = 1:k
		Bname{i} = sprintf('scan: %i ',i);
	end
end

% finished
%---------------------------------------------------------------------------
X.dt    = dt;
X.RT    = RT;
X.DesN  = {STstr BFstr{:}};
X.xX    = xX;
X.bX    = bX;
X.Xname = Xname;
X.Bname = Bname;

save fMRIDesMtx X Sess
