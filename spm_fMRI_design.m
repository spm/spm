function [X,Sess] = spm_fMRI_design(nscan,RT)
% Assembles a deign matrix for fMRI studies
% FORMAT [X,Sess] = spm_fMRI_design(nscan,RT)
%
% nscan   - n vector {nscan(n) = number of scans in session n}
% RT      - intercans interval {seconds}
%
% X.DesN  - {1} general description, {2} basis function description
% X.xX    - regressors
% X.bX    - session effects
% X.cX    - low frequency confounds
% X.Xname - names of subpartiton columns {1xn}
% X.Bname - names of subpartiton columns {1xn}
% X.Cname - names of subpartiton columns {1xn}
%
% For session s:
% Sess{s}.row    - scan indices 
% Sess{s}.col    - zeroth index of design matrix columns
% Sess{s}.ind{i} - index of design matrix columns for event/epoch i
% Sess{s}.bf{i}  - array of basis functions       for event/epoch i
% Sess{s}.pst{i} - array of peristimulus times    for event/epoch i
%_______________________________________________________________________
% %W% Karl Friston %E%

% construct Design matrix {X} - cycle over sessions
%===========================================================================

% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
nsess  = length(nscan);
dt     = 0.1;						% time bin {secs}
T      = RT/dt;						% bins per scan


% default cutoff period (bin) for High-pass filtering
%---------------------------------------------------------------------------
CUT    = max(nscan)*RT/dt;

% separate speciofications for nonrelicated session
%--------------------------------------------------------------------------
rep   = 0;
if nsess > 1 & ~any(nscan - nscan(1))
	rep = spm_input(['are sessions replicated exactly'],2,'yes|no',[1 0]);
end
xX    = [];
bX    = [];
cX    = [];
Xname = {};
Bname = {};
Cname = {};
Sess  = {};
for s = 1:nsess


	% specify design for this session?
	%-------------------------------------------------------------------
	set(Finter,'Name',sprintf('Session %0.0f',s));
        if   (s == 1) | ~rep
	X       = [];
	C       = [];
	D       = [];
	qx      = 1;
	qc      = 1;
	Xn      = {};
	Cn      = {};
	PST     = {};
	BF      = {};
	IND     = [];

	% Event/epoch related responses			
	%===================================================================
	k       = nscan(s);

	% Study type - epoch or event-related fMRI
	%-------------------------------------------------------------------
	str     = 'epoch or event-related fMRI';
	ER      = spm_input(str,1,'epoch|event',{'epoch' 'event'});


	% event/epoch onsets {ons} and window lengths {W}
	%-------------------------------------------------------------------
	[ons W] = spm_get_ons(ER,k,T);


	% get basis functions for responses
	%-------------------------------------------------------------------
	[BF BFstr]    = spm_get_bf(ER,W,dt);


	% peri-stimulus time {PST} in seconds
	%-------------------------------------------------------------------
	ntrial  = length(W);
	for   i = 1:ntrial

		e     = find(ons(:,i))*dt;
		pst   = [1:k]*RT - e(1);			
		for j = 1:length(e)
			u      = [1:k]*RT - e(j);
			v      = find(u >= -1);
			pst(v) = u(v);
		end
		PST{i} = pst;
	end

	% cutoff period {bins}
	%-------------------------------------------------------------------
	for   i = 1:ntrial
		CUT   = min([CUT max(diff(find(ons(:,i))))]);
	end


	% convolve with basis functions to create design matrix {X}
	%-------------------------------------------------------------------
	for i = 1:ntrial
		IND{i} = [];
		for  j = 1:size(BF{i},2)
			X      = [X conv(ons(:,i),BF{i}(:,j))];
			Xn{qx} = [ER{1} sprintf(': %i bf: %i',i,j)];
			IND{i} = [IND{i} qx];
			qx     = qx + 1;
		end
	end

	% resample X at times of acquisition
	%-------------------------------------------------------------------
	X     = X([1:k]*T,:);


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


	% time x trial interactions
	%===================================================================
	str  = 'time x response interactions';
	if spm_input(str,1,'yes|no',[1 0]);

		% basis functions - Type
		%-----------------------------------------------------------
		Ttype = str2mat(...
			'Exponential',...
			'Linear');
		str   = 'Form of adaptation';
		Tov   = spm_input(str,2,'m',Ttype,[1:size(Ttype,1)]);


		% get paramters
		%-----------------------------------------------------------
		if Tov == 1
			tau = round(k*RT/4);
			tau = spm_input('time constant {secs}',3,'e',tau);
		end

		% compute interaction effects
		%-----------------------------------------------------------
		str   = 'model for which ';
		str   = [str ER{1} '[s] 1 to ' sprintf('%d',ntrial)];
		for i = spm_input(str,3,'e',1)

			% basis functions - create
			%---------------------------------------------------
			t     = [];
			if     Tov == 1
				t  = exp(-[0:(k - 1)]*RT/tau)';

			elseif Tov == 3
				t  = [0:(k - 1)]'/(k - 1);
			end
			
			% append
			%---------------------------------------------------
			t     = spm_detrend(t);
			Q     = X(:,IND{i});
			Q     = spm_detrend(Q);
			for j = 1:size(Q,2);
				X      = [X Q(:,j).*t];
				Xn{qx} = [ER{1} sprintf(' %d x Time',i)];
				qx     = qx + 1;
			end
		end
	end

	% Confounds: Session effects 
	%===================================================================
	B      = ones(k,1);
	Bn{1}  = sprintf('constant');

	% and high pass filter {discrete cosine set}
	%-------------------------------------------------------------------
	if spm_input('high pass filter',1,'yes|no',[1 0])
		str   = 'cut-off period {secs}';
		CUT   = spm_input(str,2,'e',2*CUT*dt);
		u     = find([1:k/2] <= 2*(k*RT)/CUT);
		for i = 1:length(u)
			C      = [C cos(pi*[0:(k - 1)]'*u(i)/(k - 1))];
			Cn{qc} = sprintf('Low Hz: %i',i);
			qc     = qc + 1;

		end
	end

	% end specification
	%-------------------------------------------------------------------
	end

	% Session structure
	%-------------------------------------------------------------------
	Sess{s}.row = [1:k] + size(xX,1);
	Sess{s}.col = size(xX,2);
	Sess{s}.ind = IND;
	Sess{s}.bf  = BF;
	Sess{s}.pst = PST;

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
	q     = length(Cname);
	for i = 1:length(Cn) 
		Cname{q + i}  = [sprintf('Sn: %i ',s) Cn{i}];
	end


	% append into xX and cX
	%===================================================================
	[x y]   = size(xX);
	[i j]   = size(X);
	xX(([1:i] + x),([1:j] + y)) = X;
	[x y]   = size(bX);
	[i j]   = size(B);
	bX(([1:i] + x),([1:j] + y)) = B;
	[x y]   = size(cX);
	[i j]   = size(C);
	cX(([1:i] + x),([1:j] + y)) = C;

end

% finished
%---------------------------------------------------------------------------
X.DesN  = {[ER{1} '-related fMRI'] {BFstr}};
X.xX    = spm_detrend(xX);
X.bX    = bX;
X.cX    = spm_detrend(cX);
X.Xname = Xname;
X.Bname = Bname;
X.Cname = Cname;
