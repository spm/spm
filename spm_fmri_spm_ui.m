function spm_fmri_spm_ui
% Setting up the general linear model for fMRI time-series
% FORMAT spm_fmri_spm_ui
%____________________________________________________________________________
%
% spm_fmri_spm_ui configures the design matrix, data specification and
% thresholds that specify the ensuing statistical analysis. These
% arguments are passed to spm_spm that then performs the actual analysis.
%
% This routine has two uses:
%
% i)  To implement a unvariate analysis testing for effects of interest
%     at each and every voxel.  Effects of interest are specified as
%     covariates of interest and confounds are specified as covariates
%     of no interest.
%
% ii) To remove any confounds as a preprocessing step for further analysis
%      e.g. SVD.  In this instance only confounds are specified and the
%      adjusted data will be found in XA.mat
%
% The design matrix defines the experimental design and the nature of
% hypothesis testing to be implemented.  The design matrix has one row
% for each scan and one column for each effect or parameter.  These
% parameters (e.g. reference waveform, subject or block effect, 
% regression slope of voxel on a confound etc) are estimated in a
% least squares sense using the general linear model.  Specific profiles
% within these parameters are tested using a linear compound or CONTRAST
% with the t statistic.  The resulting map of t values constitutes the
% SPM{t}.  The SPM{t} is then characterized in terms of focal or regional
% differences by assuming that (under the null hypothesis) the SPM{t}
% behaves as a smooth stationary Gaussian field.
%
%     From the user's perspective it is important to specify the design
% matrix and contrasts correctly.  The design matrix is built when you
% specify the number of sessions/subjects and conditions.  The covariates 
% (that constitute the columns of the design matrix) can be thought of as
% reference vectors and can be specified as such.  Alternatively one
% can specify reference vectors in terms of response functions or
% waveforms: Waveforms are specified for each EPOCH of scans that
% constitute a particular condition.  Note that if there are two
% conditions, then two covariates are specified each expressing the
% same waveform[s] every time the condition occurs.  A waveform is
% simply a transient response to the onset of an epoch, that lasts for
% the duration of that epoch.  The form of this response waveform can
% be fixed (e.g. a box-car or half sine wave) or allowed to vary
% between conditions by using two (exponentially modulated sine)
% waveforms.  In the latter case there are two covariates for each
% condition (and the CONTRAST must be specified with this in mind).
% These two response functions correspond to the early and late
% components of a modeled response and differential adaptation of the
% hemodynamics during condition-specific epochs can be tested with the
% appropriate CONTRAST (see the final reference below) 
% 
% If you use a box-car (i.e. square wave) function you can optionally
% include its temporal derivative as an additional covariate of
% interest.  The box-car function is delayed by 6 seconds, however this
% may be two much or too little for some brain regions.  The timing
% covariate (temporal derivative) models a small shift in time that best
% fits the data.  Inferences about whether the observed delay is
% significantly different than 6 seconds can be tested using contrasts in
% the usual way (a positive parameter estimate means a shorter delay).
% 
% Covariates of no interest (called confounds) can also be specfied.  You
% will be prompted for some specific confounds such as low frequency
% artifacts (and whole brain activity).
%
% Epochs can vary in length (and order) within and between subjects or runs.
% If multiple subjects or sessions are specified, then subject or run-specific
% waveforms are used.  This means that main effects of conditions and
% interactions between conditions and subjects (or runs) can be evaluated
% with the appropriate contrast.  If you want to treat all your sessions (or
% subjects) as one then specify just one session/subject.
%
% The way that epochs or successive conditions are specified is now more
% intuitive and flexible.  If there are 3 conditions just type in the
% conditions in the order they were presented i.e. 1 2 3 3 2 1 ....
% Later you will be asked to specify the number of scans for each epoch,
% again as a vector (list of numbers).  If the epochs were all the same
% length, then just type in that length once.
%
% The CONTRAST is simply a list or vector of coefficients that are used
% to test for a pattern of effects.  The number of coefficients (length
% of the CONTRAST) should be the same as the number of covariates of interest
% By specifying different contrasts one can effect a wide variety of analyses.
%
% ER.mat is saved if event-related designs ares specified
%-----------------------------------------------------------------------
% variables pertaining to event-related fMRI in ER.mat
%
%	ERI 	- indices of DER pertaing to each event type
%	DER 	- Design matrix containing the basis functions for each event
%	PST	- Peri-stimulus time for each scan
%
% Refs:
%
% Friston KJ, Holmes A, Poline J-B, Grasby PJ, Williams SCR, Frackowiak
% RSJ & Turner R (1995) Analysis of fMRI time-series revisited. NeuroImage
% 2:45-53
%
% Worsley KJ and Friston KJ (1995) Analysis of fMRI time-series revisited -
% again. NeuroImage 2:178-181
%
% Friston KJ, Frith CD, Frackowiak RSJ, & Turner R (1995) Characterising
% dynamic brain responses with fMRI: A multivariate approach NeuroImage -
% 2:166-172
%
% Frith CD, Turner R & Frackowiak RSJ (1995) Characterising evoked 
% hemodynamics with fMRI Friston KJ, NeuroImage 2:157-165
%
% Josephs O, Turner R and Friston KJ (1997) Event-related fMRI, Hum. Brain
% Map. 0:00-00
%
%___________________________________________________________________________
% %W% Karl Friston, Jean-Baptiste Poline %E%


%-Delete files from previous analyses, if they exist
%---------------------------------------------------------------------------
spm_unlink ER.mat


% get filenames and other user specified parameters
%===========================================================================

% Initialize variables
%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');

Q      = [];				% matrix of filename strings
H      = [];				% Factors of interest
C      = [];				% covariates of interest
B      = [];				% Factors of no interest
G      = [];				% covariates of no interest
F      = [];				% Low frequency confounds
Hnames = ' ';				% design matrix effect labels
Bnames = ' ';				% design matrix effect labels
Cnames = ' ';				% design matrix effect labels
Gnames = ' ';				% design matrix effect labels
GM     = 100;				% grand mean for fMRI
TD     = 0;				% switch for temporal differences
REP    = 0;				% switch for session replications
ons    = [];				% epoch or event onset times
PST    = [];				% peri-stimulus times

CONTRAST = [];				% row matrix of contrasts



% get filenames
%---------------------------------------------------------------------------
set(Finter,'Name','fMRI analysis');
nsess  = spm_input(['# of sessions or subjects'],1,'e',1);
nscan  = zeros(1,nsess);
if nsess > 1
	REP = spm_input(['are these replications'],'!+0','yes|no',[1 0]);
end
for i  = 1:nsess
	str      = sprintf('select scans for session %0.0f',i);
	P        = spm_get(Inf,'.img',str);
 	Q        = str2mat(Q,P);
	nscan(i) = size(P,1);
end
Q(1,:) = []; P = Q;
q      = size(P,1);


% get ORIGIN
%---------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));

% Threshold for F statistic
%---------------------------------------------------------------------------
global UFp
UFp   = spm_input('threshold for F, p = ','!+1','e',0.001);

% get Repeat time
%---------------------------------------------------------------------------
RT    = spm_input('Interscan interval {secs}','!+0');


% default cutoff period (scans) for Hi-pass filtering
%---------------------------------------------------------------------------
CUT   = max(nscan);

% Study type - epoch or event-related fMRI
%---------------------------------------------------------------------------
ER    = spm_input('epoch or event-related fMRI','!+1','epoch|event',[0 1]);

% estimated hemodynamic response function
%---------------------------------------------------------------------------
HRF   = 1;					% convolve with hrf {default}
hrf   = spm_hrf(RT);

% construct responses C [and G]
%---------------------------------------------------------------------------
if ER

	% model event-related responses
	%-------------------------------------------------------------------
	Ctype = str2mat(...
		'basis functions (Fourier set)',...
		'basis functions (Windowed Fourier set)',...
		'basis functions (Gamma functions with derivatives)',...
		'basis functions (Gamma functions)',...
		'hrf (with derivative)',...
		'hrf (alone)');
	str   = 'Select basis set';
	Cov   = spm_input(str,'!+1','m',Ctype,[1:size(Ctype,1)]);


	% create small design matrix of basis functions
	%-------------------------------------------------------------------
	dt     = 0.1;					% time step {sec}
	if any(Cov == [1 2])

		% Windowed (Hanning) Fourier set
		%-----------------------------------------------------------
		str   = 'window length {secs}';
		pst   = spm_input(str,'!+0','e',32);
		pst   = [0:dt:pst]';
		str   = sprintf('order {max = %i}',fix(32/RT/2));
		h     = spm_input(str,'!+0','e',4);

		% hanning window
		%-----------------------------------------------------------
		L     = max(pst);
		if Cov == 1
			g = ones(size(pst));
		else
			g = (1 - cos(2*pi*(pst)/L))/2;
		end

		% zeroth and higher terms
		%-----------------------------------------------------------
		DER   = g;
		for i = 1:h
			DER = [DER g.*sin(i*2*pi*pst/L)];
			DER = [DER g.*cos(i*2*pi*pst/L)];	
		end

	elseif Cov == 3

		% Gamma functions and derivatives
		%-----------------------------------------------------------
		pst   = 0:dt:32;
		dx    = 0.01;
		D     = spm_Volt_W(pst);
		DER   = [D (spm_Volt_W(pst - dx) - D)/dx];

	elseif Cov == 4

		% Gamma functions alone
		%-----------------------------------------------------------
		pst   = 0:dt:32;
		DER   = spm_Volt_W(pst);

	elseif Cov == 5

		% hrf and derivatives
		%-----------------------------------------------------------
		D     = spm_hrf(dt);
		pst   = [0:(length(D) - 1)]*dt;
		DER   = [D' gradient(D)'];

	elseif Cov == 6

		% hrf and derivatives
		%-----------------------------------------------------------
		D     = spm_hrf(dt);
		pst   = [0:(length(D) - 1)]*dt;
		DER   = D';

	end

	% scans (model the same event for all sessions)
	%-------------------------------------------------------------------
	h     = size(DER,2);

	% for each session
	%-------------------------------------------------------------------
	for s = 1:nsess


	    % reset name and session partitions
	    %---------------------------------------------------------------
	    set(Finter,'Name',sprintf('Session or subject %0.0f',s));
	    E     = [];
	    CD    = [];
	    GD    = [];
	    PSTD  = [];

	    % scans for this session
	    %---------------------------------------------------------------
	    k     = nscan(s);

	    % cycle over events
	    %---------------------------------------------------------------
	    for v = 1:spm_input('number of event types',4,'e',1)
	
		str   = 'inter-stimulus intervals';
		if spm_input(str,'!+1','fixed|variable',[1 0]);
			isi = spm_input('interstimulus interval (scans)','!+0');
			ons = spm_input('first stimulus (scans)','!+0','e',1);
			ons = ons:isi:k;
		else
			ons = spm_input('stimulus onset times (scans)','!+0');
		end

		% sort in ascending order
		%-----------------------------------------------------------
		ons   = sort(ons);

		% cutoff period
		%-----------------------------------------------------------
		CUT   = min([CUT max(diff(ons))]);

		% peri-stimulus time {PST} in seconds
		%-----------------------------------------------------------
		pst   = ([1:k] - ons(1))*RT;			
		for i = 1:length(ons)
			u  = ([1:k] - ons(i));
			j  = find(u >= -1);
			pst(j) = u(j)*RT;
		end

		% create stimulus functions
		%-----------------------------------------------------------
		D     = [];
		SF    = full(sparse(ons*RT/dt,1,1,k*RT/dt,1));
		j     = round([1:k]*RT/dt);
		for i = 1:h
			d = conv(SF,DER(:,i));
			D = [D d(j)];
		end

		% append to E and PSTD
		%-----------------------------------------------------------
		PSTD  = [PSTD pst(:)];
		E     = [E D];

	    end


	    % modeling evoked responses
	    %---------------------------------------------------------------
	    nevent = v;
	    if nevent > 1
		Cinf  = str2mat(...
			'all events',...
			'some events (others = confounds)',...
			'pairwise differences');
		str   = 'Make inferences about';
		Cf    = spm_input(str,'!+1','m',Cinf,[1:size(Cinf,1)]);
	    else
		Cf    = 1;
	    end


	    % all events
	    %---------------------------------------------------------------
	    if Cf == 1

		% append to CD
		%-----------------------------------------------------------
		CD    = [CD E];

		% append labels
		%-----------------------------------------------------------
		for i = 1:nevent
			for j = 1:h
				str    = sprintf('Sess %0.0f Event %0.0f',s,i);
				Cnames = str2mat(Cnames,str);
			end
		end
	    end

	    % some events
	    %---------------------------------------------------------------
	    if Cf == 2

		d     = 1:nevent;
		u     = spm_input('events of interest {e.g 1 2}','!+0');
		d(u)  = [];

		% select relevent psts
		%-----------------------------------------------------------
		PSTD  = PSTD(:,u);
		for i = 1:length(u)

			% append to CD
			%---------------------------------------------------
			CD     = [CD E(:,([1:h] + (u(i) - 1)*h))];

			% append labels
			%---------------------------------------------------
			for j = 1:h
			    str    = sprintf('Sess %0.0f Event %0.0f',s,u(i));
			    Cnames = str2mat(Cnames,str);
			end
		end
		for i = 1:length(d)

			% append to GD
			%---------------------------------------------------
			GD     = [GD E(:,([1:h] + (d(i) - 1)*h))];

			% append labels
			%---------------------------------------------------
			for j = 1:h
			    str    = sprintf('Sess %0.0f Event %0.0f',s,d(i));
			    Gnames = str2mat(Gnames,str);
			end
		end
	    end


	    % differences
	    %---------------------------------------------------------------
	    if Cf == 3

		d     = 1:nevent;
		u     = spm_input('which pair of events {e.g 1 2}','!+1');
		d(u)  = [];

		% no relevent psts
		%-----------------------------------------------------------
		PST   = [];
		a     = E(:,([1:h] + (u(1) - 1)*h));
		b     = E(:,([1:h] + (u(2) - 1)*h));


		% append to C
		%-----------------------------------------------------------
		CD     = [CD (a - b)];

		% append to G
		%-----------------------------------------------------------
		GD     = [GD (a + b)];

		% append labels
		%-----------------------------------------------------------
		for j = 1:h
			str    = sprintf('Events %0.0f - %0.0f (%0.0f)',...
				 u(1),u(2),j);
			Cnames = str2mat(Cnames,str);
			str    = sprintf('Events %0.0f + %0.0f (%0.0f)',...
				 u(1),u(2),j);
			Gnames = str2mat(Gnames,str);
		end

		for i = 1:length(d)

			% append to G
			%---------------------------------------------------
			GD     = [GD E(:,([1:h] + (d(i) - 1)*h))];

			% append labels
			%---------------------------------------------------
			for j = 1:h
				str    = sprintf('Event %0.0f (%0.0f)',d(i),j);
				Gnames = str2mat(Gnames,str);
			end
		end
	    end

	    % append to C
	    %---------------------------------------------------------------
	    [x y]  = size(C);
	    [i j]  = size(CD);
	    d      = [1:i] + x;
	    C(d,([1:j] + y)) = CD;

	    % append to G
	    %---------------------------------------------------------------
	    [x y] = size(G);
	    [i j] = size(GD);
	    G(d,([1:j] + y)) = GD;

	    % append to PST
	    %---------------------------------------------------------------
	    [x y] = size(PST);
	    [i j] = size(PSTD);
	    PST(d,([1:j] + y)) = PSTD + 1e-8;


	end % (loop over sessions)

	% save design matrix, indices and peri-stimulus times for plotting
	%-------------------------------------------------------------------
	DER   = [zeros(4/dt,size(DER,2)); DER];
	ERI   = reshape([1:size(C,2)],h,size(C,2)/h);
	if length(ERI)
		save ER ERI DER PST
	end


end % (event-related)

% Epoch-related fMRI
%---------------------------------------------------------------------------
if ~ER

	% covariates of interest - Type
	%-------------------------------------------------------------------
	Ctype = str2mat(...
		'basis functions  (Discrete Cosine Set)',...
		'basis functions  (Two Exponential sines)',...
		'fixed response   (Half-sine)',...
		'fixed response   (Box-car)',...
		'Fourier analysis (Harmonics of stimulus function)',...
		'User specified');
	str   = 'Select type of response';
	Cov   = spm_input(str,'!+1','m',Ctype,[1:size(Ctype,1)]);

	% order of basis functions - h
	%-------------------------------------------------------------------
	if Cov == 1
		str = 'number of basis functions';
		h   = spm_input(str,'!+1','e',2);

	elseif Cov == 5
		str = 'number of harmonics';
		h   = spm_input(str,'!+1','e',3);

	elseif Cov == 2

		h   = 2;
	else
		h   = 1;
	end

	% if fixed response ask for temporal differences
	%-------------------------------------------------------------------
	if any(Cov == [3 4 6])
		str = 'add temporal derivative';
		TD  = spm_input(str,'!+1','b','no|yes',[0 1]);
	end


	% for each session
	%-------------------------------------------------------------------
	for v = 1:nsess


	% reset name
	%-------------------------------------------------------------------
	set(Finter,'Name',sprintf('Session or subject %0.0f',v));

	% scans for this session
	%-------------------------------------------------------------------
	k     = nscan(v);

	% user specified
	%-------------------------------------------------------------------
	if Cov == 5

		% if replications assume previous parameters
		%-----------------------------------------------------------
		if (v == 1) | ~REP

			% get covariates of interest
			%---------------------------------------------------
			c     = spm_input('Duty cycle (scans)','!+1');
			D     = [];
			for i = 1:h
				d    = sin(i*2*pi*[1:k]/c);
				D    = [D d'];
				d    = cos(i*2*pi*[1:k]/c);
				D    = [D d'];
			end

			% do not convolve with HRF
			%---------------------------------------------------
			HRF   = 0;
		end

		% add labels
		%-----------------------------------------------------------
	     	for i = 1:size(D,2)
			str    = sprintf('Sess %0.0f harmonic',v);
			Cnames = str2mat(Cnames,str);
		end

	elseif Cov == 6

		% if replications assume previous parameters
		%-----------------------------------------------------------
		if (v == 1) | ~REP

			% get covariates of interest
			%---------------------------------------------------
			c     = spm_input('number of covariates','!+1');
			D     = [];
			while size(D,2) < c
				u   = size(D,2) + 1;
				str = sprintf('[%d]-variate %d, sess %d',k,u,v);
				d   = spm_input(str,'!+1');
				if size(d,2) == k
					d    = d';    
				end
				if size(d,1) == k
					D    = [D d];
				end	
			end

			% convolve with HRF?
			%---------------------------------------------------
			HRF   = spm_input('convolve with hrf',...
					  '!+0','b','no|yes',[0 1]);
		end

		% add labels
		%-----------------------------------------------------------
	     	for i = 1:size(D,2)
			str    = sprintf('Sess %0.0f Cond %0.0f',v,i);
			Cnames = str2mat(Cnames,str);
		end



	else % build response partition
	%-------------------------------------------------------------------

		% if replications assume previous parameters
		%-----------------------------------------------------------
		if (v == 1) | ~REP

			% vector of conditions
			%---------------------------------------------------
			str   = 'epoch order eg 1 2 1..{0 = null}';
			a     = spm_input(str,4);
			a     = a(:); a = a';
			ncond = max(a);
	
			% vector of epoch lengths
			%---------------------------------------------------
			e     = [];
			while any(~e) | sum(e) ~= k
				str = 'scans per epoch eg 8 or 8 6... ';
				e   = spm_input(str,'!+1');
				while length(e) < length(a)
					e = [e e(:)'];
				end
				e   = e(1:length(a));
			end
			e     = e(:);
			e     = e';

			% epoch onsets (starting at 0)
			%---------------------------------------------------
			ons   = cumsum(e) - e;

			% cutoff period
			%---------------------------------------------------
			for i = 0:ncond
				CUT   = min([CUT max(diff(ons(a == i)))]);
			end

		end

		% Assemble basis functions for each epoch
		%-----------------------------------------------------------
		D     = zeros(k,ncond*h);
		for i = 1:length(e)

		    % Skip if null or rest {a(i) = 0}
		    %-------------------------------------------------------
		    if a(i)

		    	for j = 1:h

			% Discrete cos set
			%---------------------------------------------------
			if Cov == 1			
				W = cos((j - 1)*pi*[1:e(i)]/e(i));
				D([1:length(W)]+ons(i),(a(i) - 1)*h + j) = W(:);	

			% Exponential sine functions
			%---------------------------------------------------
			elseif Cov == 2			
				u = [1:e(i)];
				W = sin(pi*u/e(i)).*exp((1.5 - j)*u/2);
				W = W(:)/max(W);
				D([1:length(W)] + ons(i),(a(i) - 1)*h + j) = W;	
	

			% sine wave
			%---------------------------------------------------
			elseif Cov == 3
				W = sin(pi*[0:(e(i) + 1)]'/(e(i) + 1));
				D(([1:size(W,1)] + ons(i)),a(i)) = W/max(W);

			% box car
			%---------------------------------------------------
			elseif Cov == 4
				D([1:e(i)] + ons(i),a(i)) = ones(e(i),1);

			end

			end
		    end
		end

		% trim
		%-----------------------------------------------------------
		D     = D(1:k,:);

		% mean centre box-cars if there are more than one
		%-----------------------------------------------------------
		if (Cov == 4) & min(a)
			D     = spm_detrend(D);
		end

		

		% append labels
		%-----------------------------------------------------------
		for i = 1:ncond
		    for j = 1:h
			str    = sprintf('Sess %0.0f Cond %0.0f (%0.0f)',v,i,j);
			Cnames = str2mat(Cnames,str);
		    end
		end

		% convolve with hemodynamic response function - hrf
		%-----------------------------------------------------------
		if HRF
			d = length(hrf);
			D = [ones(d,1)*D(1,:); D];
			D = spm_sptop(hrf,k + d)*D;
			D = D([1:k] + d,:);
		end

		% add temporal differences if specified
		%-----------------------------------------------------------
		if TD
			% append labels
			%---------------------------------------------------
			for i = 1:size(D,2)
				str    = sprintf('timing (%0.0f)',i);
				Cnames = str2mat(Cnames,str);
			end

			% append to D
			%---------------------------------------------------
			D     = [D [diff(D); zeros(1,size(D,2))]];
		end


	end % (else)

	% append to C
	%-------------------------------------------------------------------
	[x y]  = size(C);
	[i j]  = size(D);
	C(([1:i] + x),([1:j] + y)) = D;

	end % (nsess)

end

% contruct block (session or subject) partition - B
%---------------------------------------------------------------------------
for v = 1:nsess

	% append to B
	%-------------------------------------------------------------------
	k        = nscan(v);
	[x y]    = size(B);
	B(([1:k] + x),(1 + y)) = ones(k,1);

	% append labels
	%-------------------------------------------------------------------
	Bnames   = str2mat(Bnames,sprintf('Session %0.0f',v));
end

% confound parition - G
%---------------------------------------------------------------------------
for v = 1:nsess


	% reset name
	%-------------------------------------------------------------------
	set(Finter,'Name',sprintf('Session or subject %0.0f',v));

	% number of scans
	%-------------------------------------------------------------------
	k     = nscan(v);

	% if replications assume previous parameters
	%-------------------------------------------------------------------
	if (v == 1) | ~REP
		g = spm_input('# of confounds','!+0','e',0);
		D = [];
		while size(D,2) < g
			u   = size(D,2) + 1;
			str = sprintf('[%d]-confound %d, session %d',k,u,v);
			d   = spm_input(str,'!+0');
			if size(d,2) == k
				d = d';
			end
			if size(d,1) == k
				D = [D d];
			end
		end
	end

	% append to G
	%-------------------------------------------------------------------
	[x y] = size(G);
	[i j] = size(D);
	G(([1:i] + x),([1:j] + y)) = D;

	% append labels
	%-------------------------------------------------------------------
	for i = 1:size(D,2)
		Gnames = str2mat(Gnames,sprintf('Confound %0.0f',i));
	end

end


% high pass filter using discrete cosine set
%---------------------------------------------------------------------------
if spm_input('high pass filter','!+1','yes|no',[1 0])

	CUT   = spm_input('cut-off period {secs}','!+0','e',2*CUT*RT);
	for v = 1:nsess
		D     = [];
		k     = nscan(v);
		u     = [1:k/2];
		u     = find(u <= 2*(k*RT)/CUT);
		for i = 1:length(u)
			d = cos(pi*[0:(k - 1)]*u(i)/(k - 1));
			D = [D d(:)];
		end

    		% append to F
   		%-----------------------------------------------------------
    		[x y] = size(F);
    		[i j] = size(D);
   		F(([1:i] + x),([1:j] + y)) = D;

	end

	% mean correct
	%-------------------------------------------------------------------
	F      = spm_detrend(F);

	%-combine confounds and append confound effect names
	%-------------------------------------------------------------------
	G      = spm_en([G F]);
	for i  = 1:size(F,2);
		Gnames = str2mat(Gnames,sprintf('Low Hz (%0.0f)',i));
	end
end


% global normalization
%---------------------------------------------------------------------------
str    = 'remove global effects';
GLOBAL = spm_input(str,'!+1','scale|ancova|none',[2 1 0]);
str    = 'remove effects on response';
GXxC   = spm_input(str,'!+1','yes|no',[1 0]);


%-Construct full design matrix and name matrices for display
%---------------------------------------------------------------------------
Hnames(1,:) = [];
Cnames(1,:) = [];
Bnames(1,:) = [];

% display design matrix partition C to help contrast specification
%---------------------------------------------------------------------------
figure(Fgraph); spm_clf; axis off
axes('Position',[0.2 0.5 0.6 0.4])
imagesc(spm_DesMtxSca(C)')
set(gca,'FontSize',8)
set(gca,'YTick',[1:size(Cnames)],'YTickLabels',Cnames)
xlabel('scan')
title(['Effects of interest'],'Fontsize',16,'Fontweight','Bold')
drawnow


% time x condition interactions - over all sessions
%---------------------------------------------------------------------------
TxC = spm_input('time x response interactions','!+1','yes|no',[1 0]);

if TxC

	% number of scans
	%-------------------------------------------------------------------
	k     = sum(nscan);

	% basis functions - Type
	%-------------------------------------------------------------------
	D     = [];
	Ttype = str2mat(...
		'basis functions (Discrete Cosine Set)',...
		'Exponential decay',...
		'Linear decay',...
		'User specified');
	str   = 'Select type of response';
	Tov   = spm_input(str,'!+1','m',Ttype,[1:size(Ttype,1)]);

	% basis functions - create
	%-------------------------------------------------------------------
	if Tov == 1
		for i = 1:spm_input('number of basis functions','!+0','e',2);
			d = cos(i*pi*[0:(k - 1)]/(k - 1));
			D = [D d(:)];
		end

	elseif Tov == 2
		d   = spm_input('time constant {secs}','!+0','e',round(k/3*RT));
		D   = exp(-[0:(k - 1)]/(d/RT))';

	elseif Tov == 3
		D   = [0:(k - 1)]'/(k - 1);

	elseif Tov == 4
		% get covariates of interest
		%----------------------------------------------------------
		t     = spm_input('number of functions','!+0');
		while size(D,2) < t
			str = sprintf('[%d]-variate %d',k,size(D,2) + 1);
			d   = spm_input(str,'!+1');
			if size(d,2) == k
				d    = d';    
			end
			if size(d,1) == k
				D    = [D d];
			end	
		end
	end

	% basis functions - select effects for interaction and apply
	%-------------------------------------------------------------------
	T     = [];
	D     = spm_detrend(D);
	t     = spm_input('which responses eg 1:2','!+0','e',1);
	for i = 1:length(t)
		d = C(:,t(i));
		d = spm_detrend(d);
		d = (d*ones(1,size(D,2))).*D;
		T = [T d];
	end

	% create labels
	%-------------------------------------------------------------------
	Tnames = [];
	for i  = 1:size(T,2);
		Tnames = str2mat(Tnames,'Interactions');
	end
	Tnames(1,:) = [];

	% determine which partition to them in
	%-------------------------------------------------------------------
	Tinf  = str2mat(...
		'as effects of interest',...
		'as the only effects of interest',...
		'as confounds');
	str   = 'treat interactions';
	Tf    = spm_input(str,'!+1','m',Tinf,[1:size(Tinf,1)]);

	% asign these effects
	%-------------------------------------------------------------------
	if Tf == 1
		C      = [C T];
		Cnames = str2mat(Cnames,Tnames);

	elseif Tf == 2

		G      = [G C];
		C      = T;
		Gnames = str2mat(Gnames,Cnames);
		Cnames = Tnames;

	elseif Tf == 3

		G      = [G T];
		Gnames = str2mat(Gnames,Tnames);

	end

	% re-display design matrix partition C
	%-------------------------------------------------------------------
	figure(Fgraph); spm_clf; axis off
	axes('Position',[0.2 0.5 0.6 0.4])
	imagesc(spm_DesMtxSca(C)')
	set(gca,'FontSize',8)
	set(gca,'YTick',[1:size(Cnames)],'YTickLabels',Cnames)
	xlabel('scan')
	title(['Effects of interest'],'Fontsize',16,'Fontweight','Bold')
	drawnow

end


% get contrasts or linear compound for parameters of interest - C
%---------------------------------------------------------------------------
t   = spm_input('# of contrasts','!+1');
while size(CONTRAST,1) ~= t
	d   = [];
        str = sprintf('[%0.0f] contrast %0.0f',size(C,2),size(CONTRAST,1) + 1);
	while size(d,2) ~= size([H C],2)
		d = spm_input(str,'!+0');
	end
     	CONTRAST = [CONTRAST; d];
end

% temporal smoothing
%---------------------------------------------------------------------------
SIGMA  = spm_input('temporal smoothing fwhm-secs','!+1','e',(6 - 2*ER));
SIGMA  = SIGMA/sqrt(8*log(2))/RT;


% the interactive parts of spm_spm_ui are now finished
%---------------------------------------------------------------------------
set(Finter,'Name','thankyou','Pointer','Watch')



% get file identifiers and global values
%===========================================================================
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end


% check for consistency of image size and voxel size
%---------------------------------------------------------------------------
if ~(all(all(~diff(V([1:3],:)')))); error('images are not compatible'); end


% adjust scaling coefficients so that Grand mean - <GM> = GM
%---------------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q; GX(i) = spm_global(V(:,i)); end
V(7,:) = V(7,:)*GM/mean(GX);
GX     = GX*GM/mean(GX);

% remove global effects
%----------------------------------------------------------------------------
if GLOBAL == 1
	for v = 1:nsess

		% append to G
		%-----------------------------------------------------------
		d        = GX.*B(:,v);
		d        = spm_detrend(d);
		G        = [G d(:)];

		% append labels
		%-----------------------------------------------------------
		Gnames   = str2mat(Gnames,sprintf('Global %0.0f',v));
	end

end

% remove global x response interactions
%---------------------------------------------------------------------------
if GXxC == 1

	% interaction terms
	%-------------------------------------------------------------------
	if REP
		d     = size(C,2)/nsess;
		I     = zeros(size(C,1),d);
		for i = 1:nsess;
			I      = I + C(:,[1:d] + (i - 1)*d);
		end
	else
		I     = C;
	end
	I     = orth(spm_detrend(I));
	I     = [I.*(spm_detrend(GX)*ones(1,size(I,2)))];

	% append to G
	%-------------------------------------------------------------------
	G     = [G I]; 

	% append names
	%-------------------------------------------------------------------
	for i = 1:size(I,2)
		Gnames = str2mat(Gnames,'Glob x response');
	end
	
end

% scale if specified 
%---------------------------------------------------------------------------
if GLOBAL == 2
	V(7,:) = GM*V(7,:)./GX';
	GX = GM*GX./GX;
end

% zero pad contrast 
%---------------------------------------------------------------------------
[i j]    = size(CONTRAST);
CONTRAST = [CONTRAST zeros(i,(size([H C B G],2) - j))];


%-centre confounds (not block, which embodies the constant term)
%---------------------------------------------------------------------------
G           = spm_detrend(G);
Gnames(1,:) = [];


[nHCBG,HCBGnames] = spm_DesMtxSca(H,Hnames,C,Cnames,B,Bnames,G,Gnames);

% display analysis parameters
%===========================================================================
figure(Fgraph); spm_clf; axis off
axes('Position',[0.2 0.5 0.6 0.4])
imagesc(spm_DesMtxSca(nHCBG)')
set(gca,'FontSize',8)
set(gca,'YTick',[1:size(nHCBG)],'YTickLabels',HCBGnames)
xlabel('scan')
title(['Design Matrix'],'Fontsize',16,'Fontweight','Bold')
axes('Visible','off')

% print details
%---------------------------------------------------------------------------
text(0.00,0.40,['AnCova for ' spm_str_manip(pwd,'a24')],'Fontsize',16);
text(0.00,0.32,'Session, scan and Filename');
y     = 0.28;
dy    = 0.03;
for i = 1:size(B,2)
	u = min(find(B(:,i) > 0));
	v = max(find(B(:,i) > 0));
	d = sprintf('%d  %d to %d ',i,u,v);
	text(0,y,d,'FontSize',10)
	text(0.2,y, [P(u,:) '...'],'FontSize',8)
	y = y - dy;
end

%---------------------------------------------------------------------------
y     = y - dy;
str   = sprintf('Temporal Smoothing: {%0.1f sec FWHM}',RT*SIGMA*sqrt(8*log(2)));
text(0,y,str);
y     = y - dy;

%---------------------------------------------------------------------------
if size(F,2)
	str = sprintf('Hi-Pass cutoff period %0.0f secs',CUT);
	text(0,y,str); y = y - dy;
end

%---------------------------------------------------------------------------
text(0,y,sprintf('Grand Mean = %0.2f a.u.',GM)); y = y - dy;
text(0,y,sprintf('Interscan Interval = %0.2f secs',RT)); y = y - dy;
text(0,y,['Response Form: ' Ctype(Cov,:)  ]); y = y - dy;
text(0,y,sprintf('SPM{F} threshold p = %0.3f',UFp)); y = y - dy;

%---------------------------------------------------------------------------
if GLOBAL == 2
	text(0,y,'Data scaled to Grand Mean'); y = y - dy;
end
if GLOBAL == 1
	text(0,y,'Global effects modelled'); y = y - dy;
end
if GXxC
	text(0,y,'Global x response interactions modelled'); y = y - dy;
end
if TxC
	text(0,y,'Time x response interactions modelled'); y = y - dy;
end

% print
%---------------------------------------------------------------------------
spm_print



% implement analysis proper
%---------------------------------------------------------------------------
spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,GX,HCBGnames,P,SIGMA,RT);



%-Clear figure
%---------------------------------------------------------------------------
spm_clf(Finter); set(Finter,'Name',' ','Pointer','Arrow');
