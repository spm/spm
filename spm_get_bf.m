function [BF,BFstr] = spm_get_bf(W,dt)
% creates basis functions for each trial type {i} in struct BF{i}
% FORMAT [BF BFstr] = spm_get_bf(W,dt)
%
% W     - W(i) window length {time bins} for trial type i
% dt    - time bin length {seconds}
%
% BF{i} - Array of basis functions for trial type {i}
% BFstr - description of basis functions specified
%___________________________________________________________________________
%
% if W(i) > 0 then temporally extended basis functions are specified
% otherwise an event_related basis set is specified
%---------------------------------------------------------------------------
% %W% Karl Friston %E%

%---------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');

% assemble basis functions {bf}
%---------------------------------------------------------------------------
if     all(W == 0)
 

	% model event-related responses
	%-------------------------------------------------------------------
	Ctype = {
		'basis functions (Fourier set)',...
		'basis functions (Windowed Fourier set)',...
		'basis functions (Gamma functions)',...
		'basis functions (Gamma functions with derivatives)',...
		'hrf (alone)',...
		'hrf (with time derivative)',...
		'hrf (with time and dispersion derivatives)'};
	str   = 'Select basis set';
	Cov   = spm_input(str,1,'m',Ctype);
	BFstr = Ctype(Cov);


	% create basis functions
	%-------------------------------------------------------------------
	if     Cov == 1 | Cov == 2

		% Windowed (Hanning) Fourier set
		%-----------------------------------------------------------
		str   = 'window length {secs}';
		pst   = spm_input(str,2,'e',32);
		pst   = [0:dt:pst]';
		pst   = pst/max(pst);
		h     = spm_input('order',3,'e',4);

		% hanning window
		%-----------------------------------------------------------
		if Cov == 1
			g = ones(size(pst));
		else
			g = (1 - cos(2*pi*pst))/2;
		end

		% zeroth and higher terms
		%-----------------------------------------------------------
		bf    = g;
		for i = 1:h
			bf = [bf g.*sin(i*2*pi*pst)];
			bf = [bf g.*cos(i*2*pi*pst)];	
		end

	elseif Cov == 3 | Cov == 4


		% Gamma functions alone
		%-----------------------------------------------------------
		pst   = [0:dt:32]';
		dx    = 0.01;
		bf    = spm_gamma_bf(pst);

		% Gamma functions and derivatives
		%---------------------------------------------------
		if Cov == 4
			bf    = [bf (spm_gamma_bf(pst - dx) - bf)/dx];
		end


	elseif Cov == 5 | Cov == 6 | Cov == 7


		% hrf and derivatives
		%-----------------------------------------------------------
		[bf p] = spm_hrf(dt);

		% add time derivative
		%-----------------------------------------------------------
		if Cov == 6 | Cov == 7

			dp    = 1;
			p(6)  = p(6) + dp;
			D     = (bf(:,1) - spm_hrf(dt,p))/dp;
			bf    = [bf D(:)];
			p(6)  = p(6) - dp;

			% add dispersion derivative
			%---------------------------------------------------
			if Cov == 7

				dp    = 0.01;
				p(3)  = p(3) + dp;
				D     = (bf(:,1) - spm_hrf(dt,p))/dp;
				bf    = [bf D(:)];
			end
		end
	end


	% Orthogonalize and fill in basis function structure
	%-------------------------------------------------------------------
	for i = 1:length(W)
		BF{i}  =  spm_orth(bf);
	end


elseif all(W > 0)


	% covariates of interest - Type
	%-------------------------------------------------------------------
	Ctype = {
		'basis functions  (Discrete Cosine Set)',...
		'basis functions  (Mean & exponential decay)',...
		'fixed response   (Half-sine)',...
		'fixed response   (Box-car)'};
	str   = 'Select type of response';
	Cov   = spm_input(str,1,'m',Ctype);
	BFstr = Ctype(Cov);


	% convolve with HRF?
	%-------------------------------------------------------------------
	if Cov == 1
		str = 'number of basis functions';
		h   = spm_input(str,2,'e',2);
	end

	% convolve with HRF?
	%-------------------------------------------------------------------
	HRF   = spm_input('convolve with hrf',2,'b','yes|no',[1 0]);

	% ask for temporal differences
	%-------------------------------------------------------------------
	str   = 'add temporal derivatives';
	TD    = spm_input(str,3,'b','no|yes',[0 1]);
 


	% Assemble basis functions for each trial type
	%-------------------------------------------------------------------
	for i = 1:length(W)

		pst   = [1:W(i)]' - 1;
		pst   = pst/max(pst);

		% Discrete cosine set
		%-----------------------------------------------------------
		if     Cov == 1

			bf    = [];
			for j = 0:(h - 1)
				bf = [bf cos(j*pi*pst)];	
			end

		% Mean and exponential
		%-----------------------------------------------------------
		elseif Cov == 2
		
			bf    = [ones(size(pst)) exp(-pst/4)];

		% Half sine wave
		%-----------------------------------------------------------
		elseif Cov == 3

			bf    = sin(pi*pst);

		% Box car
		%-----------------------------------------------------------
		elseif Cov == 4

			bf    = ones(size(pst));

		end

		% convolve with hemodynamic response function - hrf
		%-----------------------------------------------------------
		if HRF
			hrf   = spm_hrf(dt);
			[p q] = size(bf);
			D     = [];
			for j = 1:q
				D = [D conv(bf(:,j),hrf)];
			end
			bf    = D;
		end

		% add temporal differences if specified
		%-----------------------------------------------------------
		if TD
			bf    = [bf [diff(bf); zeros(1,size(bf,2))]/dt];
		end

		% Orthogonalize and fill in Sess structure
		%-----------------------------------------------------------
		BF{i}         =  spm_orth(bf);
	end

else

	%-------------------------------------------------------------------
	str   = get(Finter,'Name'); 
	for i = 1:length(W)
		set(Finter,'Name',[str sprintf(': Trial type: %d',i)]); 

		[bf,bfstr] = spm_get_bf(W(i),dt);
		BF{i}      = bf{1};
		BFstr{i}   = bfstr{1};
	end
	set(Finter,'Name',str); 
end


% orthogonalize basis functions
%---------------------------------------------------------------------------
function bf = spm_orth(BF)
% recursive orthogonalization
% FORMAT bf = spm_orth(bf)
%___________________________________________________________________________
bf    = BF(:,1);
for i = 2:size(BF,2)
	D     = BF(:,i);
	D     = D - bf*(pinv(bf)*D);
	if any(D)
		bf = [bf D];
	end
end


% compute Gamma functions functions
%---------------------------------------------------------------------------
function bf = spm_gamma_bf(u)
% returns basis functions used for Volterra expansion
% FORMAT bf = spm_gamma_bf(u);
% u   - times {seconds}
% bf  - basis functions (mixture of Gammas)
%___________________________________________________________________________
u     = u(:);
bf    = [];
for i = 2:4
        m   = 2^i;
        s   = sqrt(m);
        bf  = [bf spm_Gpdf(u,(m/s)^2,m/s^2)];
end
