function [BF,BFstr] = spm_get_bf(name,T,dt,Fstr,n_s,n_c)
% creates basis functions for each trial type {i} in struct BF{i}
% FORMAT [BF BFstr] = spm_get_bf(name,T,dt,Fstr)
%
% name  - name{1 x n} name of trials or conditions
% T     - time bins per scan
% dt    - time bin length {seconds}
% Fstr  - Prompt string (usually indicates session)
%
% BF{i} - Array of basis functions for trial type {i}
% BFstr - description of basis functions specified
%_______________________________________________________________________
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%_______________________________________________________________________
% %W% Karl Friston %E%

global batch_mat;
global iA;



%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

%-Condition arguments
%-----------------------------------------------------------------------
if nargin < 6, n_c = []; end


% if no trials
%-----------------------------------------------------------------------
n      = length(name);
if ~n
	BF    = {};
	BFstr = 'none';
	return
end

% determine sort of basis functions
%-----------------------------------------------------------------------
Rtype = {'events',...
	 'epochs',...
	 'mixed'};
if n == 1
	Rtype = Rtype(1:2);
	spm_input(name{1},1,'d',Fstr,'batch',batch_mat)
else
	spm_input(Fstr,1,'d','batch',batch_mat)
end

if n > 1
	Rov   = 'mixed';
else 
	Rov   = spm_input('are these trials',2,'b',Rtype,...
	'batch',batch_mat,{'model',iA,'conditions',n_s},'types',n_c);
   Fstr='';
end

switch Rov

	% assemble basis functions {bf}
	%===============================================================
	case 'events'

	% model event-related responses
	%---------------------------------------------------------------
	Ctype = {
		'hrf (alone)',...
		'hrf (with time derivative)',...
		'hrf (with time and dispersion derivatives)',...
		'basis functions (Fourier set)',...
		'basis functions (Windowed Fourier set)',...
		'basis functions (Gamma functions)',...
		'basis functions (Gamma functions with derivatives)'};
	str   = 'Select basis set';
	Cov   = spm_input(str,2,'m',Ctype,'batch',batch_mat,...
   {'model',iA,'conditions',n_s,'bf_ev',n_c},'ev_type');
	BFstr = Ctype{Cov};


	% create basis functions
	%---------------------------------------------------------------
	if     Cov == 4 | Cov == 5

		% Windowed (Hanning) Fourier set
		%-------------------------------------------------------
		str   = 'window length {secs}';
		pst   = spm_input(str,3,'e',32,'batch',batch_mat,...
      			{'model',iA,'conditions',n_s,'bf_ev',n_c},'win_len');
		pst   = [0:dt:pst]';
		pst   = pst/max(pst);
		h     = spm_input('order',4,'e',4,'batch',batch_mat,...
      			{'model',iA,'conditions',n_s,'bf_ev',n_c},'order');


		% hanning window
		%-------------------------------------------------------
		if Cov == 4
			g = ones(size(pst));
		else
			g = (1 - cos(2*pi*pst))/2;
		end

		% zeroth and higher terms
		%-------------------------------------------------------
		bf    = g;
		for i = 1:h
			bf = [bf g.*sin(i*2*pi*pst)];
			bf = [bf g.*cos(i*2*pi*pst)];	
		end

	elseif Cov == 6 | Cov == 7


		% Gamma functions alone
		%-------------------------------------------------------
		pst   = [0:dt:32]';
		dx    = 0.01;
		bf    = spm_gamma_bf(pst);

		% Gamma functions and derivatives
		%-------------------------------------------------------
		if Cov == 7
			bf  = [bf (spm_gamma_bf(pst - dx) - bf)/dx];
		end


	elseif Cov == 1 | Cov == 2 | Cov == 3


		% hrf and derivatives
		%-------------------------------------------------------
		[bf p] = spm_hrf(dt);

		% add time derivative
		%-------------------------------------------------------
		if Cov == 2 | Cov == 3

			dp    = 1;
			p(6)  = p(6) + dp;
			D     = (bf(:,1) - spm_hrf(dt,p))/dp;
			bf    = [bf D(:)];
			p(6)  = p(6) - dp;

			% add dispersion derivative
			%-----------------------------------------------
			if Cov == 3

				dp    = 0.01;
				p(3)  = p(3) + dp;
				D     = (bf(:,1) - spm_hrf(dt,p))/dp;
				bf    = [bf D(:)];
			end
		end
	end


	% Orthogonalize and fill in basis function structure
	%---------------------------------------------------------------
	bf    =  spm_orth(bf);
	for i = 1:n
		BF{i}  =  bf;
	end


	% assemble basis functions {bf}
	%===============================================================
	case 'epochs'


	% covariates of interest - Type
	%---------------------------------------------------------------
	Ctype = {'basis functions  (Discrete Cosine Set)',...
		 'basis functions  (Mean & exponential decay)',...
		 'fixed response   (Half-sine)',...
		 'fixed response   (Box-car)'};
	str   = 'Select type of response';
	Cov   = spm_input(str,2,'m',Ctype,'batch',batch_mat,...
   {'model',iA,'conditions',n_s,'bf_ep',n_c},'ep_type');

   BFstr = Ctype{Cov};


	% convolve with HRF?
	%---------------------------------------------------------------
	if Cov == 1
		str = 'number of basis functions';
		h   = spm_input(str,3,'e',2,'batch',batch_mat,...
  		{'model',iA,'conditions',n_s,'bf_ep',n_c},'fct_nb');
	end

	% convolve with HRF?
	%---------------------------------------------------------------
	HRF   = spm_input('convolve with hrf',3,'b','yes|no',[1 0],...
   'batch',batch_mat,{'model',iA,'conditions',n_s,'bf_ep',n_c},'conv');

	% ask for temporal differences
	%---------------------------------------------------------------
	str   = 'add temporal derivatives';
	TD    = spm_input(str,4,'b','yes|no',[1 0],'batch',batch_mat,...
  	{'model',iA,'conditions',n_s,'bf_ep',n_c},'deriv');
 

	% Assemble basis functions for each trial type
	%---------------------------------------------------------------
	for i = 1:n

		str   = ['epoch length {scans} for ' name{i}];
		W     = spm_input(str,'+1','r','batch',batch_mat,...
 	   {'model',iA,'conditions',n_s,'bf_ep',n_c},'length');
		pst   = [1:W*T]' - 1;
		pst   = pst/max(pst);

		% Discrete cosine set
		%-------------------------------------------------------
		if     Cov == 1

			bf    = [];
			for j = 0:(h - 1)
				bf = [bf cos(j*pi*pst)];	
			end

		% Mean and exponential
		%-------------------------------------------------------
		elseif Cov == 2
		
			bf    = [ones(size(pst)) exp(-pst/4)];

		% Half sine wave
		%-------------------------------------------------------
		elseif Cov == 3

			bf    = sin(pi*pst);

		% Box car
		%-------------------------------------------------------
		elseif Cov == 4

			bf    = ones(size(pst));

		end

		% convolve with hemodynamic response function - hrf
		%-------------------------------------------------------
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
		%-------------------------------------------------------
		if TD
			bf    = [bf [diff(bf); zeros(1,size(bf,2))]/dt];
		end

		% Orthogonalize and fill in Sess structure
		%-------------------------------------------------------
		BF{i}         =  spm_orth(bf);
	end


	% mixed event and epoch model
	%===============================================================
	case 'mixed'
	for i = 1:n
		
		BF(i)  = spm_get_bf(name(i),T,dt,'',n_s,i);

	end
	BFstr = 'mixed';

end



%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================


function bf = spm_orth(BF)
% recursive orthogonalization of basis functions
% FORMAT bf = spm_orth(bf)
%_______________________________________________________________________
bf    = BF(:,1);
bf    = bf/sqrt(mean(bf.^2));
for i = 2:size(BF,2)
	D     = BF(:,i);
	D     = D - bf*(pinv(bf)*D);
	if any(D)
		bf = [bf D/sqrt(mean(D.^2))];
	end
end


% compute Gamma functions functions
%-----------------------------------------------------------------------
function bf = spm_gamma_bf(u)
% returns basis functions used for Volterra expansion
% FORMAT bf = spm_gamma_bf(u);
% u   - times {seconds}
% bf  - basis functions (mixture of Gammas)
%_______________________________________________________________________
u     = u(:);
bf    = [];
for i = 2:4
        m   = 2^i;
        s   = sqrt(m);
        bf  = [bf spm_Gpdf(u,(m/s)^2,m/s^2)];
end
