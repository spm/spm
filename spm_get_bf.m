function [bf,Bname] = spm_get_bf(dt)
% return hemodynamic basis functions
% FORMAT [bf Bname] = spm_get_bf(dt);
%
% dt    - time bin length {seconds}
%
% bf    - Matrix [hemodyanic] basis functions
% Bname - description of basis functions specified
%_______________________________________________________________________
%
% spm_get_bf prompts for basis functions to model event or epoch-related
% responses.  The basis functions returned are unitary and orthonormal
% when defined as a function of peri-stimulus time in time-bins.
% It is at this point that the distinction between event and epoch-related 
% responses enters.
%_______________________________________________________________________
% %W%  Karl Friston %E%
%
% Programmers Guide
% Batch system implemented on this routine. See spm_bch.man
% If inputs are modified in this routine, try to modify spm_bch.man
% and spm_bch_bchmat (if necessary) accordingly. 
% Calls to spm_input in this routine use the BCH gobal variable.  
%    BCH.bch_mat 
%    BCH.index0  = {'model',index_of_Analysis};
%_______________________________________________________________________

global BCH

%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)
spm_input('Hemodynamic Basis functions...',1,'d')


% assemble basis functions
%=======================================================================

% model event-related responses
%-----------------------------------------------------------------------
Ctype = {
		'hrf (alone)',...
		'hrf (with time derivative)',...
		'hrf (with time and dispersion derivatives)',...
		'basis functions (Fourier set)',...
		'basis functions (Fourier set with Hanning)',...
		'basis functions (Gamma functions)',...
		'basis functions (Gamma functions with derivatives)',...
		'basis functions (Finite Impulse Response)'};
str   = 'Select basis set';
Sel   = spm_input(str,2,'m',Ctype);
Bname = Ctype{Sel};


% create basis functions
%-----------------------------------------------------------------------
if     Sel == 4 | Sel == 5

	% Windowed (Hanning) Fourier set
	%---------------------------------------------------------------
	str   = 'window length {secs}';
	pst   = spm_input(str,3,'e',32);
	pst   = [0:dt:pst]';
	pst   = pst/max(pst);
	h     = spm_input('order',4,'e',4);


	% hanning window
	%---------------------------------------------------------------
	if Sel == 4
		g = ones(size(pst));
	else
		g = (1 - cos(2*pi*pst))/2;
	end

	% zeroth and higher terms
	%---------------------------------------------------------------
	bf    = g;
	for i = 1:h
		bf = [bf g.*sin(i*2*pi*pst)];
		bf = [bf g.*cos(i*2*pi*pst)];	
	end

elseif Sel == 6 | Sel == 7


	% Gamma functions alone
	%---------------------------------------------------------------
	pst   = [0:dt:32]';
	dx    = 0.01;
	bf    = spm_gamma_bf(pst);

	% Gamma functions and derivatives
	%---------------------------------------------------------------
	if Sel == 7
		bf  = [bf (spm_gamma_bf(pst - dx) - bf)/dx];
	end


elseif Sel == 8


	% Finite Impulse Response
	%---------------------------------------------------------------
	bin   = spm_input('bin size (seconds)',3,'e',2);	
	nb    = spm_input('number of bins',4,'e',8);
	bf    = kron(eye(nb),ones(round(bin/dt),1));


elseif Sel == 1 | Sel == 2 | Sel == 3


	% hrf and derivatives
	%---------------------------------------------------------------
	[bf p] = spm_hrf(dt);

	% add time derivative
	%---------------------------------------------------------------
	if Sel == 2 | Sel == 3

		dp    = 1;
		p(6)  = p(6) + dp;
		D     = (bf(:,1) - spm_hrf(dt,p))/dp;
		bf    = [bf D(:)];
		p(6)  = p(6) - dp;

		% add dispersion derivative
		%--------------------------------------------------------
		if Sel == 3

			dp    = 0.01;
			p(3)  = p(3) + dp;
			D     = (bf(:,1) - spm_hrf(dt,p))/dp;
			bf    = [bf D(:)];
		end
	end
end


% Orthogonalize and fill in basis function structure
%------------------------------------------------------------------------
bf    =  spm_orth(bf);


%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

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
