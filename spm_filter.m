function [vargout] = spm_filter(Action,K,Y)
% filter routine
% FORMAT [K] = spm_filter('set',K)
% FORMAT [Y] = spm_filter('apply',K,Y)
%
% Action    - 'set'   fills in filter structure K
% Action    - 'apply' applies K to Y = K*Y
% K         - filter convolution matrix or:
% K{s}      - cell of structs containing session-specific specifications
%
% K{s}.RT       - repeat time in seconds
% K{s}.row      - row of Y constituting session s
% K{s}.LChoice  - Low-pass  filtering {'hrf' 'Gaussian' 'none'}
% K{s}.LParam   - Gaussian parameter in seconds
% K{s}.HChoice  - High-pass filtering {'specify' 'none'}
% K{s}.HParam   - cut-off period in seconds
% 
% Y         - data matrix
%
% K         - filter structure
% Y         - filtered data K.K*Y
%___________________________________________________________________________
%
% spm_filter implements band pass filtering in an efficient way by
% using explicitly the projector matrix form of the High pass
% component.  spm_filter also configures the filter structure in
% accord with the specification fields if required
%___________________________________________________________________________
% %W% Karl Friston %E%


% set or apply
%---------------------------------------------------------------------------
switch Action

	case 'set'
	%-------------------------------------------------------------------
	for s = 1:length(K)

		% matrix order
		%-----------------------------------------------------------
		k     = length(K{s}.row);

		% make low pass filter
		%-----------------------------------------------------------
		switch K{s}.LChoice

			case 'none'
			%---------------------------------------------------
			K{s}.KL = speye(k);

			case 'hrf'
			%---------------------------------------------------
			h       = spm_hrf(K{s}.RT);
			h       = [h; zeros(size(h))];
			g       = abs(fft(h));
			h       = real(ifft(g));
			n       = length(h);
			h       = h([1:n/2]);
			FIL     = [h' zeros(1,k - n/2)];
			K{s}.KL = sparse(toeplitz(FIL));

			case 'Gaussian'
			%---------------------------------------------------
			sigma   = K{s}.LParam/K{s}.RT;
			FIL     = exp(-[0:(k-1)].^2/(2*sigma^2));
			FIL     = FIL.*(FIL > 1e-4);
			K{s}.KL = sparse(toeplitz(FIL));

			otherwise
			%---------------------------------------------------
			warning('Low pass Filter option unknown');

		end

		% normalize low pass
		%-----------------------------------------------------------
		K{s}.KL = spdiags(1./sum(K{s}.KL')',0,k,k)*K{s}.KL;


		% make high pass filter
		%-----------------------------------------------------------
		switch K{s}.HChoice

			case 'none'
			%---------------------------------------------------
			K{s}.KH = [];

			case 'specify'
			%---------------------------------------------------
			n       = fix(2*(k*K{s}.RT)/K{s}.HParam + 1);
			X       = spm_dctmtx(k,n);
			K{s}.KH = X(:,[2:n]);

			otherwise
			%---------------------------------------------------
			warning('High pass Filter option unknown');

		end

	end

	% return structure
	%-------------------------------------------------------------------
	vargout = K;


	case 'apply'
	%-------------------------------------------------------------------
	if iscell(K)


		% ensure requisite feild are present
		%-----------------------------------------------------------
		if ~isfield(K{1},'KL')
			K = spm_filter('set',K);
		end

		for s = 1:length(K)

			% select data
			%---------------------------------------------------
			y = Y(K{s}.row,:);

			% apply low pass filter
			%---------------------------------------------------
			y = K{s}.KL*y;

			% apply high pass filter
			%---------------------------------------------------
			if ~isempty(K{s}.KH)
				y = y - K{s}.KH*(K{s}.KH'*y);
			end

			% reset filtered data in Y
			%---------------------------------------------------
			Y(K{s}.row,:) = y;

		end

	% K is simply a convolution matrix
	%-------------------------------------------------------------------
	else
		Y = K*Y;
	end

	% return filtered data
	%-------------------------------------------------------------------
	vargout   = Y;


	otherwise
	%-------------------------------------------------------------------
	warning('Filter option unknown');


end
