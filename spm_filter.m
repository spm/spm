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
% K{s}.HP       - low frequencies to be removed
% K{s}.LP       - sparse toepltz low-pass convolution matrix
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
			h       = 1;
			d       = 0;

			case 'hrf'
			%---------------------------------------------------
			h       = spm_hrf(K{s}.RT);
			h       = [h; zeros(size(h))];
			g       = abs(fft(h));
			h       = real(ifft(g));
			h       = fftshift(h)';
			n       = length(h);
			d       = [1:n] - n/2 - 1;

			case 'Gaussian'
			%---------------------------------------------------
			sigma   = K{s}.LParam/K{s}.RT;
			h       = exp(-[-4*sigma:4*sigma].^2/(2*sigma^2));
			n       = length(h);
			d       = [1:n] - (n + 1)/2;

			otherwise
			%---------------------------------------------------
			error('Low pass Filter option unknown');
			return

		end

		% create and normalize low pass filter
		%-----------------------------------------------------------
		K{s}.KL = spdiags(ones(k,1)*h,d,k,k);
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
