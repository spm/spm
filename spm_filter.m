function [vargout] = spm_filter(K,Y)
% High and low-pass filtering
% FORMAT [Y] = spm_filter(K,Y)
% FORMAT [K] = spm_filter(K)
%
% K         - filter matrix or:
% K{s}      - cell of structs containing block-specific specifications
%
% K{s}.RT       - observation interval in seconds
% K{s}.row      - row of Y constituting block s
% K{s}.LParam   - Gaussian parameter in seconds
% K{s}.HParam   - cut-off period in seconds
%
% K{s}.HP       - low frequencies to be removed (DCT without constant)
% K{s}.LP       - sparse toeplitz low-pass convolution matrix
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
% accord with the specification fields if called with one argument
%___________________________________________________________________________
% %W% Karl Friston %E%

% set or apply
%---------------------------------------------------------------------------
if nargin == 1 & iscell(K)

	% set K.HP and/or K.LP
	%-------------------------------------------------------------------
	for s = 1:length(K)

		% block size
		%-----------------------------------------------------------
		k     = length(K{s}.row);

		% create and normalize low-pass filter
		%-----------------------------------------------------------
		if isfield(K{s},'LParam')

			sigma   = K{s}.LParam/K{s}.RT;
			h       = round(4*sigma);
			h       = exp(-[-h:h].^2/(2*sigma^2));
			n       = length(h);
			d       = [1:n] - (n + 1)/2;

			if      n == 1, h = 1; end

			K{s}.KL = spdiags(ones(k,1)*h,d,k,k);
			K{s}.KL = spdiags(1./sum(K{s}.KL')',0,k,k)*K{s}.KL;

		end

		% make high pass filter
		%-----------------------------------------------------------
		if isfield(K{s},'HParam')

			n       = fix(2*(k*K{s}.RT)/K{s}.HParam + 1);
			X       = spm_dctmtx(k,n);
			K{s}.KH = sparse(X(:,[2:n]));
		end

	end

	% return structure
	%-------------------------------------------------------------------
	vargout = K;

else
	% apply
	%-------------------------------------------------------------------
	if iscell(K)

		% ensure requisite feilds are present
		%-----------------------------------------------------------
		if ~isfield(K{1},'KL') & ~isfield(K{1},'HL')
			K = spm_filter(K);
		end

		for s = 1:length(K)

			% select data
			%---------------------------------------------------
			y = Y(K{s}.row,:);

			% apply low pass filter
			%---------------------------------------------------
			if isfield(K{s},'KL')
				y = K{s}.KL*y;
			end

			% apply high pass filter
			%---------------------------------------------------
			if isfield(K{s},'KH')
				y = y - K{s}.KH*(K{s}.KH'*y);
			end

			% reset filtered data in Y
			%---------------------------------------------------
			Y(K{s}.row,:) = y;

		end

	% K is simply a filter matrix
	%-------------------------------------------------------------------
	else
		Y = K*Y;
	end

	% return filtered data
	%-------------------------------------------------------------------
	vargout   = Y;

end
