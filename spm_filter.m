function varargout = spm_filter(K,Y)
% Removes low frequency confounds X0
% FORMAT [Y] = spm_filter(K,Y)
% FORMAT [K] = spm_filter(K)
%
% K           - filter matrix or:
% K(s)        - struct array containing partition-specific specifications
%
% K(s).RT     - observation interval in seconds
% K(s).row    - row of Y constituting block/partition s
% K(s).HParam - cut-off period in seconds
%
% K(s).X0     - low frequencies to be removed (DCT)
%
% Y           - data matrix
%
% K           - filter structure
% Y           - filtered data
%__________________________________________________________________________
%
% spm_filter implements high-pass filtering in an efficient way by
% using the residual forming matrix of X0 - low frequency confounds
%.spm_filter also configures the filter structure in accord with the
% specification fields if called with one argument
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_filter.m 4191 2011-02-03 13:30:02Z guillaume $


%-Configure filter
%==========================================================================
if nargin == 1 && isstruct(K)
    
    % set K.X0
    %----------------------------------------------------------------------
    for s = 1:length(K)
        
        % make high pass filter
        %------------------------------------------------------------------
        k       = length(K(s).row);
        n       = fix(2*(k*K(s).RT)/K(s).HParam + 1);
        X0      = spm_dctmtx(k,n);
        K(s).X0 = X0(:,2:end);
        
    end
    
    % return filter structure
    %----------------------------------------------------------------------
    varargout = { K };
    
%-Apply filter
%==========================================================================
else
    
    % K is a filter structure
    %----------------------------------------------------------------------
    if isstruct(K)
        
        % ensure requisite fields are present
        %------------------------------------------------------------------
        if ~isfield(K(1),'X0')
            K = spm_filter(K);
        end
        
        for s = 1:length(K)
            
            % select data
            %--------------------------------------------------------------
            y = Y(K(s).row,:);
            
            % apply high pass filter
            %--------------------------------------------------------------
            y = y - K(s).X0*(K(s).X0'*y);
            
            % reset filtered data in Y
            %--------------------------------------------------------------
            Y(K(s).row,:) = y;
            
        end
        
    % K is simply a filter matrix
    %----------------------------------------------------------------------
    else
        
        Y = K * Y;
        
    end
    
    % return filtered data
    %----------------------------------------------------------------------
    varargout = { Y };
end
