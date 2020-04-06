function S = spm_sum(varargin)
% Sum of elements
% FORMAT S = spm_sum(X,vecdim)
%
% Compatibility layer for SUM for MATLAB < R2018b
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging

% Guillaume Flandin
% $Id: spm_sum.m 7814 2020-04-06 11:47:47Z guillaume $

persistent usebuiltin
if isempty(usebuiltin)
   usebuiltin = strcmp(spm_check_version,'matlab') && ...
                spm_check_version('matlab','9.5') >= 0;
end

if usebuiltin
    S = sum(varargin{:});
else
    if nargin == 2 && isnumeric(varargin{2}) && numel(varargin{2}) > 1
        vecdim = varargin{2};
        S = sum(varargin{1},vecdim(1),varargin{3:end});
        for i=2:numel(vecdim)
            S = sum(S,vecdim(i),varargin{3:end});
        end
    else
        S = sum(varargin{:});
    end
end
