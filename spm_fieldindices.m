function [ix] = spm_fieldindices(X,varargin)
% Return the indices of fields in a structure (and vice versa)
% FORMAT [i]     = spm_fieldindices(X,field1,field2,...)
% FORMAT [field] = spm_fieldindices(X,i1,i2,...)
%
% X         - structure
% field1,.. - fields
%
% i         - vector of indices or feildname{s}
%
%__________________________________________________________________________
% Copyright (C) 2010-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fieldindices.m 5769 2013-11-27 19:37:01Z karl $


% if varargin is a vector simply return fieldnames
%--------------------------------------------------------------------------
if nargin == 2
   if isnumeric(varargin{1})
       if numel(varargin{1}) > 1
           for j = 1:length(varargin{1})
               ix{j} = spm_fieldindices(X,varargin{1}(j));
           end
           return
       end
   end
end


% create structure of zeros
%--------------------------------------------------------------------------
X0    = spm_vec(X)*0;
ix    = X0;
X0    = spm_unvec(X0,X);

% and add one to specified fields
%--------------------------------------------------------------------------
for i = 1:length(varargin)
    
    if ischar(varargin{i}) && isfield(X0,varargin{i})
        
        x  = X0;
        f  = x.(varargin{i});
        f  = spm_unvec(spm_vec(f) + 1,f);
        x.(varargin{i}) = f;
        ix = ix + spm_vec(x);
    
    % or return the name of the field
    %----------------------------------------------------------------------
    elseif isnumeric(varargin{1})
        name  = fieldnames(X);
        for j = 1:length(name)
            k = spm_fieldindices(X,name{j});
            if any(ismember(varargin{i},k))
                ix = name{j};
                return
            end
        end
    end
end

% find indices
%--------------------------------------------------------------------------
ix = find(ix);
