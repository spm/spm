function [i] = spm_fieldindices(X,varargin)
% returns the indices of fields in a structure
% FORMAT [i] = spm_fieldindices(X,feild1,feild2,...)
%
% X         - structure
% feild1,.. – fields
%
% i         - vector of indices
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fieldindices.m 4100 2010-10-22 19:49:17Z karl $
 
 
% create structure of zeros
%--------------------------------------------------------------------------
X0    = spm_vec(X)*0;
ix    = X0;
X0    = spm_unvec(X0,X);
 
% and add one to specified fields
%--------------------------------------------------------------------------
for i = 1:length(varargin)
    
    x  = X0;
    f  = getfield(x,varargin{i});
    f  = spm_unvec(spm_vec(f) + 1,f);
    x  = setfield(x,varargin{i},f);
    ix = ix + spm_vec(x);
    
end
 
% find indices
%--------------------------------------------------------------------------
i = find(ix);
