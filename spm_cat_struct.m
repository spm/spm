function s = spm_cat_struct(s1, s2, varargin)
% Concatenates structure arrays with possibly different fields
% FORMAT s = spm_cat_struct(s1, s2, ...)
%__________________________________________________________________________
% Copyright (C) 2013-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cat_struct.m 7074 2017-05-05 10:51:08Z guillaume $


if ~nargin
    s1 = struct([]);
    s2 = struct([]);
elseif nargin == 1
    s2 = struct([]);
elseif nargin > 2
    s2 = spm_cat_struct(s2,varargin{:});
end

if isempty(s1)
    s = s2;
elseif isempty(s2)
    s = s1;
else
    addfields1 = setdiff(fieldnames(s2), fieldnames(s1));
    for i = 1:numel(addfields1)
        [s1.(addfields1{i})] = deal([]);
    end
    
    addfields2 = setdiff(fieldnames(s1), fieldnames(s2));
    for i = 1:numel(addfields2)
        [s2.(addfields2{i})] = deal([]);
    end
    
    s = [s1(:); s2(:)];
    
    if size(s1, 1) == 1
        s = s';
    end
end
