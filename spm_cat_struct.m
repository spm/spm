function s = spm_cat_struct(s1, s2)
% Concatenates two structure arrays with possibly different fields
% FORMAT s = spm_cat_struct(s1, s2)
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cat_struct.m 5592 2013-07-24 16:25:55Z vladimir $

if isempty(s1)
    s = s2;
elseif isempty(s2)
    s = s1;
else    
    addfields1 = setdiff(fieldnames(s2), fieldnames(s1));
    for i = 1:numel(addfields1)
        [s1(:).(addfields1{i})] = deal([]);
    end
    
    addfields2 = setdiff(fieldnames(s1), fieldnames(s2));
    for i = 1:numel(addfields2)
        [s2(:).(addfields2{i})] = deal([]);
    end
    
    s = [s1(:); s2(:)];
    
    if size(s1, 1) == 1
        s = s';
    end
end