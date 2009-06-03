function ID = spm_data_id(varargin)
% generates a specific real number in a deterministic way
% from any data structure
% FORMAT ID = spm_data_id(X);
% X  - numeric, character, cell or stucture array[s]
% ID - specific ID
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak (based on Karl's spm_vec)
% $Id: spm_data_id.m 3177 2009-06-03 08:47:41Z vladimir $

X     = varargin;

if length(X) == 1
    X = X{1};
end

if ischar(X) % For now strings are not taken into account
    ID = 0;
elseif isnumeric(X) 
    Y = double(X(:));
    ID = sum(abs(Y(~isnan(Y) & ~isinf(Y)))); 
elseif isstruct(X) || isobject(X)
    X = struct(X);
    f = fieldnames(X);
    X = X(:);
    ID = 0;
    for i = 1:length(f)
        ID = ID + spm_data_id({X.(f{i})});      
    end
elseif iscell(X)
    X     = X(:);
    ID = 0;
    for i = 1:length(X)
        ID = ID + spm_data_id(X{i});
    end   
end

if isempty(ID)
    ID = 0;
end

if ID > 0
    ID = 10^-(floor(log10(ID))-2)*ID;
end