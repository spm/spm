function fields = bf_std_fields(sel)
%
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2013-2023 Wellcome Centre for Human Neuroimaging


fields = {
    'data'
    'sources'
    'features'
    'inverse'
    'output'
    'write'
    };

if nargin > 0
    fields = fields(sel);
end
    