function fields = bf_std_fields(sel)

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
    