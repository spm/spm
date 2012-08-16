function fields = bf_std_fields(sel)

fields = {
    'data'
    'sources'
    'features'
    'inverse'
    'postprocessing'
    'output'
    };

if nargin > 0
    fields = fields(sel);
end
    