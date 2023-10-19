function BF = bf_load(file, fields)
% Load BF data into memory with just the requested fields
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2012-2023 Wellcome Centre for Human Neuroimaging


% get filename
%-------------------------------------------------------------------------
if nargin == 0
    [file, sts] = spm_select(1, '^BF.mat$', 'Select BF.mat file');
    if ~sts, BF = []; return; end
end

vars =  whos('-file', file);
if numel(intersect({vars(:).name}, bf_std_fields)) == 0
    error('No BF fields in the file');
end 

% load MAT file
%--------------------------------------------------------------------------
if nargin <= 1
    BF = load(file);
else
    if ~iscell(fields)
        fields = {fields};
    end
    BF = load(file, fields{:});
end

[sel1, sel2] = spm_match_str(bf_std_fields, fieldnames(BF));

[other_fields, sel3] = setdiff(fieldnames(BF), bf_std_fields);
[other_fields_sorted, sel4] = sort(other_fields);

tmpcell = struct2cell(BF);
BF = cell2struct(tmpcell([sel2 sel3(sel4)]), [bf_std_fields(sel1) other_fields_sorted], 1);
