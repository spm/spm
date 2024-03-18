function bf_save_path(BF,path, overwrite)
% Saves BF data in a MAT file
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2012-2023 Wellcome Centre for Human Neuroimaging


if nargin == 2 && exist(path, 'file')
    % Orignal append option appears to make the structure exponentially
    % grow with multiple calls, trying to mitgate the difference below
    %     save(path, '-struct', 'BF', '-append');
    
    % More efficient storage solution to append - GCO
    BF0 = bf_load(path);
    f1 = fields(BF);
    
    % overwrite fields
    for ii = 1:numel(f1)
        BF0.(f1{ii}) = BF.(f1{ii});
    end
    save(path, '-struct', 'BF0', '-v7.3');
    
else
    save(path, '-struct', 'BF', '-v7.3');
end
