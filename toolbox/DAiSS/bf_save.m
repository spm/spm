function bf_save(BF, overwrite)
% Saves BF data in a mat file
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_save.m 8223 2022-02-17 16:01:40Z george $

if nargin == 1 && exist(fullfile(pwd, 'BF.mat'), 'file')
    
    % Orignal append option appears to make the structure exponentially
    % grow with multiple calls, trying to mitgate the difference below
    %     save('BF.mat', '-struct', 'BF', '-append');
    
    % More efficient storage solution to append - GCO
    BF0 = bf_load(fullfile(pwd, 'BF.mat'));
    f0 = fields(BF0);
    f1 = fields(BF);
    
    missing = setdiff(f1,f0);
    
    if ~isempty(missing)
        for ii = 1:numel(missing)
            BF0.(missing{ii}) = BF.(missing{ii});
        end
    end
    save('BF.mat', '-struct', 'BF0', '-v7.3');
else
    save('BF.mat', '-struct', 'BF', '-v7.3');
end