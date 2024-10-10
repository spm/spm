function datadir = spm_shp_install(datadir,force)
% Download files required for Shape PCA model
%
% FORMAT datadir = spm_shp_install(datadir,[force=true])
%__________________________________________________________________________

% John Ashburner, Yael Balbastre
% Copyright (C) 2023-2024 Wellcome Centre for Human Neuroimaging

default_datadir = fullfile(spm('Dir'), 'tpm', 'shp');

if nargin < 2,                     force   = false;             end
if nargin < 1 || isempty(datadir), datadir = default_datadir;   end

if ~exist(datadir,'dir')
    try
        fprintf('Create directory "%s"...', datadir);
        mkdir(datadir);
        fprintf(' done.\n');
    catch
        fprintf('\n');
        error(['Failed to make a directory for the template and labelling files (' datadir ').']);
    end
end

filename{1} = 'subspace_scaled.nii';
url{1}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535025';
filename{2} = 'model_variables.mat';
url{2}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535022';
filename{3} = 'Template_0.nii';
url{3}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49534998';
filename{4} = 'Template_1.nii';
url{4}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535001';
filename{5} = 'Template_2.nii';
url{5}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535004';
filename{6} = 'Template_3.nii';
url{6}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535007';
filename{7} = 'Template_4.nii';
url{7}      = 'https://rdr.ucl.ac.uk/ndownloader/files/49535013';

for i=1:numel(filename)
    path_out = fullfile(datadir, filename{i});

    if force || ~exist(path_out,'file')
        try
            fprintf('Download file "%s"...', filename{i});
            websave(path_out, url{i});
            fprintf(' done.\n');
        catch
            fprintf('\n');
            error('Download failed.');
        end
    end
end

