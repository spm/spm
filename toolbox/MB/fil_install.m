function varargout = fil_install(datadir)
% Download files required for Factorisation-based Image Labelling
% FORMAT [mufile,filfile] = fil_install(datadir)
%
% https://figshare.com/projects/Factorisation-based_Image_Labelling/128189
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if ~nargin
    datadir = spm_get_defaults('tbx.mb.data');
end

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

% Head Tissue Template
filename{1} = 'mu_X.nii';
url{1} = 'https://figshare.com/ndownloader/files/31699187';

% Trained FIL Model (15 training + 15 test subjects)
filename{2} = 'fil30-nuNaN-v1-d4-K24-r2-sd1.5.mat';
%filename{2} = 'fil15-nuNaN-v1-d4-K24-r3-sd2.mat';
%url{2} = 'https://figshare.com/ndownloader/files/31784579'; % Old HDF5 version
url{2} = 'https://figshare.com/ndownloader/files/39594484'; % MATLAB -v6 file format

% Deformation between ICBM-space and X
filename{3} = 'y_icbm_to_X.nii';
url{3} = 'https://figshare.com/ndownloader/files/31700369';


for i=1:numel(filename)
    varargout{i} = fullfile(datadir, filename{i});

    if ~exist(varargout{i},'file')
        try
            fprintf('Download file "%s"...', filename{i});
            websave(varargout{i}, url{i});
            fprintf(' done.\n');
        catch
            fprintf('\n');
            error('Download failed.');
        end
    end
end

