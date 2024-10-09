function path = spm_shp_get_model(name,datadir)
% Get the path to a Shape PCA model file, install/download if needed.
%
% FORMAT path = spm_get_model(name)
%
% name     - Model variable to return
%            {'subspace_scaled', 'model_variables', 'Template_{01234}'}
% dartadir - Data directory [spm('Dir')/tpl/shp]
% path     - Path to model file
%__________________________________________________________________________

% Yael Balbastre
% Copyright (C) 2024 Wellcome Centre for Human Neuroimaging

name2fname.subspace_scaled = 'subspace_scaled.nii';
name2fname.model_variables = 'model_variables.mat';
name2fname.Template_0      = 'Template_0.nii';
name2fname.Template_1      = 'Template_1.nii';
name2fname.Template_2      = 'Template_2.nii';
name2fname.Template_3      = 'Template_3.nii';
name2fname.Template_4      = 'Template_4.nii';

if nargin < 2 || isempty(datadir)
    datadir = fullfile(spm('Dir'), 'tpm', 'shp');
end

path = fullfile(datadir, name2fname.(name));

if ~exist(path, 'file')
    spm_shp_install(datadir);
end

