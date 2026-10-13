function spm_get_test_data
% Obtains a copy of the SPM test data
%
%__________________________________________________________________________

% URL to test data
url = 'https://www.fil.ion.ucl.ac.uk/spm/download/data/tests/tests_data.zip';

% Where to extract the contents
output_file = fullfile(spm('dir'), 'tests.zip');

% Download the ZIP file
fprintf('Downloading test data...\n');
websave(output_file, url);
fprintf('Download complete\n');

% Unzip
fprintf('Unzipping...\n')
unzip(output_file, spm('dir'));
fprintf('Unzip complete\n');