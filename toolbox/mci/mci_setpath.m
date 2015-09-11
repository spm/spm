
% Add subdirectories to path
addpath(genpath(pwd));

% Set work_dir so as to put SPM12 and Sundials on your search path
%work_dir='C:\Users\admin\work\';
work_dir='D:\home\wpenny\';

disp('Adding SPM12 and Sundials to path ');
disp(' ');

add_dir{1}=[work_dir,'spm12'];
add_dir{2}=[work_dir,'spm12\toolbox\dcm_meeg'];
add_dir{3}=[work_dir,'spm12\toolbox\spectral'];
add_dir{4}=[work_dir,'fil\sampling\sundials-2.5.0\sundialsTB\cvodes'];
add_dir{5}=[work_dir,'fil\sampling\sundials-2.5.0\sundialsTB\cvodes\cvm'];

for i=1:length(add_dir)
    addstr=['addpath ',add_dir{i}];
    disp(addstr);
    eval(addstr);
end

