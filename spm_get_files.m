function Q = spm_get_files(dir,fil)
% FORMAT Q spm_get_files(dir,fil)
% dir    - string : the directory name
% fil    - string : filter, like 'sn*.img'
%
% Based on spm_list_files.m
% used in m-files for batch mode to help the users
% in entering the file names.
% Not used by any spm routine.
% %W% %E%

[P,tmp] = spm_list_files(dir,fil);
for k = 1:size(P,1)
  Q(k,:) = fullfile(dir,P(k,:));
end

