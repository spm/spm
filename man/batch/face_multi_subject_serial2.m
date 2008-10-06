% Collect missing inputs
subj1dir = cellstr(spm_select([1 1],'dir','Subject Dir'));
subj1epi = cellstr(spm_select([1 inf],'image','Raw EPI images'));
subj1ana = cellstr(spm_select([1 1],'image','Anatomy image'));
subj1con = cellstr(spm_select([1 1],'mat','Multiple conditions'));
subj2dir = cellstr(spm_select([1 1],'dir','Subject Dir'));
subj2epi = cellstr(spm_select([1 inf],'image','Raw EPI images'));
subj2ana = cellstr(spm_select([1 1],'image','Anatomy image'));
subj2con = cellstr(spm_select([1 1],'mat','Multiple conditions'));
% Run batch, assuming face_single_subject_template.m is in your 
% MATLAB path or working directory
% If it is not, then a full path and file name can be used instead.
spm_jobman('serial', ...
    {'face_single_subject_template.m', 'face_single_subject_template.m'}, ...
    '', subj1dir, subj1epi, subj1ana, subj1con, ...
    subj2dir, subj2epi, subj2ana, subj2con);
