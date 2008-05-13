%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 264M $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio{1}.runjobs.jobs = {fullfile(spm('dir'),'man','batch','face_template_single_subject.m')};
matlabbatch{1}.cfg_basicio{1}.runjobs.inputs{1}{1}.indir = '<UNDEFINED>';
matlabbatch{1}.cfg_basicio{1}.runjobs.inputs{1}{2}.innifti = '<UNDEFINED>';
matlabbatch{1}.cfg_basicio{1}.runjobs.inputs{1}{3}.innifti = '<UNDEFINED>';
matlabbatch{1}.cfg_basicio{1}.runjobs.inputs{1}{4}.inmat = '<UNDEFINED>';
matlabbatch{1}.cfg_basicio{1}.runjobs.save.savejobs.outstub = 'face_single_subject_run';
matlabbatch{1}.cfg_basicio{1}.runjobs.save.savejobs.outdir = {fullfile(spm('dir'),'man','batch')};
matlabbatch{1}.cfg_basicio{1}.runjobs.missing = 'skip';
