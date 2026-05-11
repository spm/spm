%-----------------------------------------------------------------------
% Job saved on 28-Jan-2026 16:02:43 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (00.00)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.preproc.merge.D = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.preproc.merge.recode.file = '.*';
matlabbatch{1}.spm.meeg.preproc.merge.recode.labelorg = '.*';
matlabbatch{1}.spm.meeg.preproc.merge.recode.labelnew = '#labelorg#';
matlabbatch{1}.spm.meeg.preproc.merge.prefix = 'c';
matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Merging: Merged Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.avref.fname = 'avref_montage.mat';
matlabbatch{2}.spm.meeg.preproc.prepare.task{2}.sortconditions.label = {
                                                                        'Famous'
                                                                        'Unfamiliar'
                                                                        'Scrambled'
                                                                        }';
matlabbatch{3}.spm.meeg.preproc.montage.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.montspec.montage.montagefile(1) = cfg_dep('Prepare: Average reference montage', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','avrefname'));
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.montspec.montage.keepothers = 1;
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.blocksize = 655360;
matlabbatch{3}.spm.meeg.preproc.montage.mode.write.prefix = 'M';
matlabbatch{4}.spm.meeg.other.delete.D(1) = cfg_dep('Prepare: Prepared Datafile', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
