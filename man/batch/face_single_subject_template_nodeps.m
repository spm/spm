%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 325M $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.cfg_named_dir.name = 'Subject directory';
matlabbatch{1}.cfg_basicio.cfg_named_dir.dirs = {'<UNDEFINED>'};
matlabbatch{2}.cfg_basicio.cfg_cd.dir = '<UNDEFINED>';
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent = '<UNDEFINED>';
matlabbatch{3}.cfg_basicio.cfg_mkdir.name = 'categorical';
matlabbatch{4}.spm.spatial.realign.estwrite.data = {'<UNDEFINED>'};
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realign.estwrite.eoptions.weight = {};
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{4}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
matlabbatch{5}.spm.temporal.st.scans = {'<UNDEFINED>'};
matlabbatch{5}.spm.temporal.st.nslices = 24;
matlabbatch{5}.spm.temporal.st.tr = 2;
matlabbatch{5}.spm.temporal.st.ta = 1.91666666666667;
matlabbatch{5}.spm.temporal.st.so = [24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
matlabbatch{5}.spm.temporal.st.refslice = 12;
matlabbatch{5}.spm.temporal.st.prefix = 'a';
matlabbatch{6}.spm.spatial.coreg.estimate.ref = '<UNDEFINED>';
matlabbatch{6}.spm.spatial.coreg.estimate.source = '<UNDEFINED>';
matlabbatch{6}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{7}.spm.spatial.preproc.data = '<UNDEFINED>';
matlabbatch{7}.spm.spatial.preproc.output.GM = [0 0 1];
matlabbatch{7}.spm.spatial.preproc.output.WM = [0 0 1];
matlabbatch{7}.spm.spatial.preproc.output.CSF = [0 0 0];
matlabbatch{7}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{7}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{7}.spm.spatial.preproc.opts.tpm = {
                                               '/export/spm-devel/spm/trunk/tpm/grey.nii'
                                               '/export/spm-devel/spm/trunk/tpm/white.nii'
                                               '/export/spm-devel/spm/trunk/tpm/csf.nii'
                                               };
matlabbatch{7}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{7}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{7}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{7}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{7}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{7}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{7}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{7}.spm.spatial.preproc.opts.msk = {''};
matlabbatch{8}.spm.spatial.normalise.write.subj.matname = '<UNDEFINED>';
matlabbatch{8}.spm.spatial.normalise.write.subj.resample = '<UNDEFINED>';
matlabbatch{8}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{8}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{8}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
matlabbatch{8}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{8}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{8}.spm.spatial.normalise.write.roptions.prefix = 'w';
matlabbatch{9}.spm.spatial.smooth.data = '<UNDEFINED>';
matlabbatch{9}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{9}.spm.spatial.smooth.dtype = 0;
matlabbatch{9}.spm.spatial.smooth.prefix = 's';
matlabbatch{10}.spm.stats.fmri_spec.dir = '<UNDEFINED>';
matlabbatch{10}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{10}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t = 24;
matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t0 = 12;
matlabbatch{10}.spm.stats.fmri_spec.sess.scans = '<UNDEFINED>';
matlabbatch{10}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{10}.spm.stats.fmri_spec.sess.multi = '<UNDEFINED>';
matlabbatch{10}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg = '<UNDEFINED>';
matlabbatch{10}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{10}.spm.stats.fmri_spec.fact(1).name = 'Fam';
matlabbatch{10}.spm.stats.fmri_spec.fact(1).levels = 2;
matlabbatch{10}.spm.stats.fmri_spec.fact(2).name = 'Rep';
matlabbatch{10}.spm.stats.fmri_spec.fact(2).levels = 2;
matlabbatch{10}.spm.stats.fmri_spec.bases.hrf.derivs = [1 1];
matlabbatch{10}.spm.stats.fmri_spec.volt = 1;
matlabbatch{10}.spm.stats.fmri_spec.global = 'None';
matlabbatch{10}.spm.stats.fmri_spec.mask = {''};
matlabbatch{10}.spm.stats.fmri_spec.cvi = 'AR(1)';
matlabbatch{11}.spm.stats.fmri_est.spmmat = '<UNDEFINED>';
matlabbatch{11}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{12}.spm.stats.con.spmmat = '<UNDEFINED>';
matlabbatch{12}.spm.stats.con.consess{1}.fcon.name = 'Effects of interest';
%%
matlabbatch{12}.spm.stats.con.consess{1}.fcon.convec = {
                                                        [1 0 0 0 0 0 0 0 0 0 0 0
                                                        0 1 0 0 0 0 0 0 0 0 0 0
                                                        0 0 1 0 0 0 0 0 0 0 0 0
                                                        0 0 0 1 0 0 0 0 0 0 0 0
                                                        0 0 0 0 1 0 0 0 0 0 0 0
                                                        0 0 0 0 0 1 0 0 0 0 0 0
                                                        0 0 0 0 0 0 1 0 0 0 0 0
                                                        0 0 0 0 0 0 0 1 0 0 0 0
                                                        0 0 0 0 0 0 0 0 1 0 0 0
                                                        0 0 0 0 0 0 0 0 0 1 0 0
                                                        0 0 0 0 0 0 0 0 0 0 1 0
                                                        0 0 0 0 0 0 0 0 0 0 0 1]
                                                        }';
%%
matlabbatch{12}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{12}.spm.stats.con.delete = 0;
matlabbatch{13}.spm.stats.results.spmmat = '<UNDEFINED>';
matlabbatch{13}.spm.stats.results.conspec.titlestr = '';
matlabbatch{13}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{13}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{13}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{13}.spm.stats.results.conspec.extent = 0;
matlabbatch{13}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{13}.spm.stats.results.print = 1;
