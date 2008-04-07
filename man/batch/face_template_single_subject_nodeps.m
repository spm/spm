%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 185M $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio{1}.cfg_named_dir.name = 'Subject directory';
matlabbatch{1}.cfg_basicio{1}.cfg_named_dir.dirs = {'<UNDEFINED>'};
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir = '<UNDEFINED>';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent = '<UNDEFINED>';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.name = 'categorical';
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.data = {'<UNDEFINED>'};
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.quality = double(0.900000000000000022);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.sep = double(4);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.fwhm = double(5);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.rtm = double(0);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.interp = double(2);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.wrap = double([0 0 0]);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.eoptions.weight = {};
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.roptions.which = double([2 1]);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.roptions.interp = double(4);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.roptions.wrap = double([0 0 0]);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.roptions.mask = double(1);
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.roptions.prefix = 'r';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans = {'<UNDEFINED>'};
matlabbatch{5}.spmjobs{1}.temporal{1}.st.nslices = double(24);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.tr = double(2);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.ta = double(1.91666666666666674);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.so = double([24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1]);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.refslice = double(12);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.prefix = 'a';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref = '<UNDEFINED>';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.source = '<UNDEFINED>';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.other = {''};
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = double([4 2]);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = double([0.0200000000000000004 0.0200000000000000004 0.0200000000000000004 0.00100000000000000002 0.00100000000000000002 0.00100000000000000002 0.0100000000000000002 0.0100000000000000002 0.0100000000000000002 0.00100000000000000002 0.00100000000000000002 0.00100000000000000002]);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = double([7 7]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data = '<UNDEFINED>';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.GM = double([0 0 1]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.WM = double([0 0 1]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.CSF = double([0 0 0]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.biascor = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.cleanup = double(0);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.tpm = {
                                                         '/afs/fbi.ukl.uni-freiburg.de/apps/spm8/tpm/grey.nii'
                                                         '/afs/fbi.ukl.uni-freiburg.de/apps/spm8/tpm/white.nii'
                                                         '/afs/fbi.ukl.uni-freiburg.de/apps/spm8/tpm/csf.nii'
};
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.ngaus = double([2
                                                                  2
                                                                  2
                                                                  4]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.regtype = 'mni';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.warpreg = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.warpco = double(25);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.biasreg = double(0.000100000000000000005);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.biasfwhm = double(60);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.samp = double(3);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.msk = {''};
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname = '<UNDEFINED>';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample = '<UNDEFINED>';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.preserve = double(0);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.bb = double([-78 -112 -50
                                                                              78 76 85]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.vox = double([3 3 3]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.interp = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.wrap = double([0 0 0]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.prefix = 'w';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data = '<UNDEFINED>';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.fwhm = double([8 8 8]);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.dtype = double(0);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.prefix = 's';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir = '<UNDEFINED>';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.units = 'scans';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.RT = double(2);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.fmri_t = double(24);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = double(12);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans = '<UNDEFINED>';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi = '<UNDEFINED>';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg = '<UNDEFINED>';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.hpf = double(128);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.fact(1).name = 'Fam';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.fact(1).levels = double(2);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.fact(2).name = 'Rep';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.fact(2).levels = double(2);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.bases.hrf.derivs = double([1 1]);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.volt = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.global = 'None';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.mask = {''};
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.cvi = 'AR(1)';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat = '<UNDEFINED>';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.method.Classical = double(1);
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat = '<UNDEFINED>';
matlabbatch{12}.spmjobs{1}.stats{1}.con.consess{1}.fcon.name = 'Effects of interest';
matlabbatch{12}.spmjobs{1}.stats{1}.con.consess{1}.fcon.convec{1} = double([1 0 0 0 0 0 0 0 0 0 0 0
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
                                                                            0 0 0 0 0 0 0 0 0 0 0 1]);
matlabbatch{12}.spmjobs{1}.stats{1}.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{12}.spmjobs{1}.stats{1}.con.delete = double(0);
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat = '<UNDEFINED>';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.titlestr = '';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.contrasts = double(Inf);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.threshdesc = 'FWE';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.thresh = double(0.0500000000000000028);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.extent = double(0);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{13}.spmjobs{1}.stats{1}.results.print = double(1);
