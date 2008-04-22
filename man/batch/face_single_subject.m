%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 187 $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio{1}.cfg_named_dir.name = 'Subject directory';
matlabbatch{1}.cfg_basicio{1}.cfg_named_dir.dirs{1} = {'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/'};
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1) = cfg_dep;
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tname = 'Directory';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_spec{1}(1).name = 'class';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).tgt_spec{1}(2).value = 'e';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).sname = 'Subject directory(1)';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(1).type = '.';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(1).subs = 'val';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(2).type = '{}';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(2).subs{1} = double(1);
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(3).type = '.';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(3).subs = 'val';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(4).type = '{}';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_output(1).type = '.';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_output(1).subs = 'dirs';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_output(2).type = '{}';
matlabbatch{2}.cfg_basicio{1}.cfg_cd.dir(1).src_output(2).subs{1} = double(1);
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1) = cfg_dep;
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tname = 'Parent Directory';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_spec{1}(1).name = 'class';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).tgt_spec{1}(2).value = 'e';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).sname = 'Subject directory(1)';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(1).type = '.';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(1).subs = 'val';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(2).type = '{}';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(2).subs{1} = double(1);
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(3).type = '.';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(3).subs = 'val';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(4).type = '{}';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_output(1).type = '.';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_output(1).subs = 'dirs';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_output(2).type = '{}';
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.parent(1).src_output(2).subs{1} = double(1);
matlabbatch{3}.cfg_basicio{1}.cfg_mkdir.name = 'categorical';
matlabbatch{4}.spmjobs{1}.spatial{1}.realign{1}.estwrite.data{1} = {
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0006.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0007.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0008.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0009.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0010.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0011.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0012.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0013.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0014.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0015.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0016.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0017.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0018.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0019.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0020.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0021.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0022.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0023.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0024.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0025.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0026.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0027.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0028.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0029.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0030.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0031.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0032.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0033.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0034.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0035.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0036.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0037.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0038.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0039.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0040.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0041.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0042.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0043.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0044.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0045.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0046.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0047.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0048.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0049.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0050.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0051.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0052.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0053.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0054.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0055.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0056.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0057.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0058.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0059.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0060.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0061.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0062.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0063.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0064.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0065.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0066.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0067.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0068.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0069.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0070.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0071.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0072.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0073.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0074.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0075.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0076.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0077.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0078.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0079.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0080.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0081.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0082.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0083.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0084.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0085.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0086.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0087.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0088.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0089.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0090.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0091.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0092.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0093.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0094.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0095.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0096.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0097.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0098.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0099.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0100.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0101.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0102.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0103.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0104.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0105.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0106.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0107.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0108.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0109.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0110.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0111.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0112.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0113.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0114.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0115.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0116.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0117.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0118.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0119.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0120.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0121.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0122.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0123.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0124.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0125.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0126.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0127.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0128.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0129.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0130.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0131.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0132.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0133.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0134.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0135.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0136.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0137.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0138.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0139.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0140.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0141.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0142.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0143.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0144.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0145.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0146.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0147.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0148.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0149.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0150.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0151.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0152.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0153.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0154.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0155.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0156.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0157.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0158.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0159.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0160.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0161.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0162.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0163.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0164.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0165.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0166.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0167.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0168.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0169.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0170.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0171.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0172.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0173.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0174.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0175.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0176.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0177.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0178.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0179.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0180.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0181.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0182.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0183.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0184.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0185.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0186.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0187.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0188.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0189.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0190.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0191.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0192.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0193.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0194.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0195.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0196.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0197.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0198.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0199.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0200.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0201.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0202.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0203.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0204.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0205.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0206.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0207.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0208.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0209.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0210.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0211.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0212.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0213.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0214.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0215.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0216.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0217.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0218.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0219.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0220.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0221.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0222.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0223.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0224.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0225.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0226.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0227.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0228.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0229.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0230.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0231.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0232.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0233.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0234.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0235.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0236.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0237.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0238.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0239.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0240.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0241.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0242.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0243.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0244.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0245.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0246.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0247.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0248.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0249.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0250.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0251.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0252.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0253.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0254.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0255.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0256.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0257.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0258.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0259.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0260.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0261.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0262.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0263.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0264.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0265.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0266.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0267.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0268.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0269.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0270.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0271.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0272.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0273.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0274.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0275.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0276.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0277.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0278.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0279.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0280.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0281.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0282.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0283.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0284.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0285.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0286.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0287.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0288.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0289.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0290.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0291.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0292.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0293.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0294.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0295.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0296.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0297.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0298.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0299.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0300.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0301.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0302.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0303.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0304.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0305.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0306.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0307.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0308.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0309.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0310.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0311.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0312.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0313.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0314.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0315.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0316.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0317.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0318.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0319.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0320.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0321.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0322.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0323.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0324.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0325.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0326.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0327.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0328.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0329.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0330.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0331.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0332.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0333.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0334.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0335.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0336.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0337.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0338.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0339.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0340.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0341.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0342.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0343.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0344.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0345.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0346.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0347.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0348.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0349.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0350.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0351.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0352.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0353.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0354.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0355.img,1'
                                                                    '/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/RawEPI/sM03953_0005_0356.img,1'
};
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
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1) = cfg_dep;
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tname = 'Session';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_spec{1}(1).name = 'class';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).tgt_spec{1}(2).value = 'e';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).sname = 'Resliced Images (Sess 1)';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(1).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(1).subs = 'val';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(2).type = '{}';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(2).subs{1} = double(4);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(3).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(3).subs = 'val';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(4).type = '{}';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(5).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(5).subs = 'val';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(6).type = '{}';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(7).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(7).subs = 'val';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(8).type = '{}';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_exbranch(8).subs{1} = double(1);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(1).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(1).subs = 'sess';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(2).type = '()';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(2).subs{1} = double(1);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(3).type = '.';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.scans{1}(1).src_output(3).subs = 'rfiles';
matlabbatch{5}.spmjobs{1}.temporal{1}.st.nslices = double(24);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.tr = double(2);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.ta = double(1.91666666666666674);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.so = double([24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1]);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.refslice = double(12);
matlabbatch{5}.spmjobs{1}.temporal{1}.st.prefix = 'a';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1) = cfg_dep;
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tname = 'Reference Image';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_spec{1}(1).name = 'class';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).tgt_spec{1}(2).value = 'e';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).sname = 'Mean Image';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(1).type = '.';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(1).subs = 'val';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(2).type = '{}';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(2).subs{1} = double(4);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(3).type = '.';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(3).subs = 'val';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(4).type = '{}';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(5).type = '.';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(5).subs = 'val';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(6).type = '{}';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(7).type = '.';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(7).subs = 'val';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(8).type = '{}';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_exbranch(8).subs{1} = double(1);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_output(1).type = '.';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.ref(1).src_output(1).subs = 'rmean';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.source = {'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/Structural/sM03953_0007.img,1'};
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.other = {''};
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = double([4 2]);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = double([0.0200000000000000004 0.0200000000000000004 0.0200000000000000004 0.00100000000000000002 0.00100000000000000002 0.00100000000000000002 0.0100000000000000002 0.0100000000000000002 0.0100000000000000002 0.00100000000000000002 0.00100000000000000002 0.00100000000000000002]);
matlabbatch{6}.spmjobs{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = double([7 7]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1) = cfg_dep;
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tname = 'Data';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_spec{1}(1).name = 'class';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).sname = 'Coregistered Images';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(1).type = '.';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(1).subs = 'val';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(2).type = '{}';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(2).subs{1} = double(6);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(3).type = '.';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(3).subs = 'val';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(4).type = '{}';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(5).type = '.';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(5).subs = 'val';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(6).type = '{}';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(7).type = '.';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(7).subs = 'val';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(8).type = '{}';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_exbranch(8).subs{1} = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_output(1).type = '.';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.data(1).src_output(1).subs = 'cfiles';
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.GM = double([0 0 1]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.WM = double([0 0 1]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.CSF = double([0 0 0]);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.biascor = double(1);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.output.cleanup = double(0);
matlabbatch{7}.spmjobs{1}.spatial{1}.preproc.opts.tpm = {
  fullfile(spm('dir'),'tpm','grey.nii')
  fullfile(spm('dir'),'tpm','white.nii')
  fullfile(spm('dir'),'tpm','csf.nii')
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
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1) = cfg_dep;
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_spec{1}(1).name = 'class';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).tgt_spec{1}(2).value = 'e';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).sname = 'Norm Params File Subj->MNI (Subj 1)';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(1).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(1).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(2).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(2).subs{1} = double(7);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(3).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(3).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(4).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(5).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(5).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(6).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_output(1).type = '()';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_output(1).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_output(2).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.matname(1).src_output(2).subs = 'snfile';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1) = cfg_dep;
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_spec{1}(1).name = 'class';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).tgt_spec{1}(2).value = 'e';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).sname = 'Slice Timing (Sess 1)';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(1).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(1).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(2).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(2).subs{1} = double(5);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(3).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(3).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(4).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(5).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(5).subs = 'val';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(6).type = '{}';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_output(1).type = '()';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_output(1).subs{1} = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_output(2).type = '.';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.subj.resample(1).src_output(2).subs = 'files';
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.preserve = double(0);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.bb = double([-78 -112 -50
                                                                              78 76 85]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.vox = double([3 3 3]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.interp = double(1);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.wrap = double([0 0 0]);
matlabbatch{8}.spmjobs{1}.spatial{1}.normalise{1}.write.roptions.prefix = 'w';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1) = cfg_dep;
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tname = 'Images to Smooth';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_spec{1}(1).name = 'class';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).tgt_spec{1}(2).value = 'e';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).sname = 'Normalised Images Subj 1';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(1).type = '.';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(1).subs = 'val';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(2).type = '{}';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(2).subs{1} = double(8);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(3).type = '.';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(3).subs = 'val';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(4).type = '{}';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(5).type = '.';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(5).subs = 'val';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(6).type = '{}';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(7).type = '.';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(7).subs = 'val';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(8).type = '{}';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_exbranch(8).subs{1} = double(1);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_output(1).type = '()';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_output(1).subs{1} = double(1);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_output(2).type = '.';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.data(1).src_output(2).subs = 'files';
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.fwhm = double([8 8 8]);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.dtype = double(0);
matlabbatch{9}.spmjobs{1}.spatial{1}.smooth.prefix = 's';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1) = cfg_dep;
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tname = 'Directory';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_spec{1}(1).name = 'filter';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_spec{1}(1).value = 'dir';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).tgt_spec{1}(2).value = 'e';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).sname = 'Make Directory ''categorical''';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(1).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(2).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(2).subs{1} = double(3);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(3).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(3).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(4).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_output(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.dir(1).src_output(1).subs = 'dir';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.units = 'scans';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.RT = double(2);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.fmri_t = double(24);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.timing.fmri_t0 = double(12);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1) = cfg_dep;
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tname = 'Scans';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_spec{1}(1).name = 'class';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).tgt_spec{1}(2).value = 'e';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).sname = 'Smoothed Images';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(1).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(2).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(2).subs{1} = double(9);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(3).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(3).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(4).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(5).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(5).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(6).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_output(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.scans(1).src_output(1).subs = 'files';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi = {'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/all-conditions.mat'};
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1) = cfg_dep;
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tname = 'Multiple regressors';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).name = 'class';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).tgt_spec{1}(2).value = 'e';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).sname = 'Realignment Param File (Sess 1)';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(1).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(2).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(2).subs{1} = double(4);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(3).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(3).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(4).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(5).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(5).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(6).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(7).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(7).subs = 'val';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(8).type = '{}';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_exbranch(8).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(1).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(1).subs = 'sess';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(2).type = '()';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(2).subs{1} = double(1);
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(3).type = '.';
matlabbatch{10}.spmjobs{1}.stats{1}.fmri_spec.sess.multi_reg(1).src_output(3).subs = 'rpfile';
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
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_spec{1}(1).name = 'class';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).sname = 'SPM.mat File (fMRI Design & Data)';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(1).type = '.';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(1).subs = 'val';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(2).type = '{}';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(2).subs{1} = double(10);
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(3).type = '.';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(3).subs = 'val';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(4).type = '{}';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(5).type = '.';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(5).subs = 'val';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(6).type = '{}';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_output(1).type = '.';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.spmmat(1).src_output(1).subs = 'spmmat';
matlabbatch{11}.spmjobs{1}.stats{1}.fmri_est.method.Classical = double(1);
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1) = cfg_dep;
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_spec{1}(1).name = 'class';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).sname = 'SPM.mat File (Estimation)';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(1).type = '.';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(1).subs = 'val';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(2).type = '{}';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(2).subs{1} = double(11);
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(3).type = '.';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(3).subs = 'val';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(4).type = '{}';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(5).type = '.';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(5).subs = 'val';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(6).type = '{}';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_output(1).type = '.';
matlabbatch{12}.spmjobs{1}.stats{1}.con.spmmat(1).src_output(1).subs = 'spmmat';
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
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1) = cfg_dep;
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_exbranch = struct('type', {}, 'subs', {});
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_input = struct('type', {}, 'subs', {});
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_spec{1}(1).name = 'class';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_spec{1}(1).value = 'cfg_files';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_spec{1}(2).name = 'strtype';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).tgt_spec{1}(2).value = 'e';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).jtsubs = struct('type', {}, 'subs', {});
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).sname = 'SPM.mat File (Contrast Estimation)';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(1).type = '.';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(1).subs = 'val';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(2).type = '{}';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(2).subs{1} = double(12);
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(3).type = '.';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(3).subs = 'val';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(4).type = '{}';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(4).subs{1} = double(1);
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(5).type = '.';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(5).subs = 'val';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(6).type = '{}';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_exbranch(6).subs{1} = double(1);
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_output(1).type = '.';
matlabbatch{13}.spmjobs{1}.stats{1}.results.spmmat(1).src_output(1).subs = 'spmmat';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.titlestr = '';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.contrasts = double(Inf);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.threshdesc = 'FWE';
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.thresh = double(0.0500000000000000028);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.extent = double(0);
matlabbatch{13}.spmjobs{1}.stats{1}.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{13}.spmjobs{1}.stats{1}.results.print = double(1);
