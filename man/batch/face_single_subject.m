%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 325M $)
%-----------------------------------------------------------------------
matlabbatch{1}.cfg_basicio.cfg_named_dir.name = 'Subject directory';
matlabbatch{1}.cfg_basicio.cfg_named_dir.dirs = {{'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/'}};
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1) = cfg_dep;
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tname = 'Directory';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).tgt_spec = {};
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).sname = 'Subject directory(1)';
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{2}.cfg_basicio.cfg_cd.dir(1).src_output = substruct('.','dirs', '{}',{1});
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1) = cfg_dep;
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1).tname = 'Parent Directory';
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1).tgt_spec = {};
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1).sname = 'Subject directory(1)';
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1).src_exbranch = substruct('.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{3}.cfg_basicio.cfg_mkdir.parent(1).src_output = substruct('.','dirs', '{}',{1});
matlabbatch{3}.cfg_basicio.cfg_mkdir.name = 'categorical';
%%
matlabbatch{4}.spm.spatial.realign.estwrite.data = {
                                                    {
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
                                                    }
                                                    }';
%%
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
matlabbatch{5}.spm.temporal.st.scans{1}(1) = cfg_dep;
matlabbatch{5}.spm.temporal.st.scans{1}(1).tname = 'Session';
matlabbatch{5}.spm.temporal.st.scans{1}(1).tgt_spec = {};
matlabbatch{5}.spm.temporal.st.scans{1}(1).sname = 'Resliced Images (Sess 1)';
matlabbatch{5}.spm.temporal.st.scans{1}(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{5}.spm.temporal.st.scans{1}(1).src_output = substruct('.','sess', '()',{1}, '.','rfiles');
matlabbatch{5}.spm.temporal.st.nslices = 24;
matlabbatch{5}.spm.temporal.st.tr = 2;
matlabbatch{5}.spm.temporal.st.ta = 1.91666666666667;
matlabbatch{5}.spm.temporal.st.so = [24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1];
matlabbatch{5}.spm.temporal.st.refslice = 12;
matlabbatch{5}.spm.temporal.st.prefix = 'a';
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1) = cfg_dep;
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1).tname = 'Reference Image';
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1).tgt_spec = {};
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1).sname = 'Mean Image';
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{6}.spm.spatial.coreg.estimate.ref(1).src_output = substruct('.','rmean');
matlabbatch{6}.spm.spatial.coreg.estimate.source = {'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/Structural/sM03953_0007.img,1'};
matlabbatch{6}.spm.spatial.coreg.estimate.other = {''};
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{6}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
matlabbatch{7}.spm.spatial.preproc.data(1) = cfg_dep;
matlabbatch{7}.spm.spatial.preproc.data(1).tname = 'Data';
matlabbatch{7}.spm.spatial.preproc.data(1).tgt_spec = {};
matlabbatch{7}.spm.spatial.preproc.data(1).sname = 'Coregistered Images';
matlabbatch{7}.spm.spatial.preproc.data(1).src_exbranch = substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{7}.spm.spatial.preproc.data(1).src_output = substruct('.','cfiles');
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
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1) = cfg_dep;
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1).tname = 'Parameter File';
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1).tgt_spec = {};
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1).sname = 'Norm Params File Subj->MNI (Subj 1)';
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1).src_exbranch = substruct('.','val', '{}',{7}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.spm.spatial.normalise.write.subj.matname(1).src_output = substruct('()',{1}, '.','snfile');
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep;
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1).tname = 'Images to Write';
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1).tgt_spec = {};
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1).sname = 'Slice Timing (Sess 1)';
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1).src_exbranch = substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{8}.spm.spatial.normalise.write.subj.resample(1).src_output = substruct('()',{1}, '.','files');
matlabbatch{8}.spm.spatial.normalise.write.roptions.preserve = 0;
matlabbatch{8}.spm.spatial.normalise.write.roptions.bb = [-78 -112 -50
                                                          78 76 85];
matlabbatch{8}.spm.spatial.normalise.write.roptions.vox = [3 3 3];
matlabbatch{8}.spm.spatial.normalise.write.roptions.interp = 1;
matlabbatch{8}.spm.spatial.normalise.write.roptions.wrap = [0 0 0];
matlabbatch{8}.spm.spatial.normalise.write.roptions.prefix = 'w';
matlabbatch{9}.spm.spatial.smooth.data(1) = cfg_dep;
matlabbatch{9}.spm.spatial.smooth.data(1).tname = 'Images to Smooth';
matlabbatch{9}.spm.spatial.smooth.data(1).tgt_spec = {};
matlabbatch{9}.spm.spatial.smooth.data(1).sname = 'Normalised Images Subj 1';
matlabbatch{9}.spm.spatial.smooth.data(1).src_exbranch = substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{9}.spm.spatial.smooth.data(1).src_output = substruct('()',{1}, '.','files');
matlabbatch{9}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{9}.spm.spatial.smooth.dtype = 0;
matlabbatch{9}.spm.spatial.smooth.prefix = 's';
matlabbatch{10}.spm.stats.fmri_spec.dir(1) = cfg_dep;
matlabbatch{10}.spm.stats.fmri_spec.dir(1).tname = 'Directory';
matlabbatch{10}.spm.stats.fmri_spec.dir(1).tgt_spec = {};
matlabbatch{10}.spm.stats.fmri_spec.dir(1).sname = 'Make Directory ''categorical''';
matlabbatch{10}.spm.stats.fmri_spec.dir(1).src_exbranch = substruct('.','val', '{}',{3}, '.','val', '{}',{1});
matlabbatch{10}.spm.stats.fmri_spec.dir(1).src_output = substruct('.','dir');
matlabbatch{10}.spm.stats.fmri_spec.timing.units = 'scans';
matlabbatch{10}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t = 24;
matlabbatch{10}.spm.stats.fmri_spec.timing.fmri_t0 = 12;
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep;
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1).tname = 'Scans';
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1).tgt_spec = {};
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1).sname = 'Smoothed Images';
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1).src_exbranch = substruct('.','val', '{}',{9}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{10}.spm.stats.fmri_spec.sess.scans(1).src_output = substruct('.','files');
matlabbatch{10}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{10}.spm.stats.fmri_spec.sess.multi = {'/afs/fbi.ukl.uni-freiburg.de/projects/spm-kurs/data/face/all-conditions.mat'};
matlabbatch{10}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1) = cfg_dep;
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1).tname = 'Multiple regressors';
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1).tgt_spec = {};
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1).sname = 'Realignment Param File (Sess 1)';
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1).src_exbranch = substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{10}.spm.stats.fmri_spec.sess.multi_reg(1).src_output = substruct('.','sess', '()',{1}, '.','rpfile');
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
matlabbatch{11}.spm.stats.fmri_est.spmmat(1) = cfg_dep;
matlabbatch{11}.spm.stats.fmri_est.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{11}.spm.stats.fmri_est.spmmat(1).tgt_spec = {};
matlabbatch{11}.spm.stats.fmri_est.spmmat(1).sname = 'SPM.mat File (fMRI Design & Data)';
matlabbatch{11}.spm.stats.fmri_est.spmmat(1).src_exbranch = substruct('.','val', '{}',{10}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{11}.spm.stats.fmri_est.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{11}.spm.stats.fmri_est.method.Classical = 1;
matlabbatch{12}.spm.stats.con.spmmat(1) = cfg_dep;
matlabbatch{12}.spm.stats.con.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{12}.spm.stats.con.spmmat(1).tgt_spec = {};
matlabbatch{12}.spm.stats.con.spmmat(1).sname = 'SPM.mat File (Estimation)';
matlabbatch{12}.spm.stats.con.spmmat(1).src_exbranch = substruct('.','val', '{}',{11}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{12}.spm.stats.con.spmmat(1).src_output = substruct('.','spmmat');
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
matlabbatch{13}.spm.stats.results.spmmat(1) = cfg_dep;
matlabbatch{13}.spm.stats.results.spmmat(1).tname = 'Select SPM.mat';
matlabbatch{13}.spm.stats.results.spmmat(1).tgt_spec = {};
matlabbatch{13}.spm.stats.results.spmmat(1).sname = 'SPM.mat File (Contrast Estimation)';
matlabbatch{13}.spm.stats.results.spmmat(1).src_exbranch = substruct('.','val', '{}',{12}, '.','val', '{}',{1}, '.','val', '{}',{1});
matlabbatch{13}.spm.stats.results.spmmat(1).src_output = substruct('.','spmmat');
matlabbatch{13}.spm.stats.results.conspec.titlestr = '';
matlabbatch{13}.spm.stats.results.conspec.contrasts = Inf;
matlabbatch{13}.spm.stats.results.conspec.threshdesc = 'FWE';
matlabbatch{13}.spm.stats.results.conspec.thresh = 0.05;
matlabbatch{13}.spm.stats.results.conspec.extent = 0;
matlabbatch{13}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
matlabbatch{13}.spm.stats.results.print = 1;
