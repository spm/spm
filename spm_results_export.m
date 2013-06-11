function spm_results_export(SPM,xSPM,TabDat)
% Export SPM results in HTML
% FORMAT spm_results_export(SPM,xSPM,TabDat)
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_results_export.m 5543 2013-06-11 17:48:18Z guillaume $


if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3
    TabDat = spm_list('Table',xSPM);
end

MIP     = spm_mip(xSPM.Z,xSPM.XYZmm,xSPM.M,xSPM.units);
DesMtx  = (SPM.xX.nKX + 1)*32;

outdir  = pwd;
fHTML   = spm_file(fullfile(outdir,'index.html'),'unique');
fMIP    = spm_file(fullfile(outdir,'MIP.png'),'unique');
fDesMtx = spm_file(fullfile(outdir,'DesMtx.png'),'unique');
fcon    = spm_file(fullfile(outdir,'contrast.png'),'unique');
fcursor = spm_file(fullfile(outdir,'cursor.png'),'unique');

%-Save images as PNG files
%--------------------------------------------------------------------------
imwrite(MIP,gray(64),fMIP,'png');

ml = floor(size(DesMtx,1)/size(DesMtx,2));
DesMtx = reshape(repmat(DesMtx,ml,1),size(DesMtx,1),[]);
imwrite(DesMtx,gray(64),fDesMtx,'png');

con = [SPM.xCon(xSPM.Ic).c]';
con = (con/max(abs(con(:)))+1)*32;
con = kron(con,ones(ml,10));
imwrite(con,gray(64),fcon,'png');

cursor = zeros(11,8,3);
cursor([5 6 7 16 17 18 26 27 29 30 37 38 40 41 47 48 52 53 57 58 64 65 68 69 75 76 78 79 87 88]) = 1;
imwrite(cursor,fcursor,'png','Transparency',[0 0 0]);

%-Save results as HTML file
%==========================================================================
tpl = spm_file_template(fullfile(spm('Dir'),'help'),'keep');
tpl = tpl.file('TPL_RES','spm_results.tpl');
tpl = tpl.block('TPL_RES','resftr','resftrs');
tpl = tpl.block('TPL_RES','cursor','cursors');
tpl = tpl.block('TPL_RES','resrow','resrows');
tpl = tpl.var('DATE',[datestr(now,8) ' ' datestr(now,1) ' ' datestr(now,13)]);
tpl = tpl.var('SPM',spm('Version'));
tpl = tpl.var('CON_TITLE',xSPM.title);
tpl = tpl.var('RES_TITLE',TabDat.tit);
tpl = tpl.var('IMG_MIP',spm_file(fMIP,'filename'));
tpl = tpl.var('IMG_CURSOR',spm_file(fcursor,'filename'));
tpl = tpl.var('IMG_CON',spm_file(fcon,'filename'));
tpl = tpl.var('IMG_X',spm_file(fDesMtx,'filename'));
tpl = tpl.var('cursors','');
XYZmm = [0 0 0];
XYZ = spm_XYZreg('RoundCoords',XYZmm,xSPM.M,xSPM.DIM);
%-Copied from spm_mip_ui.m
DXYZ = [182 218 182];
CXYZ = [091 127 073];
%-Coordinates of Talairach origin in MIP image
% Transverse: [Po(1), Po(2)]
% Saggital  : [Po(1), Po(3)]
% Coronal   : [Po(4), Po(3)]
% 4 voxel offsets in Y since using image '<' as a pointer.
Po(1)  =                  CXYZ(2) -2;
Po(2)  = DXYZ(3)+DXYZ(1) -CXYZ(1) -4;
Po(3)  = DXYZ(3)         -CXYZ(3) -4;
Po(4)  = DXYZ(2)         +CXYZ(1) -2;
x(1)=Po(1)+XYZ(2); y(1)=Po(2)+XYZ(1);
x(2)=Po(1)+XYZ(2); y(2)=Po(3)-XYZ(3);
x(3)=Po(4)+XYZ(1); y(3)=Po(3)-XYZ(3);
for i=1:3
    tpl = tpl.var('CS_ID',sprintf('cs%d',i));
    tpl = tpl.var('CS_X',sprintf('%d',x(i)));
    tpl = tpl.var('CS_Y',sprintf('%d',y(i)));
    tpl = tpl.parse('cursors','cursor',1);
end
tpl = tpl.var('RES_STR',TabDat.str);
tpl = tpl.var('STAT_STR',strrep(strrep(xSPM.STATstr,'_{','<sub>'),'}','</sub>'));
tpl = tpl.var('CON_STAT',xSPM.STAT);
tpl = tpl.var('resrows','');
for i=1:size(TabDat.dat,1)
    tpl = tpl.var('RES_COL1','1');
    for j=1:size(TabDat.dat,2)-1
        tpl = tpl.var(sprintf('RES_COL%d',j),sprintf(TabDat.fmt{j},TabDat.dat{i,j}));
    end
    tpl = tpl.var('RES_COL12',strrep(sprintf(TabDat.fmt{end},TabDat.dat{i,end}),' ','&nbsp;'));
    tpl = tpl.var('RES_XYZ',sprintf('%d,%d,%d',TabDat.dat{i,end}));
    tpl = tpl.parse('resrows','resrow',1);
end
tpl = tpl.var('resftrs','');
for i=1:size(TabDat.ftr,1)
	tpl = tpl.var('RES_FTR',sprintf(TabDat.ftr{i,1},TabDat.ftr{i,2}));
    tpl = tpl.parse('resftrs','resftr',1);
end
tpl = tpl.parse('OUT','TPL_RES');
fid = fopen(fHTML,'wt');
fprintf(fid,'%c',get(tpl,'OUT'));
fclose(fid);

%-Display webpage
%==========================================================================
web(fHTML);
