function  spm_projectionsF_ui()
% used to review results of statistical analysis (SPM{F})
% FORMAT spm_projectionsF_ui
%_______________________________________________________________________
%
% spm_projectionsF_ui allows the SPM{F} created by spm_spm.m to be re-displayed
% and characterized in terms of regionally significant effects. This
% is based on K. Worsley results for the expected maximum value in
% an F-field and requires the smoothness estimation of the original   
% component fields.
%
%
% see spm_projectionsF.m for further details
%
%_______________________________________________________________________
% %W% JBP %E%


%-Get SPMt.mat for analysis, and load results
%-----------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
spm_clf(Finter)
set(Finter,'Name','SPM{F} projections')

tmp = spm_get(1,'.mat','select SPMF.mat for analysis','SPMF');
global CWD
CWD = strrep(tmp,'/SPMF.mat','');
K   = [];

load([CWD,'/SPM'])
load([CWD,'/XYZ'])
load([CWD,'/SPMF'])

if isempty([B G]) | isempty([H C])
	df = [rank([H C B G]), df];
	if ~isempty(H), df=df-[1,0]; end %-See spm_spm!
else
	df   = [rank([H C]) df];
end

U    = spm_input(' F threshold ? ',2,'e',3.2);
if U < 1 U = spm_invFcdf(1 - U,[df]); end

%-Pass arguments to spm_projections
%-----------------------------------------------------------------------
set(Finter,'Name','Thankyou','Pointer','watch')

spm_projectionsF(SPMF,XYZ,U,V,Wresid,S,[K H C B G],df,size([K H C],2));

spm_clf(Finter); return;
