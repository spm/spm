function spm_adjmean_ui
% 
% FORMAT [P,GX] = spm_adjmean_ui
%
% No global scaling for proportional scaling scales to block grand mean.
% spm_spm_ui effectively scales to 1 with no Grand Mean Scaling
%
%_______________________________________________________________________
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%

fprintf('\nSPM: spm_adjmean_ui\n')

%-Design parameters
%=======================================================================
Designs = str2mat(...
	'PropSca & Average',...
	'Single-block : General',...
	'Multi-block  : General');
DefDesType = 1;

DesPrams = str2mat(...
'bMBlok',...
	'iHForm',...
		'iGloNorm',...
			'iAdjTo',...
				'iGMsca',...
					'iWrite');

DesDefaults = [ ...
0,	1,	2,	1,	1,	1;...
0,	1,	123,	1,	12,	1;...
1,	2,	1234,	12,	123,	1	];

HForms = str2mat(...
	'iBlok,''-'',''AdjMean''',...
	'iBlok,''-'',''Block''',...
	'iCond,''-'',''Cond''',...
	'[iBlok'',iCond''],''-'',str2mat(''Blok'',''Cond'')');

sGloNorm = str2mat(...
	'No Global Normalisation',...
	'Proportional scaling',...
	'AnCova',...
	'AnCova {block-specific}');

sAdjTo = str2mat(...
	'Grand mean (mean of all globals)',...
	'Block grand mean (mean of block globals)');

sGMsca = str2mat(...
	'No Grand Mean Scaling',...
	'Scaling of overall Grand Mean',...
	'Scaling of block Grand Means');
%-NB: Grand mean scaling by block is redundent for proportional scaling

sWrite = str2mat(...
	'mean only',...
	'adjusted data only',...
	'mean & adjusted data');

if ( size(Designs,1) ~= size(DesDefaults,1)) | ...
	(size(DesPrams,1) ~= size(DesDefaults,2))
	fprintf('%cSize mismatch in design parameter specification\n',7)
	return
end

%-Initialise indicies
%-----------------------------------------------------------------------
bMCond  = 0;		% Multi-Cond not supported (yet)
iBlok   = [];		% block index
iCond   = [];		% condition (or scan) (per subject) index
iRepl   = [];		% replication (per condition) index
P       = [];		% string matrix of filenames

%-Get parameters
%=======================================================================

%-Initialise Figure windows
%-----------------------------------------------------------------------
spm_clf('Interactive')
% set(Finter,'Name','AdjMean')

%-Select design
%-----------------------------------------------------------------------
DesType = spm_input('Select mean type...','+1','m',Designs,[],DefDesType);
DesName = deblank(Designs(DesType,:));
for p   = 1:size(DesPrams,1)
    eval([deblank(DesPrams(p,:)),' = DesDefaults(DesType,p);']), end
HForm   = HForms(iHForm,:);


%-Get filenames, accrue block, condition & replication indicies
%-----------------------------------------------------------------------
nBlok     = 1;
if bMBlok
	nBlok = spm_input('Number of blocks ?','+1');
	bMBlok = nBlok > 1;
end
for blok  = 1:nBlok
	sBlok = []; if bMBlok, sBlok = ['Block ',int2str(blok),': ']; end
	if bMCond, nCond = spm_input([sBlok,'# of conditions ? '],'0');
		else, nCond = 1; end
	for cond  = 1:nCond
	    if nCond > 1, t_str = sprintf('Cond %d: ',cond);
	    	else t_str=''; end
	    tP    = spm_get(Inf,'.img',[sBlok,t_str,'Select scans...']);
	    nRepl = size(tP,1);
	    P     = str2mat(P,tP);
	    iBlok = [iBlok; blok*ones(nRepl,1)];
	    iCond = [iCond; cond*ones(nRepl,1)];
	    iRepl = [iRepl; [1:nRepl]'];
	end
end
P(1,:)  = [];

%-Total #observations
%-----------------------------------------------------------------------
q       = length(iBlok);

%-H partition
%-----------------------------------------------------------------------
eval(['[H,Hnames] = spm_DesMtx(',HForm,');'])

%-Create iCOND indicators, indicating conditions
% uniquely across blocks. Watch out for unbalanced designs!
%-Construct nCOND as total number of conditions across blocks.
%-----------------------------------------------------------------------
nCond   = []; for blok=1:nBlok, nCond=[nCond;max(iCond(iBlok==blok))]; end
nCOND	= sum(nCond);               			%-#conditions in total
tmp	= cumsum([0;nCond]);
iCOND   = iCond+tmp(cumsum([1;diff(iBlok)])); 		%-Index to conditions


%=======================================================================

%-Global normalization options
%-----------------------------------------------------------------------
if iGloNorm>9
	%-User has a choice from the options in iGloNorm.
	%-iGloNorm contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGloNorm);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Don't offer block specific AnCova if not bMBlok
	if ~bMBlok, tmp(find(tmp==4))=[]; end
	iGloNorm=spm_input('Select global normalisation','+1','m',...
	    	sGloNorm(tmp,:),tmp);
end

%-Adjustment options for AnCova designs
%-----------------------------------------------------------------------
%-If doing proportional global normalisation, or none at all, then
% adjustment is effectively to the grand mean (after scaling).
if any(iGloNorm==[1,2]), iAdjTo=1; end
if iAdjTo>9
	%-User has a choice from the options in iAdjTo.
	%-iAdjTo contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iAdjTo);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Don't offer block specifics if not bMBlok
	if ~bMBlok, tmp(find(tmp==3))=[]; end
	iAdjTo=spm_input('Adjust to','+1','m',...
	    	sAdjTo(tmp,:),tmp);
end

%-Grand mean scaling options
%-----------------------------------------------------------------------
if iGMsca>9
	%-User has a choice from the options in iGMsca.
	%-iGMsca contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGMsca);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Scaling by block redundent if proportional scaling,
	% don't offer block specifics if not bMBlok
	if (iGloNorm==2 | ~bMBlok) & any(tmp==3), tmp(find(tmp==3))=[]; end
	iGMsca=spm_input...
	    ('Grand mean scaling','+1','m',sGMsca(tmp,:),tmp);
end
if iGMsca>1,	GM = spm_input('Value for grand mean ?','+1','e',50);
		if GM==0, iGMsca=1; end
else, GM=0; end


%-Image writing options
%-----------------------------------------------------------------------
%if iWrite>9
%	%-User has a choice from the options in iWrite.
%	%-iWrite contains an integer, each digit specifies an option
%	%---------------------------------------------------------------
%	str = int2str(iWrite);
%	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
%	iWrite=spm_input('Images to write','+1','m',...
%	    	sWrite(tmp,:),tmp);
%end


%-Computation - Design matrix
%=======================================================================

%-Get file identifiers
%-----------------------------------------------------------------------
V     = zeros(12,q);
for i = 1:q; V(:,i) = spm_map(P(i,:));  end

%-Check for consistency of image size and voxel size
%-----------------------------------------------------------------------
if ~all(all(~diff(V([1:6],:)')))
	error('data do not have the same image and voxel size'); end

%-Get ORIGIN from first image
%-----------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(P(1,:));

%-Compute global values
%-----------------------------------------------------------------------
GX     = zeros(q,1);
for i  = 1:q
	GX(i) = spm_global(V(:,i)); end

%-Save scalefactors, unmap files and canonicalise V
%-----------------------------------------------------------------------
SF = V(7,:)';
for v = V; spm_unmap(v); end
V = [V(1:6,1); ORIGIN(:)];


%-Grand mean scaling for AnCova / no normalisation options (if required)
%-Scale scaling coefficients so that Grand mean (mean of global means) is GM
%-----------------------------------------------------------------------
rGX = GX; rSF = SF;
if iGMsca==1
	%-No grand mean scaling
elseif iGMsca==2
	%-Grand mean scaling (overall)
	SF = SF*GM/mean(rGX);
	GX = GX*GM/mean(rGX);
elseif iGMsca==3
	%-Grand mean scaling by block
	SF = SF.*GM./spm_meanby(rGX,iBlok);
	GX = GX.*GM./spm_meanby(rGX,iBlok);
else, error('Invalid iGMsca option'), end
sGX = GX;

%-Construct Global covariates of no interest partition.
%-Centre global means if included in AnCova models, by block if requested
%-----------------------------------------------------------------------
if iGloNorm == 1
	%-No global adjustment -----------------------------------------
	G = []; GNames = '';
elseif iGloNorm == 2
	%-Proportional scaling -----------------------------------------
	G = []; GNames = '';
	if (GM ~= 0) 	SF = rSF./rGX'*GM; GX = ones(size(rGX))*GM;
	else
		SF = rSF.*spm_meanby(GX,iBlok)./rGX;
		GX = spm_meanby(GX,iBlok);
	end
elseif iGloNorm == 3
	%-AnCova -------------------------------------------------------
	if iAdjTo==1
		G = GX - mean(GX);
		Gnames = 'cGlobal';
	elseif iAdjTo==2
		G = GX - spm_meanby(GX,iBlok);
		Gnames = 'bcGlobal';
	else, error('Illegal iAdjTo'), end
	Gnames = 'Global';
elseif iGloNorm == 5
	%-AnCova by block ----------------------------------------------
	if iAdjTo==1
		[G,Gnames] = spm_DesMtx(...
			[iBlok',GX-mean(GX)],...
				'FxC',['Block ';'cGlobal']);
	elseif iAdjTo==2
		[G,Gnames] = spm_DesMtx(...
			[iBlok',GX-spm_meanby(GX,iBlok)],...
				'FxC',['Block ';'bcGlobal']);
	else, error('Illegal iAdjTo'), end
else, error('Invalid iGloNorm option\n',7), end


%-Computation - calculations handled by spm_mean.c
%=======================================================================
X        = [H G];
XTXinvX  = inv(X'*X)*X';
%Hat      = X*XTXinvX;
%R        = eye(size(Hat))-Hat;

%-Save parameters to SPMadj.mat in current directory
%-----------------------------------------------------------------------
%save SPMadj ...
%	DesType DesName ...
%	iHForm iGloNorm iAdjTo iGMsca ...
%	iBlok nBlok iCond nCond iCOND nCOND iRepl P ...
%	H Hnames G Gnames rGX sGX GX
	

%-Parameter images (of interest) - Adjusted mean images
%-----------------------------------------------------------------------
fprintf('\tWriting parameter images...\n')
guiPos = '+1'; %**** guiPos = spm_input('!NextPos');
for i = 1:size(H,2)
	Q = deblank(Hnames(i,:));
	Q = spm_input(sprintf('%s: file to write?',Q),guiPos,'s',lower(Q));
	sf = spm_mean(prod(DIM),TYPE,[Q,'.img'],P,XTXinvX(i,:).*SF');
	DESCRIP = sprintf('Adjusted mean (spm_adjmean) - %s',Q);
	spm_hwrite([Q,'.hdr'],DIM,VOX,sf,TYPE,OFFSET,ORIGIN,DESCRIP);
	spm_get_space(Q,spm_get_space(P(1,:)));
	fprintf('\t...written parameter image %d: %s\n',i,Q)
end


%-Adjusted images
%-----------------------------------------------------------------------
%G_XTXinvX = XTXinvX(size())
%for i = 1:q
%	Q = [deblank(Hnames(i,:)),'.img'];
%	sf = spm_mean(prod(DIM),TYPE,Q,P,XTXinvX(i,:).*SF');
%	DESCRIP = sprintf('Adjusted mean (spm_adjmean) - %s',Q);
%	spm_hwrite(Q,DIM,VOX,sf,TYPE,OFFSET,ORIGIN,DESCRIP);
%	spm_get_space(Q,spm_get_space(P(1,:)));
%end



%-Diagnostic output
%=======================================================================
% plot(iSubj,GX,'x'), xlabel('Subject'), ylabel('Global')



%-END
%=======================================================================
fprintf('\n\n')
