function spm_adjmean_ui
% Scaled (for grand mean) & adjusted (for global) means via General Linear Model
% FORMAT spm_adjmean_ui
%_______________________________________________________________________
%
% spm_adjmean_ui uses the General Linear Model to produce mean images
% adjusted for global effects.
%
% This program is designed for collapsing data within condition to give
% a single adjusted mean scan per condition per subject, suitable for a
% (2nd level) random effects analysis.
%
% See spm_RandFX.man for further details on implementing random effects
% analyses in SPM96 using a multi-level approach.
%
% Overview
% ----------------------------------------------------------------------
% The program supports multiple conditions, multiple subjects, Grand
% Mean (GM) scaling by subject or overall grand mean, proportional
% scaling global normalisation; and AnCova (regression) global
% normalisation, both overall and subject specific, with adjustment to
% subject or overall grand means (after scaling). The availability of
% these options is customised for various designs.
%
% Adjustment is performed via the General Linear Model. The interface
% is similar to SPM-PET, and the adjusted means are the parameter
% estimates from the model. Having chosen a design, the user is
% prompted for scans (by subject and/or condition where appropriate),
% and then asked to set scaling/normalisation/adjustment options as
% appropriate. With the design specified, the model is constructed, and
% the user prompted to enter/confirm filenames for the adjusted mean
% images, which are written to the current working directory (pwd).
%
% The model, filenames, global values and options are saved to a MatLab
% *.mat file named SPMadj.mat in the current working directory.
%
% No masking is carried out: The adjusted mean is calculated at all
% voxels. Voxels which are zero in *all* the input images pertaining to
% an adjusted mean (usually those from the appropriate subject), will
% remain zero, since the computation is only a weighted mean. Modelling
% (& computation) proceeds assumming that at each voxel the data is
% either all zero, or that there is usable data from all images. (This
% is *not* a softmean.) Data realigned in a single session with SPM'96
% (or later) are automatically consistently masked in this way.
%
% GM, the value for grand mean scaling, is user specified.
% The default value is 50.
%
% If computing adjusted means for subsequent (2nd level) modelling, as
% with a random effects analysis, then it is important to use a
% seperable model, such that the adjustment for one subject is
% independent of other subjects entered into the model. Thus,
% proportional scaling or subject-specific AnCova adjustment must be
% used. Further, multiple runs *must* use the same GM value, and should
% scale Grand mean *by subject*.
%
% ( A separate program (spm_adjmean_fmri_ui) is available for computing  )
% ( adjusted condition means of fMRI data, adjusting for global effects  )
% ( and removing low frequency drifts with a high-pass filer (discrete   )
% ( cosine basis set). The functionality is similar to this code, but    )
% ( the two routines have been separated for algorithmic clarity.        )
%
% Diagnostic output
% ----------------------------------------------------------------------
% Diagnostic output consists of two sections:
%
% The first is a list of the image filenames; their global values; and
% the respective subject (block), condition and replication indices; in
% the order in which the the data are entered into the model. This is
% followed by a brief description of appropriate parameters (Grand mean
% scaling etc.)
%
% The second part is a single page depicting the design matrix, effect
% names, parameter contrasts used, and the corresponding image files
% written.
%
% As always, look at the resulting mean images to make sure they look OK!
%
%
% AdjMean "recipies"...
% ----------------------------------------------------------------------
% Rather than offer a bewildering array of options, various
% pre-configured recipies are offered for common scenarios:
%
% * Basic means: +/- grand mean scaling; +/- global normalisation
%
%  1) Straight mean
%       - as the neame suggests! Prompts for files and writes their mean.
%  2) PropSca & Average
%       - Average of images adjusted for global differences by proportional
%         scaling: Scales all images to have global mean of GM, and then
%         writes their mean.
%  3) Linear (AnCova) adjusted mean (scaled mean)
%       - Data scaled so that grand mean is GM. Single mean of images
%         adjusted for global effects by linear regression to mean global.
%         (Actually, this turns out to be a straight mean of the grand mean
%         scaled data, hence the description "scaled mean".)
%  4) Multi linear adjusted means (scaled means)
%       - Multiple block means. Data scaled within blocks to (block) Grand
%         Means of GM. Linear global adjustment within block to block grand
%         mean, and computation of adjusted block means. It's like running
%         option (3) multiple times. Since this is equivalent to grand mean
%         scaling within block and then writing out the block means, it's
%         also tagged "scaled means".
%
% * The "condition" recipies: Adjusted condition means, computed within subj.
%
%  5) SingleSubj: Condition means (PropSca)
%       - Proportional scaling global normalisation of image global means
%         to GM. Computation of means of adjusted data within condition.
%  6) SingleSubj: Condition means (AnCova)
%       - Grand mean scaling of mean global to GM. AnCova global
%         normalisation (parallel lines for each condition).
%         Condition means of AnCova adjusted data written. These are the
%         condition effects of the standard SPM single subject activation
%         AnCova.
%  7) MultiSubj: Condition means (PropSca)
%       - Multiple subject version of option (5).
%         It's like running option (5) repeatedly.
%  8) MultiSubj: Condition means (AnCova by subject)
%       - Multiple subject version of option (6):
%         Grand mean scaling by subject, AnCova by subject.
%         It's like running option (6) repeatedly.
%
%
% Algorithm
% ----------------------------------------------------------------------
% The model at each voxel is Y = X*B + e, with least squares estimates
% (for full rank X) for the vector B of parameters as b =
% inv(X'*X)*X'*Y. For c a vector of contrast weights extracting the
% appropriate parameter, the contrast of the parameter estimates is
% c'*b = c'* inv(X'*X)*X' * Y, a weighted sum (or weighted mean) of the
% data at that voxel. These weights are identical for all voxels, so
% the image of the parameter estimate can be computed as a weighted
% mean of the images.
%
% Once the weights have been worked out for each adjusted mean image,
% computation proceeds by passing appropriate weights and image
% filenames to spm_mean, which writes out the appropriate parameter
% image as an Analyze format image of the same type (see spm_type) as
% the input images.
%
% Variables saved in SPMadj.mat data file
% ----------------------------------------------------------------------
% iDesType	Design type index
% DesName	Design name
% iHForm	Form of DesMtx H partition
% iGMsca	Grand mean scaling option
% sGMsca	Grand mean scaling description
% iGloNorm	Global normalisation option
% sGloNorm	Global normalisation description
% iAdjTo	Adjustment (for AnCova) option
% sAdjTo	Adjustment (for AnCova) description
% iSubj		Subject indicator vector
% nSubj		number of subjects
% iCond		Condition indicator vector
% iRepl		Replication indicator vector
% P		String matrix of filenames
% H		H (condition effects) partition of design matrix
% Hnames	Names of parameters corresponding to columns of H
% G		G (confounding covariates) partition of design matrix
% Gnames	Names of parameters corresponding to columns of G
% c		Matrix of contrasts, contrasts in rows
% cNames	Names associated with contrasts
% W		Weights for images corresponding to contrasts
% CWD		Current Working Directory (when run)
% Fnames	Filenames of adjusted mean images written
% rGX		raw global means (before any scaling)
% gSF		Image scale factors for global scaling (inc. grand mean scaling)
% GX		Global means after scaling
% GM		Grans Mean used for scaling
%
%
% Platform
% ----------------------------------------------------------------------
% This version was written for MatLab4.2c with SPM'97d (spm_mean.m v1.11)
%_______________________________________________________________________
% %W% Andrew Holmes %E%

%SPM96:This version was written for MatLab4.2c with SPM'96 (spm_mean.m v1.10)

%=======================================================================
% - S E T U P
%=======================================================================
fprintf('\nSPM: spm_adjmean_ui\n'),fprintf('%c','='*ones(1,72)),fprintf('\n')
spm_clf('Interactive')
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), Finter=spm('CreateIntWin'); end
set(Finter,'Name','AdjMean')
spm_help('!ContextHelp','spm_adjmean_ui.m')



%=======================================================================
% - D E S I G N   P A R A M E T E R S
%=======================================================================
DesNames = str2mat(...
	'Straight mean',...					%-1
	'PropSca & Average',...					%-2
	'Linear (AnCova) adjusted mean (scaled mean)',...	%-3
	'Multi linear adjusted means (scaled means)',...	%-4
	'SingleSubj: Condition means (PropSca)',...		%-5
	'SingleSubj: Condition means (AnCova)',...		%-6
	'MultiSubj: Condition means (PropSca)',...		%-7
	'MultiSubj: Condition means (AnCova by subject)');	%-8
iDesType = 2;	%-Default design type

DesPrams = str2mat(...
'bMSubj',...
	'bBlock',...
		'bMCond',...
			'iHForm',...
				'iGloNorm',...
					'iGMsca',...
						'iAdjTo');

DesDefaults = [ ...
1,	1,	0,	2,	1,	1,	4;...		%-1
1,	1,	0,	1,	2,	4,	4;...		%-2
0,	1,	0,	1,	3,	2,	2;...		%-3
1,	1,	0,	1,	4,	3,	3;...		%-4
0,	0,	1,	3,	2,	4,	4;...		%-5
0,	0,	1,	3,	3,	2,	2;...		%-6
1,	0,	1,	4,	2,	4,	4;...		%-7
1,	0,	1,	4,	4,	3,	3;	];	%-8

HForms = str2mat(...
	'iSubj,''-'',''AdjMean''',...
	'iSubj,''-'',''mean''',...
	'iCond,''-'',''Cond''',...
	'[iSubj,iCond],''-'',str2mat(''Subj'',''Cond'')');

%-Global normalization options
sGloNorm = str2mat(...
	'No Global Normalisation',...
	'Proportional scaling',...
	'AnCova',...
	'AnCova {subject-specific}');

%-Grand mean scaling options
sGMsca = str2mat(...
	'No Grand Mean Scaling',...
	'Scaling of overall Grand Mean',...
	'Scaling of subject Grand Means',...
	'(Implicit in PropSca global normalisation)');
%-NB: Grand mean scaling by subject is redundent for proportional scaling

%-Adjustment options for AnCova designs (for centering of globals)
%-If Grand mean scaling, then would usually AnCova adjust in a similar
% fashion, i.e. to GM.
sAdjTo = str2mat(...
	'Specify...',...
	'Grand mean (mean of all globals)',...
	'Subject grand mean (mean of subjects globals)',...
	'(redundant: not doing AnCova)');

if ( size(DesNames,1) ~= size(DesDefaults,1)) | ...
	(size(DesPrams,1) ~= size(DesDefaults,2))
	error('Size mismatch in design parameter specification')
end

%-Initialise indicies
%-----------------------------------------------------------------------
iSubj   = [];		% Subject (block) index
iCond   = [];		% condition (or scan) (per subject) index
iRepl   = [];		% replication (per condition) index
P       = [];		% string matrix of filenames


%-Get parameters
%=======================================================================

%-Select design
%-----------------------------------------------------------------------
iDesType = spm_input('Select mean type...',1,'m',DesNames,[],iDesType);
DesName = deblank(DesNames(iDesType,:));
for p   = 1:size(DesPrams,1)
    eval([deblank(DesPrams(p,:)),' = DesDefaults(iDesType,p);']), end
HForm   = HForms(iHForm,:);
if bBlock, sSubBlk='block'; else, sSubBlk='subject'; end


%-Get filenames, accrue subject, condition & replication indicies
%-----------------------------------------------------------------------
guiPos = '+1'; %*** guiPos = spm_input('!NextPos');
nSubj     = 1;
if bMSubj
	nSubj = spm_input(['number of ',sSubBlk,'s ?'],guiPos);
	bMSubj = nSubj > 1;
end
for subj  = 1:nSubj
	sSubj = ''; if bMSubj, sSubj = [sSubBlk,' ',int2str(subj),': ']; end
	if bMCond, nCond = spm_input([sSubj,'# of conditions ? '],'0');
		else, nCond = 1; end
	for cond  = 1:nCond
	    if nCond > 1, t_str = sprintf('cond %d: ',cond);
	    	else t_str=''; end
	    tP    = spm_get(Inf,'.img',[sSubj,t_str,'Select scans...']);
	    nRepl = size(tP,1);
	    P     = str2mat(P,tP);
	    iSubj = [iSubj; subj*ones(nRepl,1)];
	    iCond = [iCond; cond*ones(nRepl,1)];
	    iRepl = [iRepl; [1:nRepl]'];
	end
end
P(1,:)  = [];

%-Total #observations
%-----------------------------------------------------------------------
q       = length(iSubj);

%-H partition
%-----------------------------------------------------------------------
eval(['[H,Hnames] = spm_DesMtx(',HForm,');'])

%-Global normalization options
%-----------------------------------------------------------------------
if iGloNorm>9
	%-User has a choice from the options in iGloNorm.
	%-iGloNorm contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGloNorm);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Don't offer subject specific AnCova if not bMSubj
	if ~bMSubj, tmp(find(tmp==4))=[]; end
	iGloNorm=spm_input('Select global normalisation','+1','m',...
	    	sGloNorm(tmp,:),tmp);
end
sGloNorm = deblank(sGloNorm(iGloNorm,:));

%-Grand mean scaling options
%-----------------------------------------------------------------------
%-Grand mean scaling is implicit in PropSca global normalisation
if iGloNorm==2, iGMsca=4; end
if iGMsca>9
	%-User has a choice from the options in iGMsca.
	%-iGMsca contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iGMsca);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Scaling by subject redundent if proportional scaling,
	%-Don't offer subject specifics if not bMSubj
	if ~bMSubj, tmp(find(tmp==3))=[]; end
	iGMsca=spm_input...
	    ('Grand mean scaling','+1','m',sGMsca(tmp,:),tmp);
end
if iGMsca>1
	%-Get value for grand mean scaling
	if iGloNorm==2, str='GM: PropSca global mean to ?';
		else, str='GM: Scale grand mean to ?'; end
	GM = spm_input(str,'+1','e',50);
	if GM==0, iGMsca=0; end
else
	GM=0;
end
sGMsca = deblank(sGMsca(iGMsca,:));


%-Adjustment options for AnCova designs (for centering of globals)
%-----------------------------------------------------------------------
if any(iGloNorm==[1,2]), iAdjTo=4; end
if iAdjTo>9
	%-User has a choice from the options in iAdjTo.
	%-iAdjTo contains an integer, each digit specifies an option
	%---------------------------------------------------------------
	str = int2str(iAdjTo);
	tmp = []; for i = 1:length(str), tmp = [tmp, eval(str(i))]; end
	%-Don't offer subject specifics if not bMSubj
	if ~bMSubj, tmp(find(tmp==2))=[]; end
	iAdjTo=spm_input...
		('AnCova adjust (centre globals), after any scaling to',...
		 '+1','m',sAdjTo(tmp,:),tmp);
end
if iAdjTo==1, aGM = spm_input('AnCova adjust to ?','+1','e',GM); end
sAdjTo = deblank(sAdjTo(iAdjTo,:));



%=======================================================================
% - C O M P U T A T I O N
%=======================================================================
set(Finter,'Name','AdjMean - computing','Pointer','watch')
fprintf('\tcomputing: ')

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
fprintf('(globals)')
GX     = zeros(q,1);
for i  = 1:q, GX(i) = spm_global(V(:,i)); end
fprintf('\b - done)\n')

%-Save scalefactors
%-----------------------------------------------------------------------
iSF = V(7,:)';


%-Scaling: compute global scaling factors required to implement proportional
% scaling global normalisation or Grand mean scaling, as requested
%-----------------------------------------------------------------------
rGX = GX;
if iGloNorm==2
	%-Proportional scaling global normalisation
	gSF = GM./GX;
	GX  = GM*ones(q,1);
	%** scale rGX for printout? ...or just graph them?
elseif iGMsca==2
	%-Grand mean scaling (overall)
	gSF = GM/mean(GX);
	GX  = GX*gSF;
	sGMsca = sprintf('%s to %g',sGMsca,GM);
elseif iGMsca==3
	%-Grand mean scaling by subject
	gSF = GM./spm_meanby(GX,iSubj);
	GX  = GX.*gSF;
	sGMsca = sprintf('%s to %g',sGMsca,GM);
else	%-No scaling
	gSF = ones(q,1);
end

%-AnCova options: Construct Global covariates of no interest partition
%-----------------------------------------------------------------------
if any(iGloNorm==[3,4])
	if iAdjTo==1
		%-aGM set by user
	elseif iAdjTo==2
		aGM = mean(GX);
	elseif iAdjTo==3
		aGM = spm_meanby(GX,iSubj);
	else,	error('Illegal iAdjTo')
	end

	if iGloNorm == 3	%-AnCova
	    G  = GX - aGM; Gnames = 'Global';
	elseif iGloNorm == 4	%-AnCova by subject
	    [G,Gnames] = ...
		spm_DesMtx([iSubj,GX-aGM],'FxC',str2mat(sSubBlk,'global'));
	end
else
	G = []; Gnames = '';
end


%-Design matrix, parameter estimation matrix, "hat" matrix, "residual" matrix
%-----------------------------------------------------------------------
X        = [H G];
[nX,Xnames] = spm_DesMtxSca(H,Hnames,G,Gnames);
XTXinvX  = inv(X'*X)*X';

%-Contrasts (c)
%-----------------------------------------------------------------------
c  = [eye(size(H,2)), zeros(size(H,2),size(G,2))];
nc = size(c,1);
W  = c * XTXinvX;

cNames = Hnames;

%-Parameter images (of interest) - Adjusted mean images
%-----------------------------------------------------------------------
set(Finter,'Name','AdjMean - writing','Pointer','arrow')
fprintf('\twriting parameter images... enter filenames...\n')

%-Computation - calculations handled by spm_mean.c
% Relying on all images used to make each contrast image being masked
% identically, with out of mask values set to zero. This is the case
% for sets of images realigned with SPM'96 realignment.

Fnames = []; CWD = pwd;
guiPos = '+1'; %*** guiPos = spm_input('!NextPos');
for i = 1:nc
	Fn = deblank(cNames(i,:)); Fn=Fn(Fn~=' ');
	Fn = spm_input(sprintf('Contrast %d: filename?',i),guiPos,'s',Fn);
	%*** make sure Fn is valid? Allow overwriting?
	Fnames = str2mat(Fnames,Fn);
	w  = W(i,:).*gSF'.*iSF';
	Q  = find(abs(w)>0);
	wV = V(:,Q);	wV(7,:) = w(Q);		%SPM97:
	fprintf('\t...writing image %d: %-20s',i,Fn)
	sf = spm_mean(wV,[Fn,'.img']);		%SPM97:
	%SPM96:sf = spm_mean(prod(DIM),TYPE,[Fn,'.img'],P(Q,:),w(Q));
	str = sprintf('Adjusted mean (spm_adjmean) - %s',Fn);
	spm_hwrite([Fn,'.hdr'],DIM,VOX,sf,spm_type('int16'),OFFSET,ORIGIN,str);	%SPM97:
	%SPM96:spm_hwrite([Fn,'.hdr'],DIM,VOX,sf,TYPE,OFFSET,ORIGIN,str)
	spm_get_space(Fn,spm_get_space(P(1,:)));
	fprintf(' (done)\n')
	guiPos = '0';
end
Fnames(1,:)=[];

%-Unmap files and canonicalise V
%-----------------------------------------------------------------------
for v = V; spm_unmap(v); end
V = [V(1:6,1); ORIGIN(:)];

%-Save parameters to SPMadj.mat in current directory
%-----------------------------------------------------------------------
save SPMadj ...
	iDesType DesName ...
	iHForm iGloNorm iAdjTo iGMsca ...
	sGMsca sGloNorm sAdjTo ...
	iSubj nSubj iCond iRepl P ...
	H Hnames G Gnames c cNames W ...
	CWD Fnames ...
	rGX gSF GX GM



%=======================================================================
% - D I A G N O S T I C   O U T P U T
%=======================================================================
set(Finter,'Name','AdjMean - done')
fprintf('\tdisplaying diagnostic output...\n')

Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph = spm_figure('CreateWin','Graphics','Graphics'); end

% plot(iSubj,GX,'x'), xlabel('Subject'), ylabel('Global')


%-Display files and variables
%=======================================================================

%-Compute common path components - all paths will begin with '/'
%-----------------------------------------------------------------------
d     = max(find(P(1,1:min(find(~all(P == ones(q,1)*P(1,:))))-1)=='/')) - 1;
CPath = P(1,1:d);
Q     = P(:,[(d + 2):size(P,2)]);

%-Display
%-----------------------------------------------------------------------
figure(Fgraph); spm_clf; axis off
text(0.30,1.02,'Adjusted means','Fontsize',16,'Fontweight','Bold')
text(-0.10,0.85,'Scan Index','Rotation',90)
if bMSubj, text(+0.00,0.85, sSubBlk,     'Rotation',90), end
if bMCond, text(+0.05,0.85,'condition',  'Rotation',90), end
text(+0.10,0.85,'Replication','rotation',90)
x0    = 0.15; y0 = 0.83;
dx    = 0.10; dy = 0.02;
x     = x0;
text(x + 0.02,0.85,'global','Rotation',90), x = x + dx;

text(x,0.92,'Base directory:','FontSize',10,'Fontweight','Bold')
text(x,0.90,CPath,'FontSize',10)
text(x,0.87,'Filename Tails')
y     = y0;
for i = 1:q
	text(-0.12,y,sprintf('%02d :',i));
	if bMSubj, text(-0.01,y,sprintf('%2d',iSubj(i))), end
	if bMCond, text(+0.04,y,sprintf('%2d',iCond(i))), end
	text(+0.09,y,sprintf('%2d',iRepl(i)))
	x     = x0;
	text(x,y,sprintf('%-8.6g',GX(i)),'FontSize',10), x = x + 2*dx;
	text(x,y,Q(i,:),'FontSize',10)
	y     = y - dy;
	if y < 0;
		spm_print
		spm_clf(Fgraph); axis off
		y = y0;
		text(0.16,1.02,['Adjusted means (continued)'],...
		    'Fontsize',16,'Fontweight','Bold')
	end
end

y      = y - dy;
dy     = dy*1.2;
text(0,y,['Grand mean scaling: ',sGMsca]), y = y - dy;
text(0,y,['Global normalisation: ',sGloNorm]), y = y - dy;
text(0,y,['AnCova adjustment to: ',sAdjTo]), y = y - dy;
text(0,y,sprintf('Parameters saved to: %s/SPMadj.mat',pwd),'FontSize',10)

spm_print


%-Depict and label design matrix, depict & label contrasts
%=======================================================================
spm_clf(Fgraph); axis off
text(0.1,1.02,'Design Matrix, contrasts & contrast files',...
	'Fontsize',16,'Fontweight','Bold');

%-Image scaled design matrix & label the effects
%-----------------------------------------------------------------------
hDesMtx = axes('Position',[0.1 0.5 0.6 0.3]);
image((nX + 1)*32);
ylabel('scans')
xlabel('parameters')
hEfLabs = axes('Position',[0.1 0.8 0.6 0.1],'Visible','off');
y     = 0.1;
dx    = 1/size(nX,2);
for i = 1:size(nX,2)
	text((i - 0.5)*dx,y,deblank(Xnames(i,:)),'Fontsize',8,'Rotation',90)
end

%-Depict contrasts and associated filenames
%-----------------------------------------------------------------------
dy = 0.4/nc;
axes('Position',[0.025 0.05 0.05 0.4],'Visible','off')
text(0,0.5,'contrasts','HorizontalAlignment','Center','Rotation',90,...
	'FontSize',14,'FontWeight','Bold')
axes('Position',[0.6 0.44 0.40 0.02],'Visible','off')
text(0,1,'Contrast files...','FontSize',10,'FontWeight','Bold')
text(0,0,sprintf('...in %s',CWD),'FontSize',8)
for i = 1:nc
	axes('Position',[0.1 (0.45 -dy*i) 0.6 0.9*dy])
	[tx ty] = bar(c(i,:));
	fill(tx,ty,[1 1 1]*.8);
	% h = bar(c(i,:),1);
	% set(h,'FaceColor',[1 1 1]*.8)
	set(gca,'XLim',[0.5,size(c,2)+0.5],'Visible','off')
	text(0,0,num2str(i),'HorizontalAlignment','Right','FontSize',10)
	text(size(c,2)+.55,0,sprintf('%s',deblank(Fnames(i,:))),'FontSize',10)
end

spm_print


%-END
%=======================================================================
fprintf('\n\n')
