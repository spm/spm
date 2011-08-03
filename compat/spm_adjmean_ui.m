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
% Implicit masking is carried out: Zero voxels are implicitly assummed
% to be masked out. Thus, the adjusted mean is calculated at voxels
% which are non-zero in *all* the input images pertaining to the
% adjusted mean (usually those from the appropriate subject). (This is
% *not* a softmean.) Data realigned in a single session with SPM'96 (or
% later) are automatically implicitly zero masked with a consistent
% mask in this way.
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
% filenames to spm_add, which writes out the appropriate parameter
% image as an Analyze format image of the same type (see spm_type) as
% the input images.
%
% Variables saved in SPMadj.mat data file
% ----------------------------------------------------------------------
% DesDef        Structure containing defaults for chosen design
% Des           Structure containing design parameters & specification
%   .DesName    Design name
%   .HForm      Form of DesMtx H partition
%   .iSubj      Subject indicator vector
%   .iCond      Condition indicator vector
%   .iRepl      Replication indicator vector
%   .iGloNorm   Global normalisation option
%   .sGloNorm   Global normalisation description
%   .iGMsca     Grand mean scaling option
%   .sGMsca     Grand mean scaling description
%   .GM         Grand Mean used for scaling
%   .iAdjTo     Adjustment (for AnCova) option
%   .sAdjTo     Adjustment (for AnCova) description
%   .aGM        AnCova adjustment value (subtracted from GX before AnCova)
%   .gSF        Image scale factors for global scaling
%   .X          Design matrix
%   .nX         Normalised (for imaging) design matrix
%   .Xnames     Effects corresponding to cols of X (cellstr)
%   .aPMap      Additional parameter to effect name mappings (see spm_desMtx)
%   .EXnames    English effect names corresponding to TeX parameters of Xnames 
%   .iX         Structure defining design matrix subpartitions
%       .H      Columns of X corresponding to H partition
%       .C      Columns of X corresponding to C partition
%       .B      Columns of X corresponding to B partition
%       .G      Columns of X corresponding to G partition
% c             Matrix of contrasts, contrasts in rows
% cNames        Names associated with contrasts
% W             Weights for images corresponding to contrasts
% Fnames        Filenames of adjusted mean images written (cellstr)
% rGX           raw global means (before any scaling)
% GX            Global means after scaling
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_adjmean_ui.m 4418 2011-08-03 12:00:13Z guillaume $



%=======================================================================
% - S E T U P
%=======================================================================
SCCSid = '2.8';
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','AdjMean',1);
spm_help('!ContextHelp',[mfilename,'.m'])


%=======================================================================
% - D E S I G N   P A R A M E T E R S
%=======================================================================

%-Design default definitions
%-----------------------------------------------------------------------
DesDefF = {...
'DesName';                      %-Design name
    'HForm';...                 %-Form of H partition
                    'aPMap';... %-Parameter mappings 
    'bMSubj';...                    %-#subject/block 0=ask
        'bMCond';...                %-#conditions 0=ask
            'iGloNorm';...          %-GloNorm codes
                'iGMsca';...        %-GMsca codes
                    'GM';...    %-Grand Mean value
                        'iAdjTo'};  %-AdjTo codes
DesDef = {...
'Straight mean',...
    'iSubj,''-'',''\mu''',      {'\mu','mean';'\gamma','blok'},...
    1,  0,  1,  1,  [], 4;...
'PropSca & average',...
    'iSubj,''-'',''\mu''',      {'\mu','adjmean';'\gamma','blok'},...
    1,  0,  2,  4,  [], 4;...
'Linear (AnCova) adjusted mean (scaled mean)',...
    'iSubj,''-'',''\mu''',      {'\mu','scamean'},...
    0,  0,  3,  2,  [], 2;...
'Multi linear adjusted means (scaled means)',...
    'iSubj,''-'',''\mu''',      {'\mu','adjmean';'\gamma','blok'},...
    1,  0,  4,  3,  [], 3;...
'SingleSubj: Condition means (PropSca)',...
    'iCond,''-'',''\alpha''',   '',...
    0,  1,  2,  4,  [], 4;...
'SingleSubj: Condition means (AnCova)', ...
    'iCond,''-'',''\alpha''',   '',...
    0,  1,  3,  2,  [], 2;...
'MultiSubj: Condition means (PropSca)',...
    '[iSubj,iCond],''-'',{''\gamma'',''\alpha''}',  '',...
    1,  1,  2,  4,  [], 4;...
'MultiSubj: Condition means (AnCova by subject)',...
    '[iSubj,iCond],''-'',{''\gamma'',''\alpha''}',  '',...
    1,  1,  4,  3,  [], 3   };

iDefDesDef = 2;     %-Default Design definition

%-Options
%-----------------------------------------------------------------------
%-Global normalization options
sGloNorm = {    'No Global Normalisation',...               %-1
        'Proportional scaling',...              %-2
        'AnCova',...                        %-3
        'AnCova {subject-specific}'};               %-4

%-Grand mean scaling options
sGMsca = {  'No Grand Mean Scaling',...             %-1
        'Scaling of overall Grand Mean',...         %-2
        'Scaling of subject Grand Means',...            %-3
        '(Implicit in PropSca global normalisation)'};      %-4
%-NB: Grand mean scaling by subject is redundent for proportional scaling
dGM = 50;       %-Default Grand Mean value

%-Adjustment options for AnCova designs (for centering of globals)
%-If Grand mean scaling, then would usually AnCova adjust in a similar
% fashion, i.e. to GM.
sAdjTo = {  'Specify...',...                    %-1
        'Grand mean (mean of all globals)',...          %-2
        'Subject grand mean (mean of subjects globals)',... %-3
        '(redundant: not doing AnCova)'};           %-4


%=======================================================================
% - G E T   I M A G E S   &   P A R A M E T E R S
%=======================================================================

%-Initialise indicies
%-----------------------------------------------------------------------
iSubj   = [];       % Subject (block) index
iCond   = [];       % condition (or scan) (per subject) index
P       = {};       % cell array of string filenames


%-Select design & unpack design specification defaults
%-----------------------------------------------------------------------
DesDef = cell2struct(DesDef(...
    spm_input('Select mean type...',1,'m',DesDef(:,1),[],iDefDesDef)...
        ,:)',DesDefF);
DesName  = DesDef.DesName;
HForm    = DesDef.HForm;
aPMap    = DesDef.aPMap;
bMSubj   = DesDef.bMSubj;
sSubj    = spm_DesMtx('ETeXNames','\gamma',aPMap);
bMCond   = DesDef.bMCond;
iGloNorm = DesDef.iGloNorm;
iGMsca   = DesDef.iGMsca;
GM       = DesDef.GM;
iAdjTo   = DesDef.iAdjTo;

%-Get filenames, build subject & condition indicies
%-----------------------------------------------------------------------
if bMSubj
    nSubj  = spm_input(['number of ',sSubj,'s ?'],'+1','n1');
    bMSubj = nSubj > 1;
else
    nSubj  = 1;
end
guiPos = spm_input('!NextPos');
tmp = [];
for subj  = 1:nSubj
    if bMSubj, strS = [sSubj,' ',int2str(subj),': ']; else, strS = ''; end
    tP = spm_select(Inf,'image',[strS,'select scans...']);
    nt = size(tP,1);
    P     = [P;tP];
    iSubj = [iSubj; subj*ones(nt,1)];
    if bMCond
        str = sprintf('%s[%d] iCond index',strS,nt);
        tmp = spm_input(str,guiPos,'c',tmp',nt);
        iCond = [iCond; tmp];
    else
        iCond = [iCond; ones(nt,1)];
    end
end

%-Total #observations
%-----------------------------------------------------------------------
nScan = length(iSubj);
if nScan==1, error('Only one image - gimme more!'), end

%-Construct H partition
%-----------------------------------------------------------------------
eval(['[H,Hnames] = spm_DesMtx(',HForm,');'])

%-Global normalization options
%-----------------------------------------------------------------------
if length(iGloNorm)>1       %-User has a choice from the options in iGloNorm
    %-Don't offer subject specific AnCova if not bMSubj
    if ~bMSubj, iGloNorm(find(iGloNorm==4))=[]; end
    iGloNorm = spm_input('Select global normalisation','+1','m',...
            sGloNorm(iGloNorm),iGloNorm);
end
sGloNorm = sGloNorm{iGloNorm};

%-Grand mean scaling options
%-----------------------------------------------------------------------
%-Grand mean scaling is implicit in PropSca global normalisation
if iGloNorm==2, iGMsca=4; end
if length(iGMsca)>1     %-User has a choice from the options in iGMsca
    %-Don't offer subject specifics if not bMSubj
    if ~bMSubj, iGMsca(find(iGMsca==3))=[]; end
    iGMsca = spm_input('Grand mean scaling','+1','m',sGMsca(iGMsca),iGMsca);
end
if iGMsca>1 & isempty(GM)   %-Get value for grand mean scaling
    if iGloNorm==2, str='GM: PropSca global mean to ?';
        else,   str='GM: Scale grand mean to ?'; end
    GM = spm_input(str,'+1','e',dGM);
    if GM==0, iGMsca=0; end
elseif iGMsca==1        %-Set GM to zero if not GMscaling
    GM=0;
elseif GM==0;           %-Watch out for GM==0! (=>DesDef is set wrong)
    iGMsca=1;
end
sGMsca = sGMsca{iGMsca};

%-Adjustment options for AnCova designs (for centering of globals)
%-----------------------------------------------------------------------
if any(iGloNorm==[1,2]), iAdjTo=4; end
if length(iAdjTo)>1     %-User has a choice from the options in iAdjTo.
    %-Don't offer subject specifics if not bMSubj
    if ~bMSubj, iAdjTo(find(iAdjTo==2))=[]; end
    iAdjTo=spm_input...
        ('AnCova adjust (centre globals), after any scaling to',...
         '+1','m',sAdjTo(iAdjTo,:),iAdjTo);
end
if iAdjTo==1, aGM = spm_input('AnCova adjust to ?','+1','e',GM); end
sAdjTo = sAdjTo{iAdjTo};



%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm('FigName','AdjMean: configuring',Finter,CmdLine); fprintf('\tconfiguring: ')
spm('Pointer','Watch');

%-Memory map files
%-----------------------------------------------------------------------
V = spm_vol(char(P));

%-Check for consistency of image dimensions and orientation / voxel size
%-----------------------------------------------------------------------
spm_check_orientations(V);


%-Work out required Analyze header info from handles
%-----------------------------------------------------------------------
DIM    = V(1).dim(1:3);
VOX    = sqrt(sum(V(1).mat(1:3,1:3).^2));
ORIGIN = (V(1).mat\[0 0 0 1]')';
ORIGIN = round(ORIGIN(1:3));

%-Compute global values
%-----------------------------------------------------------------------
fprintf('(globals)')
GX     = zeros(nScan,1);
for i  = 1:nScan, GX(i) = spm_global(V(i)); end
fprintf('\b - done)\n')

%-Scaling: compute global scaling factors required to implement proportional
% scaling global normalisation or Grand mean scaling, as requested
%-----------------------------------------------------------------------
rGX = GX;
if iGloNorm==2
    %-Proportional scaling global normalisation
    gSF = GM./GX;
    GX  = GM*ones(nScan,1);
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
else    %-No scaling
    gSF = ones(nScan,1);
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
    else,   error('Illegal iAdjTo')
    end

    if iGloNorm == 3    %-AnCova
        G  = GX - aGM; Gnames = '\zeta';
    elseif iGloNorm == 4    %-AnCova by block/subject
        [G,Gnames] = ...
        spm_DesMtx([iSubj,GX-aGM],'FxC',{'\gamma','\zeta'});
    end
else
    G = []; Gnames = ''; aGM=[];
end

%-Design matrix, parameter estimation matrix
%-----------------------------------------------------------------------
X           = [H G];
[nX,Xnames] = spm_DesMtx('sca',H,Hnames,G,Gnames);
EXnames = spm_DesMtx('ETeXNames',Xnames,aPMap);
iX = struct(    'H',    [1:size(H,2)],  'C',    [],...
        'B',    [],     'G',    size(H,2) + [1:size(G,2)]);
XTXinvX     = inv(X'*X)*X';

%-Contrasts (c) & associated weight matrix (W)
%-----------------------------------------------------------------------
c  = [eye(size(H,2)), zeros(size(H,2),size(G,2))];
nc = size(c,1);
W  = c * XTXinvX;

cNames = spm_DesMtx('Fnames',EXnames(iX.H));
Fnames = cNames;


%=======================================================================
% - D I A G N O S T I C   O U T P U T
%=======================================================================
figure(Fgraph);


%-Display files and variables
%=======================================================================
[P,CPath] = spm_str_manip(P,'c');

%-Display
%-----------------------------------------------------------------------
axes('Position',[0 0 1 .95],'Visible','off')
text(.40,1,'Adjusted means','Fontsize',16,'Fontweight','Bold')
text(.05,.85,'Scan Index','Rotation',90)
if bMSubj, text(.10,.85, sSubj,     'Rotation',90), end
if bMCond, text(.15,.85,'condition',  'Rotation',90), end
x0    = .20; y0 = .83;
dx    = .08; dy = .018;
x     = x0;
text(x + .02,.85,'global','Rotation',90), x = x + 1.5*dx;

text(x,.92,'Base directory:','FontSize',10,'Fontweight','Bold')
text(x,.90,CPath,'FontSize',10)
text(x,.87,'Filename Tails')
y     = y0;
for i = 1:nScan
    text(.03,y,sprintf('%02d :',i));
    if bMSubj, text(.08,y,sprintf('%2d',iSubj(i))), end
    if bMCond, text(.13,y,sprintf('%2d',iCond(i))), end
    x     = x0;
    text(x,y,sprintf('%-8.6g',GX(i)),'FontSize',10), x = x + 1.5*dx;
    text(x,y,P{i},'FontSize',10)
    y     = y - dy;
    if y < 0;
        spm_print
        spm_clf(Fgraph); axis off
        y = y0;
        text(.16,1.02,['Adjusted means (continued)'],...
            'Fontsize',16,'Fontweight','Bold')
    end
end

y      = y - dy;
dy     = dy*1.2;
text(.08,y,['Grand mean scaling: ',sGMsca]), y = y - dy;
text(.08,y,['Global normalisation: ',sGloNorm]), y = y - dy;
text(.08,y,['AnCova adjustment to: ',sAdjTo]), y = y - dy;
text(.08,y,sprintf('Parameters saved to: %s/SPMadj.mat',pwd),'FontSize',8)

spm_print


%-Depict and label design matrix, depict & label contrasts
%=======================================================================
spm_clf(Fgraph); axes('Position',[0 0 1 .95],'Visible','off')
text(.2,1,'Design Matrix, contrasts & contrast files',...
    'Fontsize',16,'Fontweight','Bold');

%-Image scaled design matrix, label the effects, add filenames (if <=32)
%-----------------------------------------------------------------------
hDesMtx = axes('Position',[.07 .5 .6 .3]);
image((nX + 1)*32);
set(hDesMtx,'TickDir','out')
ylabel('scans')
xlabel('parameters')
axes('Position',[.07 .8 .6 .1],'Visible','off')
dx    = 1/size(nX,2);
tXnames = spm_DesMtx('TeXnames',Xnames);
for i = 1:size(nX,2)
    text((i - .5)*dx    ,.1,tXnames{i},'Fontsize',10,'Interpreter','TeX')
    text((i - .5)*dx+.01,.3,EXnames{i},'Fontsize',9,'Rotation',90)
end
if nScan<=32
    set(hDesMtx,'YTick',1:nScan)
    axes('Position',[.68 .5 .3 .3],'Visible','off')
    dy = 1/nScan;
    tP = spm_str_manip(P,'k40');
    for i=1:nScan
        text(0,(nScan-i+0.5)*dy,tP{i},'FontSize',8)
    end
end

%-Depict contrasts and associated (preliminary) filenames
%-----------------------------------------------------------------------
dy = .4/nc;
axes('Position',[.025 .05 .05 .4],'Visible','off')
text(0,.5,'contrasts','HorizontalAlignment','Center','Rotation',90,...
    'FontSize',14,'FontWeight','Bold')
axes('Position',[.6 .44 .40 .02],'Visible','off')
text(0,1,'Contrast files...','FontSize',10,'FontWeight','Bold')
text(0,0,sprintf('...in %s',pwd),'FontSize',8)
hFnames = zeros(1,nc);
for i = 1:nc
    axes('Position',[.1 (.45 -dy*i) .6 .9*dy])
    h = bar(c(i,:),1);
    set(h,'FaceColor',[1 1 1]*.8)
    set(gca,'XLim',[.5,size(c,2)+.5],'Visible','off')
    text(0,0,num2str(i),'HorizontalAlignment','Right','FontSize',10)
    hFnames(i) = text(size(c,2)+.55,.1,Fnames{i},'FontSize',10,...
        'Color',[1 1 1]*.5,'FontAngle','Italic');
end



%=======================================================================
% - W R I T E   C O N T R A S T   I M A G E S
%=======================================================================
fprintf('\tspecify filenames for contrast images...\n')
spm('Pointer');
guiPos = spm_input('!NextPos');
for i = 1:nc
    Fnames{i} = spm_input(sprintf('Contrast %d: filename?',i),...
            guiPos,'s',Fnames{i});
    set(hFnames(i),'string',Fnames{i},'Color','k','FontAngle','Normal')
end
spm_print

%-Parameter images (of interest) - Adjusted mean images
%-----------------------------------------------------------------------
spm('FigName','AdjMean - writing',Finter,CmdLine);
spm('Pointer','Watch');

%-Computation - calculations handled by spm_add.c
%-Using implicit zero masking feature of spm_add: The resultant image
% will be zero at voxels where *any* of the input images are zero.

%-Create handle template for output images as 16bit Int16's
Vo = struct(    'fname',    '',...
        'dim',      V(1).dim(1:3),...
        'dt',           [4 spm_platform('bigend')],...
        'mat',      V(1).mat,...
        'pinfo' ,   [1,0,0]',...
        'descrip',  '');

%-Loop over contrasts
for i = 1:nc
    fprintf('\t...writing image %d/%d: %-20s',i,nc,Fnames{i})
    %-Implement weighted sum by weighting scalefactors in image handles
    w  = W(i,:).*gSF';
    Q  = find(abs(w)>0);
    w  = w(Q); wV = V(Q);
    for j=1:length(Q), wV(j).pinfo(1:2,:)=wV(j).pinfo(1:2,:)*w(j); end
    %-Write header
    Vo.fname   = [Fnames{i},'.img'];
    Vo.descrip = sprintf('Adjusted mean (spm_adjmean) - %s',Fnames{i});
    Vo.pinfo   = [1,0,0]';
    Vo         = spm_create_vol(Vo);
    %-Compute & rewrite header scalefactor
    Vo.pinfo(1)  = spm_add(wV,Vo,'m');
    Vo           = spm_create_vol(Vo);
    fprintf(' (done)\n')
end

%-Prepend PWD to Fnames
Fnames = cellstr([repmat([pwd,filesep],nc,1),char(Fnames)]);


%-Save parameters to SPMadj.mat in current directory
%-----------------------------------------------------------------------
%-Pack design parameters up in a structure
Des = struct(...
        'DesName',  DesName,...
        'HForm',    HForm,...
        'iSubj',    iSubj,...
        'iCond',    iCond,...
        'iGloNorm', iGloNorm,...
        'sGloNorm', sGloNorm,...
        'iGMsca',   iGMsca,...
        'sGMsca',   sGMsca,...
        'GM',       GM,...
        'gSF',      gSF,...
        'iAdjTo',   iAdjTo,...
        'sAdjTo',   sAdjTo,...
        'aGM',      aGM,...
        'X',        X,...
        'nX',       nX,...
        'Xnames',   Xnames,...
        'aPMap',    aPMap,...
        'EXnames',  EXnames,...
        'iX',       iX      );
if spm_check_version('matlab','7') >= 0
    save('SPMadj.mat','-V6','SPMid','DesDef','Des','V','c','cNames','W','Fnames','rGX','GX');
else
    save('SPMadj.mat','SPMid','DesDef','Des','V','c','cNames','W','Fnames','rGX','GX');
end;


%=======================================================================
% - E N D
%=======================================================================
spm('FigName','AdjMean - done',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')
