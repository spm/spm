function spm_adjmean_fmri_ui
% Adjusted means of fMRI via box-car General Linear Model with confounds
% FORMAT spm_adjmean_fmri_ui
%_______________________________________________________________________
%
% spm_adjmean_fmri_ui uses the General Linear Model to produce adjusted
% condition images, adjusted for global effects and confounds.
%
% This program is designed for collapsing data from a single session
% fMRI epoch-related design into a set of representative condition
% images, one for each condition, adjusted for global effects and with
% low frequency drifts removed via a discrete cosine basis set high
% pass filter. The resulting data sets are suitable for a (2nd level)
% random effects analysis of multiple sessions or subjects, or group
% comparisons.
%
% See spm_RandFX.man for further details on implementing random effects
% analyses in SPM96 using a multi-level approach.
%
%
% Overview
% ----------------------------------------------------------------------
% The program works with a single fMRI session, fitting a General
% Linear Model consisting of simple box-cars (optionally convolved with
% an estimated haemodynamic response function) and an (optional) high
% pass filter of discrete cosine basis functions. The effects of
% interest (the box-cars) are orthogonalised (residualised) with
% respect to the confounds (the high pass filter) to ensure that *all*
% confounds are removed from the data, even if they're correlated with
% the effects of interest. (For well designed experiments this makes
% little difference.) Proportional scaling and AnCova global
% normalisation are supported, the latter including the option to scale
% all the images such that their the grand mean (GM) is a specified
% value.
%
% The interface is similar to a cut-down SPM-fMRI, and the adjusted
% means are the parameter estimates from the model. The user is first
% prompted to select the scans for a single fMRI session. Then the
% epoch condition order is specified. This should be a r-vector, where
% r is the number of epochs, of integers 1:n or 0:n-1 where n is the
% number of conditions (0 can be used to indicate the baseline or
% "rest" condition. Then the number of scans per epoch is specified:
% This can be a single integer (all epochs have the same number of
% scans), or an r-vector of integers specifying number of scans for the
% corresponding epoch.
%
% Once the experimental design has been specified, the user is given
% various options: The box-cars can be convolved with an approximate
% haemodynamic reponse function; High-pass filter components can be
% added to the confounds (at a user-specified cut-off which defaults to
% twice the maximum period of the experiment); Global normalisation can
% be omitted, or implemented by proportional scaling or AnCova; and the
% grand mean scaling options specified.
%
% With the design and adjustments specified, the model is constructed,
% and the user prompted to enter/confirm filenames for the adjusted
% condition images.
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
% The default value is 100.
%
% If computing adjusted means for subsequent (2nd level) modelling, as
% with a random effects analysis, then it is important to use a
% seperable model, such that the adjustment for one subject is
% independent of other subjects entered into the model. Thus,
% proportional scaling or subject-specific AnCova adjustment must be
% used. Further, multiple runs *must* use the same GM value, and should
% scale Grand mean *by subject*.
%
% ( A separate program (spm_adjmean_ui) is available for computing       )
% ( adjusted condition means of PET data. The functionality is similar   )
% ( to this code, but the two routines have been separated for           )
% ( algorithmic clarity.                                                 )
%
% Diagnostic output
% ----------------------------------------------------------------------
% Diagnostic output consists of two sections:
%
% The first page lists the filenames, various parameters (Grand mean
% scaling etc.), and gives a plot of the image global means against
% scan number, overlaid on an "image" of the condition effects. Watch
% out for condition dependent global changes!
%
% The second part is a single page depicting the design matrix, effect
% names, parameter contrasts used, and the corresponding image files
% written.
%
% As always, look at the resulting mean images to make sure they look OK!
%
%
% Algorithm
% ----------------------------------------------------------------------
% The model at each voxel is Y = X*B + e, with a set of least squares
% estimates for the vector of parameters B as b = pinv(X)*Y. For c a
% vector of contrast weights extracting the appropriate parameter, the
% contrast of the parameter estimates is c'*b = c'*pinv(X)*Y, a
% weighted sum (or weighted mean) of the data at that voxel. These
% weights are identical for all voxels, so the image of the parameter
% estimate can be computed as a weighted mean of the images.
%
% The design matrix is split into effects of interest [C], a constant
% term [B], and confounds [G]. The columns of G are centered so that
% the confound cannot model any of the mean. The effects of interest
% are orthogonalised wirit. the confounds, using C=C-G*pinv(G)*C; This
% ensures that *all* confound effects are removed from the data, even
% if they are correlated with the effects of interest.
%
% Once the weights have been worked out for each adjusted mean image,
% computation proceeds by passing appropriate weights and image
% filenames to spm_mean, which writes out the appropriate parameter
% image as an Analyze format image of the same type (see spm_type) as
% the input images.
%
%
% Variables saved in SPMadj.mat data file
% ----------------------------------------------------------------------
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
% P     String matrix of filenames
% iCond     Condition indicator vector
% iGloNorm  Global normalisation option
% sGloNorm  Global normalisation description
% iGMsca    Grand mean scaling option
% sGMsca    Grand mean scaling description
% HPFc      High pass filter cut-off period (s)
% HPF       High pass filter
% sHPF      Description of high-pass filter
% rC        raw C partition of design matrix, prior to orthogonalisation
% C     C (covariates of interest) partition of design matrix
% Cnames    Names of parameters corresponding to columns of C
% B     B (block) partition of the design matrix
% Bnames    Names of parameters corresponding to columns of B
% G     G (confounding covariates) partition of design matrix
% Gnames    Names of parameters corresponding to columns of G
% rX        raw design matrix, prior to orthogonalisation of C partition
% X     design matrix (=[C,B,G])
% nrX       raw design matrix, normalised for display
% nX        design matrix, normalised for display
% c     Matrix of contrasts, contrasts in rows
% cNames    Names associated with contrasts
% W     Weights for images corresponding to contrasts
% CWD       Current Working Directory (when run)
% Fnames    Filenames of adjusted mean images written
% rGX       raw global means (before any scaling)
% gSF       Image scale factors for global scaling (inc. grand mean scaling)
% GX        Global means after scaling
% GM        Grans Mean used for scaling
%
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_adjmean_fmri_ui.m 4418 2011-08-03 12:00:13Z guillaume $



%=======================================================================
% - S E T U P
%=======================================================================
SCCSid = '2.9';
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','AdjMean/fMRI',1);
spm_help('!ContextHelp',[mfilename,'.m'])


%=======================================================================
% - D E S I G N   P A R A M E T E R S
%=======================================================================

%-Get filenames
%-----------------------------------------------------------------------
P     = spm_select(Inf,'image','select scans (single session)');
nScan = size(P,1);
if nScan==1, error('Only one image - gimme more!'), end

%-Get condition information & construct vector of conditions indices
%-----------------------------------------------------------------------
%-Epoch condition order
eCond = spm_input('epoch order Eg 121... {0=null}',1,'c')';
Conds = unique(eCond);

%-Epoch length(s)
eLen   = [];
guiPos = spm_input('!NextPos');
while any(~eLen) | sum(eLen) ~= nScan
    eLen = spm_input('#scans/epoch Eg 8 or 8 6 8... ',guiPos);
    if length(eLen)==1, eLen=eLen*ones(1,length(eCond)); end
end
eLen    = eLen(:)';

%-Epoch onsets, in scans, starting at 1
eOns = cumsum(eLen) - eLen + 1;

%-iCond condition indicator vector
iCond = zeros(nScan,1);
for i=1:length(eLen), iCond(eOns(i)+[1:eLen(i)]-1)=eCond(i); end


%-Get Repeat time, construct approximate haemodynamic response function
%-----------------------------------------------------------------------
RT  = spm_input('interscan interval (in seconds)','+1');
hrf = spm_hrf(RT);


%-Construct design matrix C partition
% (box-cars convolved with approximate HRF if requested)
% null condition with iCond==0 is ignored by spm_DesMtx with '~' constraint
%-----------------------------------------------------------------------
[C,Cnames,Ci] = spm_DesMtx(iCond,'~','\alpha');
% C = spm_detrend(C);

%-Convolve with hemodynamic response function, if requested
bHRFconv = spm_input('Conv. box-cars w/ approx HRF?','+1','y/n',[1,0],1);
if bHRFconv
    d = length(hrf);
    C = [ones(d,1)*C(1,:); C];
    C = spm_sptop(hrf,nScan + d,1)*C;
    C = C([1:nScan]+d,:);
end


%-Construct design matrix block (B) partition
%-----------------------------------------------------------------------
B = ones(nScan,1); Bnames = '\mu';


%-Construct design matrix confound (G) partition?
%-----------------------------------------------------------------------
%*** Allow user specified confounds?
%*** ...and user specified convolution of confounds with approximate HRF?
G = []; Gnames = '';


%-Add high-pass filter using discrete cosine set
%-----------------------------------------------------------------------
bHPF = spm_input('Use high pass filter?','+1','y/n',[1 0],1);
if bHPF
    HPFc=nScan;
    for i=Conds', HPFc = min([HPFc,max(diff(eOns(eCond==i)))]); end
    HPFc=spm_input('HPF cut-off period {in seconds}','0','e',2*HPFc*RT);
    %-Find max order for discrete cosine set, HPFk
    % (from period>HPFc, period=nScan*RT/(k/2) for order k; k<nScan/2)
    HPFk = fix(min(2*(nScan*RT)/HPFc,nScan/2));
    HPF = []; HPFnames = '';
    for k = 1:HPFk
        HPF = [HPF, cos(k*pi*([1:nScan]-1)'/(nScan-1))];
        HPFnames = ...
            strvcat(HPFnames,sprintf('LowHz (%4.0fs)',2*RT*nScan/k));
    end
    HPF = HPF - ones(nScan,1)*mean(HPF);    %-Mean correct (by column)
    sHPF = sprintf('High-pass filter of %d components, cut off = %ds',...
        HPFk,HPFc);

    G = [G, HPF]; Gnames = strvcat(Gnames,HPFnames);
else
    HPFc=0; HPF=[]; sHPF='No high-pass filter';
end

%-Global normalization options
%-----------------------------------------------------------------------
sGloNorm = strvcat('None','Proportional scaling','AnCova');
iGloNorm = spm_input('Select global normalisation',...
    '+1','m',sGloNorm,[],2);
sGloNorm = deblank(sGloNorm(iGloNorm,:));


%-Grand mean scaling
%-----------------------------------------------------------------------
sGMsca = strvcat('None','Scaling of overall Grand Mean',...
    '(Implicitly via PropSca global normalisation)');
if iGloNorm==2, iGMsca=3; else, iGMsca=2; end
sGMsca = deblank(sGMsca(iGMsca,:));
GM = 100;


%-Temporal smoothing
%-----------------------------------------------------------------------
%**** What to do with temporal smoothing?
%SIGMA = spm_input('Temporal smoothing FWHM {secs}','+1','e',6);
%SIGMA = SIGMA/sqrt(8*log(2))/RT;



%=======================================================================
% - C O M P U T A T I O N
%=======================================================================
spm('FigName','AdjMean/fMRI: configuring',Finter,CmdLine);
fprintf('\tconfiguring: ')
spm('Pointer','Watch');

%-Memory map files
%-----------------------------------------------------------------------
V = spm_vol(char(P));
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
    sGloNorm = sprintf('%s, to %g',sGloNorm,GM);
elseif iGMsca==2
    %-Grand mean scaling (overall)
    gSF = GM/mean(GX);
    GX  = GX*gSF;
    sGMsca = sprintf('%s, to %g',sGMsca,GM);
else    %-No scaling
    gSF = ones(nScan,1);
end


%-AnCova options: Construct Global covariates of no interest partition
%-----------------------------------------------------------------------
if iGloNorm==3
    G      = [G, GX-mean(GX)];
    Gnames = strvcat(Gnames,'Global');
    sGloNorm = [sGloNorm,', to Grand Mean'];
end

%-Design matrix - raw & with effects of interest orthogonalised wirit G
%-----------------------------------------------------------------------
rC = C;
rX = [rC B G];
[nrX,Xnames] = spm_DesMtx('sca',C,Cnames,B,Bnames,G,Gnames);
if size(G,2)
    C  = C - G*pinv(G)*C;
end
X  = [C B G];
nX = spm_DesMtx('sca',X,Xnames);
EXnames = spm_DesMtx('ETeXNames',Xnames);
iX = struct(    'H',    [],     'C',    [1:size(C,2)],...
        'B',    size(C,2)+1,    'G',    size(C,2) +1 +[1:size(G,2)]);

%-Temporal smoothing?
%-----------------------------------------------------------------------
%**** Include temporal smoothing? Theoretically, shouldn't make any
% difference to parameter estimates, but would improve efficiency and
% match SPM-fMRI.
%**** Q:Orthogonalisation of C wirit G pre or post temporal smoothing?

%-Parameter estimation matrix & data weights matrix for contrasts
%-----------------------------------------------------------------------
pX  = pinv(X);

%-Contrasts (c) & associated weight matrix (W)
%-----------------------------------------------------------------------
nc     = size(C,2);
c      = [eye(nc), ones(nc,size(B,2)), zeros(nc,size(G,2))];
cNames = spm_DesMtx('Fnames',EXnames(iX.C));
%-if we have a null condition, then constant *is* null condition effect
if any(Conds==0)
    c      = [zeros(1,nc), 1, zeros(1,size(G,2)); c];
    cNames = [{'baseline'};cNames];
    nc     = nc+1;
end

W   = c * pX;
Fnames = cNames;


%=======================================================================
% - D I A G N O S T I C   O U T P U T
%=======================================================================
fprintf('\tdisplaying diagnostic output...\n')
figure(Fgraph);


%-Display files and variables
%=======================================================================

%-Display parameters
%-----------------------------------------------------------------------
figure(Fgraph); spm_clf; axis off
text(0.30,1.00,'Adjusted means (fMRI)','Fontsize',16,'Fontweight','Bold')
text(-.10,0.96,'Files:','Fontweight','Bold')
text(-.10,0.94,sprintf('1 to %4d: %s...',nScan,P{1}),'FontSize',10)
text(-.10,0.92,sprintf('Repeat time = %.4f seconds',RT),'FontSize',10)
text(-.10,0.88,'Image Global means (& conditions):','Fontweight','Bold')

text(-.1,0.20,['Grand mean scaling: ',sGMsca])
text(-.1,0.17,['Global normalisation: ',sGloNorm])
text(-.1,0.14,sHPF)
text(-.1,0.11,'Effects of interest orthogonalised (residualised) to confounds')
text(-.1,0.08,'No temporal smoothing (only required for inference)')
text(-.1,0.05,sprintf('Parameters saved to: %s/SPMadj.mat',pwd),'FontSize',10)

%-Global mean values plot, with underlayed conditions
%-----------------------------------------------------------------------
%-Conditions plot (for comparison with global means plot)
axes('Position',[0.1,0.5,0.8,0.3])
image(32+(rC'-min(rC(:)))*32*(max(rC(:))-min(rC(:))))
set(gca,'Visible','off')

axes('Position',[0.1,0.5,0.8,0.3])
if iGloNorm==2
    plot(rGX)
    ylabel('global mean (pre-scaling)')
else
    plot(GX)
    ylabel('global mean (post-scaling)')
end
set(gca,'XLim',[0.5,nScan+0.5],'Color','none')
xlabel('scan index')

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
    text((i - .5)*dx+.01,.3,EXnames{i},'Fontsize',9,'Rotation',90)
    if ~strcmp(Xnames{i},EXnames{i})
        text((i - .5)*dx,.1,tXnames{i},'Fontsize',10,'Interpreter','TeX')
    end
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
spm('FigName','AdjMean/fMRI - writing',Finter,CmdLine);
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
    Vo.descrip = sprintf('Adjusted mean (spm_adjmean_fmri) - %s',Fnames{i});
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
        'DesName',  'Single sess. fMRI box-cars',...
        'HForm',    [],...
        'iSubj',    [],...
        'iCond',    iCond,...
        'iGloNorm', iGloNorm,...
        'sGloNorm', sGloNorm,...
        'iGMsca',   iGMsca,...
        'sGMsca',   sGMsca,...
        'GM',       GM,...
        'gSF',      gSF,...
        'iAdjTo',   [],...
        'sAdjTo',   [],...
        'aGM',      [],...
        'X',        X,...
        'nX',       nX,...
        'Xnames',   Xnames,...
        'aPMap',    [],...
        'EXnames',  EXnames,...
        'iX',       iX      );
if spm_check_version('matlab','7') >= 0
    save('SPMadj.mat','-V6',...
        'SPMid','Des','V','c','cNames','W','Fnames','HPFc','HPF','sHPF','rX','nrX','rGX','GX');
else
    save('SPMadj.mat',...
        'SPMid','Des','V','c','cNames','W','Fnames','HPFc','HPF','sHPF','rX','nrX','rGX','GX');
end;


%=======================================================================
% - E N D
%=======================================================================
spm('FigName','AdjMean/fMRI - done',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')
