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
% **** under construction ****
%
%
% Platform
% ----------------------------------------------------------------------
% This version was written for MatLab4.2c with SPM'96 (spm_mean.m v1.10)
%
%_______________________________________________________________________
% %W% Andrew Holmes %E%


%=======================================================================
% - S E T U P
%=======================================================================
fprintf('\nSPM: spm_adjmean_fmri_ui\n')
fprintf('%c','='*ones(1,72)),fprintf('\n')
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), Finter=spm('CreateIntWin'); end
spm_clf(Finter)
set(Finter,'Name','AdjMean/fMRI')
spm_help('!ContextHelp','spm_adjmean_fmri_ui.m')


%=======================================================================
% - D E S I G N   P A R A M E T E R S
%=======================================================================

%-Get filenames
%-----------------------------------------------------------------------
P = spm_get(Inf,'.img','select scans (single session)');
n = size(P,1);

%-Get condition information & construct vector of conditions indices
%-----------------------------------------------------------------------
%-Epoch condition order
a     = spm_input('epoch order eg 1 2 1..{0 = null}',1);
a     = a(:)';
nCond = max(a);

%-Epoch length(s)
e     = [];
guiPos = '+1';
while any(~e) | sum(e) ~= n
	e   = spm_input('scans per epoch eg 8 or 8 6... ',guiPos);
	if prod(size(e))==1, e=e*ones(1,length(a)); end
	guiPos='0';
end
e     = e(:)';

% epoch onsets, in scans, starting at 0
ons = cumsum(e) - e;

%-iCond condition indicator vector
iCond = zeros(n,1);
for i = 1:length(e), iCond([1:e(i)]+ons(i)) = a(i)*ones(e(i),1); end


%-Get Repeat time, construct approximate haemodynamic response function
%-----------------------------------------------------------------------
RT  = spm_input('interscan interval (in seconds)','+1');
hrf = spm_hrf(RT);


%-Construct design matrix C partition
% (box-cars convolved with approximate HRF if requested)
%-----------------------------------------------------------------------
%-Design matrix (box-cars) - omit condition `0' as null (baseline)
[rC,Cnames,Ci] = spm_DesMtx(iCond,'-','cond');
if Ci(1)==0, rC(:,1)=[]; Cnames(1,:)=[]; Ci(1)=[]; end
% rC = spm_detrend(rC);

%-Convolve with hemodynamic response function, if requested
C = rC;
bHRFconv = spm_input('Conv. box-cars with approx HRF?','+1','y/n',[1,0],1);
if bHRFconv
	d = length(hrf);
	C = [ones(d,1)*C(1,:); C];
	C = spm_sptop(hrf,n + d,1)*C;
	C = C([1:n]+d,:);
end


%-Construct design matrix block (B) partition
%-----------------------------------------------------------------------
B = ones(n,1); Bnames = 'Constant';


%-Construct design matrix confound (G) partition?
%-----------------------------------------------------------------------
%*** Allow user specified confounds?
%*** ...and user specified convolution of confounds with approximate HRF?
G = []; Gnames = '';


%-Add high-pass filter using discrete cosine set
%-----------------------------------------------------------------------
HPFc = n; for i = 0:nCond, HPFc = min([HPFc,max(diff(ons(a == i)))]); end
bHPF = spm_input('Use high pass filter?','+1','y/n',[1 0],1);
if bHPF
	HPFc = spm_input('HPF cut-off period {in seconds}','0','e',2*HPFc*RT);
	%-Find max order for discrete cosine set, HPFk
	% (from period>HPFc, period=n*RT/(k/2) for order k; k<n/2)
	HPFk = fix(min(2*(n*RT)/HPFc,n/2));
	HPF = [];
	for k = 1:HPFk
		HPF = [HPF, cos(k*pi*([1:n]-1)'/(n-1))];
	end
	HPF = HPF - ones(n,1)*mean(HPF);	%-Mean correct (by column)
	sHPF = sprintf('High-pass filter of %d components, cut off = %ds',...
		HPFk,HPFc);

	%-Construct effect names
	HPFnames = '';
	for i = 1:size(HPF,2);
		HPFnames = str2mat(HPFnames,sprintf('Low Hz (%d)',i));
	end
	HPFnames(1,:)=[];
	G = [G, HPF]; Gnames = str2mat(Gnames,HPFnames);
else
	HPFc=0; HPF=[]; sHPF='No high-pass filter';
end

%-Global normalization options
%-----------------------------------------------------------------------
sGloNorm = str2mat('None','Proportional scaling','AnCova');
iGloNorm = spm_input('Select global normalisation',...
	'+1','m',sGloNorm,[],2);
sGloNorm = deblank(sGloNorm(iGloNorm,:));


%-Grand mean scaling
%-----------------------------------------------------------------------
sGMsca = str2mat('None','Scaling of overall Grand Mean',...
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
set(Finter,'Name','AdjMean/fMRI - computing','Pointer','watch')
fprintf('\tcomputing: ')

%-Get file identifiers
%-----------------------------------------------------------------------
V     = zeros(12,n);
for i = 1:n; V(:,i) = spm_map(P(i,:));  end

%-Check for consistency of image size and voxel size
%-----------------------------------------------------------------------
if ~all(all(~diff(V([1:6],:)')))
	error('data do not have the same image and voxel size'); end

%-Get ORIGIN from first image
%-----------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P(1,:));

%-Compute global values
%-----------------------------------------------------------------------
fprintf('(globals)')
GX     = zeros(n,1);
for i  = 1:n, GX(i) = spm_global(V(:,i)); end
fprintf('\b - done)\n')

%-Save image scalefactors
%-----------------------------------------------------------------------
iSF = V(7,:)';


%-Scaling: compute global scaling factors required to implement proportional
% scaling global normalisation or Grand mean scaling, as requested
%-----------------------------------------------------------------------
rGX = GX;
if iGloNorm==2
	%-Proportional scaling global normalisation
	gSF = GM./GX;
	GX  = GM*ones(n,1);
	sGloNorm = sprintf('%s, to %g',sGloNorm,GM);
	%** scale rGX for printout? ...or just graph them?
elseif iGMsca==2
	%-Grand mean scaling (overall)
	gSF = GM/mean(GX);
	GX  = GX*gSF;
	sGMsca = sprintf('%s, to %g',sGMsca,GM);
else	%-No scaling
	gSF = ones(n,1);
end


%-AnCova options: Construct Global covariates of no interest partition
%-----------------------------------------------------------------------
if iGloNorm==3
	G      = [G, GX-mean(GX)];
	Gnames = str2mat(Gnames,'Global');
	sGloNorm = [sGloNorm,', to Grand Mean'];
end
%-Clean blank top line from Gnames (if necessary)
if size(Gnames,1), if all(Gnames(1,:)==' '), Gnames(1,:)=[]; end, end

%-Design matrix - raw & with effects of interest orthogonalised wirit G
%-----------------------------------------------------------------------
rC = C;
rX = [rC B G];
[nrX,Xnames] = spm_DesMtxSca(C,Cnames,B,Bnames,G,Gnames);
if size(G,2)
	C  = C - G*pinv(G)*C;
end
X  = [C B G];
nX = spm_DesMtxSca(X,Xnames);

%-Temporal smoothing?
%-----------------------------------------------------------------------
%**** Include temporal smoothing? Theoretically, shouldn't make any
% difference to parameter estimates, but would improve efficiency and
% match SPM-fMRI.
%**** Q:Orthogonalisation of C wirit G pre or post temporal smoothing?


%-Contrasts (x)
%-----------------------------------------------------------------------
nx     = size(C,2);
x      = [eye(nx), ones(nx,size(B,2)), zeros(nx,size(G,2))];
xNames = Cnames;
%-if we have a null condition, then constant *is* null condition effect
if min(iCond)==0
	x      = [zeros(1,nx), 1, zeros(1,size(G,2)); x];
	xNames = str2mat('baseline',xNames);
	nx     = nx+1;
end

%-Parameter estimation matrix & data weights matrix for contrasts
%-----------------------------------------------------------------------
pX  = pinv(X);
W   = x * pX;


%-Parameter images (of interest) - Adjusted mean images
%-----------------------------------------------------------------------
set(Finter,'Name','AdjMean/fMRI - writing','Pointer','arrow')
fprintf('\twriting parameter images... enter filenames...\n')

%-Computation: calculations handled by spm_mean.c
% Relying on all images used to make each contrast image being masked
% identically, with out of mask values set to zero. This is the case
% for sets of images realigned with SPM'96 realignment.

Fnames = []; CWD = pwd;
guiPos = '+1'; %*** guiPos = spm_input('!NextPos');
for i = 1:nx
	Fn = deblank(xNames(i,:)); Fn=Fn(Fn~=' ');
	Fn = spm_input(sprintf('Contrast %d: filename?',i),guiPos,'s',Fn);
	%*** make sure Fn is valid? Allow overwriting?
	Fnames = str2mat(Fnames,Fn);
	w  = W(i,:).*gSF'.*iSF';
	Q  = find(abs(w)>0);
	%wV = V(:,Q);	wV(7,:) = w(Q);
	fprintf('\t...writing image %d: %-20s',i,Fn)
	%sf = spm_mean(wV,[Fn,'.img']);
	sf = spm_mean(prod(DIM),TYPE,[Fn,'.img'],P(Q,:),w(Q));
	str = sprintf('Adjusted mean (spm_adjmean) - %s',Fn);
	%spm_hwrite([Fn,'.hdr'],DIM,VOX,sf,spm_type('int16'),OFFSET,ORIGIN,str);
	spm_hwrite([Fn,'.hdr'],DIM,VOX,sf,TYPE,OFFSET,ORIGIN,str);
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
	P iCond ...
	iGloNorm sGloNorm iGMsca sGMsca ...
	HPFc HPF sHPF ...
	rC C Cnames B Bnames G Gnames ...
	rX X nrX nX ...
	x xNames W ...
	CWD Fnames ...
	rGX gSF GX GM


%=======================================================================
% - D I A G N O S T I C   O U T P U T
%=======================================================================
set(Finter,'Name','AdjMean/fMRI - done')
fprintf('\tdisplaying diagnostic output...\n')

Fgraph = spm_figure('FindWin','Graphics');
if isempty(Fgraph), Fgraph = spm_figure('Create','Graphics','Graphics'); end


%-Display files and variables
%=======================================================================

%-Display parameters
%-----------------------------------------------------------------------
figure(Fgraph); spm_clf; axis off
text(0.30,1.00,'Adjusted means (fMRI)','Fontsize',16,'Fontweight','Bold')
text(-.10,0.96,'Files:','Fontweight','Bold')
text(-.10,0.94,sprintf('1 to %4d: %s...',n,deblank(P(1,:))),'FontSize',10)
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
set(gca,'XLim',[0.5,n+0.5])
xlabel('scan index')

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
dy = 0.4/nx;
axes('Position',[0.025 0.05 0.05 0.4],'Visible','off')
text(0,0.5,'contrasts','HorizontalAlignment','Center','Rotation',90,...
	'FontSize',14,'FontWeight','Bold')
axes('Position',[0.6 0.44 0.40 0.02],'Visible','off')
text(0,1,'Contrast files...','FontSize',10,'FontWeight','Bold')
text(0,0,sprintf('...in %s',CWD),'FontSize',8)
for i = 1:nx
	axes('Position',[0.1 (0.45 -dy*i) 0.6 0.9*dy])
	[tx ty] = bar(x(i,:));
	fill(tx,ty,[1 1 1]*.8);
	% h = bar(x(i,:),1);
	% set(h,'FaceColor',[1 1 1]*.8)
	set(gca,'XLim',[0.5,size(x,2)+0.5],'Visible','off')
	text(0,0,num2str(i),'HorizontalAlignment','Right','FontSize',10)
	text(size(x,2)+.55,0,sprintf('%s',deblank(Fnames(i,:))),'FontSize',10)
end

spm_print


%-END
%=======================================================================
fprintf('\n\n')
