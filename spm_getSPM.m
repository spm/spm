function [SPM,VOL,xX,xCon,xSDM] = spm_getSPM
% computes a specified and thresholded SPM/PPM following parameter estimation
% FORMAT [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
%
% SPM    - structure containing SPM, distribution & filtering details
% .swd   - SPM working directory - directory containing current SPM.mat
% .title - title for comparison (string)
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X, F or P}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .Ic    - indices of contrasts (in xCon)
% .Im    - indices of masking contrasts (in xCon)
% .pm    - p-value for masking (uncorrected)
% .Ex    - flag for exclusive or inclusive masking
% .u     - height threshold
% .k     - extent threshold {voxels}
% .XYZ   - location of voxels {voxel coords}
% .XYZmm - location of voxels {mm}
% .QQ    - indices of volxes in Y.mad file
%
%
% VOL    - structure containing details of volume analysed
% .S     - search Volume {voxels}
% .R     - search Volume {resels}
% .FWHM  - smoothness {voxels}     
% .M     - voxels - > mm matrix
% .iM    - mm -> voxels matrix
% .VOX   - voxel dimensions {mm} - column vector
% .DIM   - image dimensions {voxels} - column vector
% .Vspm  - Mapped statistic image(s)
% .Msk   - Mask in 1 of 3 forms
%            o Scalar, indicating implicit mask value in statistic image(s)
%            o Vector of indicies of elements within mask
%            o Mapped mask image
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T', 'F' or 'P')
% .c     - Contrast weights (column vector contrasts)
% .X0    - Reduced design matrix data (spans design space under Ho)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X0 = spm_FcUtil('X0',...
% .iX0   - Indicates how contrast was specified:
%          If by columns for reduced design matrix then iX0 contains the
%          column indices. Otherwise, it's a string containing the
%          spm_FcUtil 'Set' action: Usuall one of {'c','c+','X0'}
% .X1o   - Remaining design space data (X1o is orthogonal to X0)
%          Stored as coordinates in the orthogonal basis of xX.X from spm_sp
%          (Matrix in SPM99b)  Extract using X1o = spm_FcUtil('X1o',...
% .eidf  - Effective interest degrees of freedom (numerator df)
%        - Or effect-size threshold for Posterior probability
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SPM image
%
% xSDM   - structure containing contents of SPM.mat file
%          ( see spm_spm.m for contents...
%          ( ...except that xX & XYZ are removed, being returned separately
%          ( ...and xSDM.Vbeta & xSDM.VResMS are re-mmapped
%
% In addition, the xCon.mat file is updated. For newly evaluated
% contrasts, SPM images (spmT_????.{img,hdr}) are written, along with
% contrast (con_????.{img,hdr}) images for SPM{T}'s, or Extra
% Sum-of-Squares images (ess_????.{img,hdr}) for SPM{F}'s.
% 
% The contrast images are the weighted sum of the parameter images,
% where the weights are the contrast weights, and are uniquely
% estimable since contrasts are checked for estimability by the
% contrast manager. These contrast images (for appropriate contrases)
% are suitable summary images of an effect at this level, and can be
% used as input at a higher level when effecting a random effects
% analysis. (Note that the ess_????.{img,hdr} and
% SPM{T,F}_????.{img,hdr} images are not suitable input for a higher
% level analysis.) See spm_RandFX.man for further details.
%
%_______________________________________________________________________
%
% spm_getSPM prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects by subsequent routines.
% 
% For general linear model Y = XB + E with data Y, desgin matrix X,
% parameter vector B, and (independent) errors E, a contrast c'B of the
% parameters (with contrast weights c) is estimated by c'b, where b are
% the parameter estimates given by b=pinv(X)*Y.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts. Contrasts are estimable linear combinations of the
% parameters, and are specified using the SPM contrast manager
% interface [spm_conman.m]. SPMs are generated for the null hypotheses
% that the contrast is zero (or zero vector in the case of
% F-contrasts). See the help for the contrast manager [spm_conman.m]
% for a further details on contrasts and contrast specification.
% 
% A conjunction assesses the conjoint expression of two or more
% effects. The conjunction SPM is the minimum of the component SPMs
% defined by the multiple contrasts. The distributional results used
% for minimum fileds require the SPMs to be identically distributed and
% independent. Thus, all component SPMs must be either SPM{t}'s, or
% SPM{F}'s with the same degrees of freedom. Independence is roughly
% guaranteed for large degrees of freedom (and independent data) by
% ensuring that the contrasts are "orthogonal". Note that it is *not*
% the contrast weight vectors themselves that are required to be
% orthogonal, but the subspaces of the data space implied by the null
% hypotheses defined by the contrasts (c'pinv(X)).
% 
% To ensure approximate independence of the component SPMs in a
% conjunction, non-orthogonal contrasts for a conjunction are
% successively enforced to be orthogonal. If this is required, you will
% be asked to specify the orthogonalisation order (as a permutation of
% the contrast indices). The contrasts are then serially orthogonalised
% in this order, possibly generating new contrasts, such that the
% second is orthogonal to the first, the third to the first two, and so
% on.
%
% Masking simply eliminates voxels from the current contrast if they
% do not survive an uncorrected p value (based on height) in one or
% more further contrasts.  No account is taken of this masking in the
% statistical inference pertaining to the masked contrast.
% 
% The SPM is subject to thresholding on the basis of height (u) and the
% number of voxels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent. If you want to see all voxels simply enter 0.  In this
% instance the 'set-level' inference can be considered an 'omnibus test'
% based on the number of clusters that obtain.
%
% BAYESIAN INFERENCE AND PPMS - POSTERIOR PROBABILITY MAPS
%
% If conditional estimates are available (and your contrast is a T
% contrast) then you are asked whether the inference should be 'Bayesian'
% or 'classical' (using GRF).  If you choose Bayesian the contrasts are of
% conditional (i.e. MAP) estimators and the inference image is a
% posterior probability map (PPM).  PPMs encode the probability that the
% contrast exceeds some specified threshold.  This threshold is stored in
% the xCon.eidf.  Subsequent plotting and tables will use the conditional
% estimates and associated posterior or conditional probabilities.
% 
% see spm_results_ui.m for further details of the SPM results section.
%_______________________________________________________________________
% %W% Andrew Holmes, Karl Friston & Jean-Baptiste Poline %E%

SCCSid   = '%I%';

%-GUI setup
%-----------------------------------------------------------------------
SPMid  = spm('SFnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',mfilename)

%-Select SPM.mat & note SPM results directory
%-----------------------------------------------------------------------
swd    = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');

%-Preliminaries...
%=======================================================================

%-Load SPM.mat
%-----------------------------------------------------------------------
xSDM   = load(fullfile(swd,'SPM.mat'));	%-Contents of SPM.mat (in a struct)


%-Get Stats data from SPM.mat
%-----------------------------------------------------------------------
xX     = xSDM.xX;			%-Design definition structure
XYZ    = xSDM.XYZ;			%-XYZ coordinates
S      = xSDM.S;			%-search Volume {voxels}
R      = xSDM.R;			%-search Volume {resels}
xSDM   = rmfield(xSDM,{'xX','XYZ'});	%-Remove (large) duplicate fields

%-Get/Compute mm<->voxel matrices & image dimensions from SPM.mat
%-----------------------------------------------------------------------
 M     = xSDM.M;
iM     = inv(M);
DIM    = xSDM.DIM;

%-Build index from XYZ into corresponding Y.mad locations
%-----------------------------------------------------------------------
if exist(fullfile(swd,'Yidx.mat'),'file') & exist(fullfile(swd,'Y.mad'),'file')
	load(fullfile(swd,'Yidx.mat'))
	QQ = sparse(1,Yidx,1:length(Yidx),1,S);
else
	QQ = sparse(1,S);
end

%-Contrast definitions
%=======================================================================

%-Load contrast definitions (if available)
%-----------------------------------------------------------------------
if exist(fullfile(swd,'xCon.mat'),'file')
	load(fullfile(swd,'xCon.mat'))
else
	xCon = [];
end

%-See if can write to current directory (by trying to open for appending)
%-----------------------------------------------------------------------
wOK = 1;
fid = fopen(fullfile(swd,'xCon.mat'),'a');
if fid < 0
   wOK = 0;
   str = {	'Can''t write to the results directory:',...
         '(problem saving xCon.mat)',...
         ['        ',swd],...
         ' ','-> results restricted to contrasts already computed'};
   spm('alert!',str,mfilename,1);
else
   fclose(fid);
end


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================

%-Get contrasts
%-----------------------------------------------------------------------
[Ic,xCon] = spm_conman(xX,xCon,'T&F',Inf,...
	'	Select contrasts...',' for conjunction',wOK);


%-Enforce orthogonality of multiple contrasts for conjunction
% (Orthogonality within subspace spanned by contrasts)
%-----------------------------------------------------------------------
% (Don't ask orthogonalisation order if not needed.)
if length(Ic) > 1 & ~spm_FcUtil('|_?',xCon(Ic), xX.xKXs)

    %-Get orthogonalisation order from user
    %-------------------------------------------------------------------
    Ic = spm_input('orthogonlization order','+1','p',Ic,Ic);
    
    %-Successively orthogonalise
    %-NB: This loop is peculiarly controlled to account for the
    %     possibility that Ic may shrink if some contrasts diasppear
    %     on orthogonalisation (i.e. if there are colinearities)
    %-------------------------------------------------------------------
    i = 1; while(i < length(Ic)), i = i + 1;
    
	%-Orthogonalise (subspace spanned by) contrast i wirit previous
	%---------------------------------------------------------------
	oxCon = spm_FcUtil('|_',xCon(Ic(i)), xX.xKXs, xCon(Ic(1:i-1)));
    
	%-See if this orthogonalised contrast has already been entered
	% or is colinear with a previous one. Define a new contrast if
	% neither is the case.
	%---------------------------------------------------------------
	d = spm_FcUtil('In',oxCon,xX.xKXs,xCon);

	if spm_FcUtil('0|[]',oxCon,xX.xKXs)

	    %-Contrast was colinear with a previous one - drop it
	    %-----------------------------------------------------------
	    Ic(i)    = [];
	    i        = i - 1;

	elseif any(d)

	    %-Contrast unchanged or already defined - note index
	    %-----------------------------------------------------------
	    Ic(i)    = min(d);

	else

	    %-Define orthogonalised contrast as new contrast
	    %-----------------------------------------------------------
	    oxCon.name = [xCon(Ic(i)).name,' (orthogonalized w.r.t {',...
    	    	sprintf('%d,',Ic(1:i-2)), sprintf('%d})',Ic(i-1))];
    	    xCon  = [xCon, oxCon];
	    Ic(i) = length(xCon); 
	end

    end % while...
end % if length(Ic)...


%-Get contrasts for masking
%-----------------------------------------------------------------------
if spm_input('mask with other contrast(s)','+1','y/n',[1,0],2)
	[Im,xCon] = spm_conman(xX,xCon,'T&F',-Inf,...
		'Select contrasts for masking...',' for masking',wOK);

	%-Threshold for mask (uncorrected p-value)
	%---------------------------------------------------------------
	pm = spm_input('uncorrected mask p-value','+1','r',0.05,1,[0,1]);

	%-Inclusive or exclusive masking
	%---------------------------------------------------------------
	Ex = spm_input('nature of mask','+1','b','inclusive|exclusive',[0,1]);
else
	Im = [];
	pm = [];
	Ex = [];
end


%-Save contrast structure (if wOK) - ensures new contrasts are saved
%-----------------------------------------------------------------------
if wOK, save(fullfile(swd,'xCon.mat'),'xCon'), end


%-Create/Get title string for comparison
%-----------------------------------------------------------------------
if length(Ic)==1
	str  = xCon(Ic).name;
else
	str  = [sprintf('contrasts {%d',Ic(1)),sprintf(',%d',Ic(2:end)),'}'];
end
if Ex
	mstr = 'masked [exclusive] by';
else
	mstr = 'masked [inclusive] by';
end
if length(Im)==1
	str  = sprintf('%s (%s %s at p=%g)',str,mstr,xCon(Im).name,pm);

elseif ~isempty(Im)
	str  = [sprintf('%s (%s {%d',str,mstr,Im(1)),...
		sprintf(',%d',Im(2:end)),...
		sprintf('} at p=%g)',pm)];
end
titlestr = spm_input('title for comparison','+1','s',str);



%-Bayesian or classical Inference?
%-----------------------------------------------------------------------
if exist(fullfile(swd,'PPM.mat')) & xCon(Ic).STAT == 'T'
    if length(Ic) == 1 & isempty(xCon(Ic(1)).Vcon)
	if spm_input('Inference',1,'b',{'Bayesian','classical'},[1 0]);
		xCon(Ic).STAT = 'P';
	end
    end
end

% map parameter and hyperarameter estimates
%-----------------------------------------------------------------------
if xCon(Ic(1)).STAT == 'P'

	% Conditional estimators and error variance hyperparameters
	%----------------------------------------------------------------
	load(fullfile(swd,'PPM.mat'));
	xSDM.Vbeta  = ...
	spm_vol([repmat([swd,filesep],length(Vbeta),1),char(Vbeta)]);
	xSDM.VHp    = ...
	spm_vol([repmat([swd,filesep],length(VHp),1),char(VHp)]);
	xSDM.PPM    = PPM;

else
	% OLS estimators and error variance estimate
	%----------------------------------------------------------------
	xSDM.Vbeta  = ...
	spm_vol([repmat([swd,filesep],length(xSDM.Vbeta),1),char(xSDM.Vbeta)]);
	xSDM.VResMS = spm_vol(fullfile(swd,xSDM.VResMS));
end

%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
spm('Pointer','Watch')
spm_progress_bar('Init',100,'computing...')                          %-#
nPar   = size(xX.X,2);
I      = unique([Ic,Im]);
for ii = 1:length(I)

    i  = I(ii);

    %-Canonicalise contrast structure with required fields
    %-------------------------------------------------------------------
    if (~isfield(xCon(i),'eidf') | isempty(xCon(i).eidf))
	X1o           = spm_FcUtil('X1o',xCon(i),xX.xKXs);
	[trMV,trMVMV] = spm_SpUtil('trMV',X1o,xX.V);
        xCon(i).eidf  = trMV^2/trMVMV;
    else
        trMV = []; trMVMV = [];
    end

    %-Write contrast/ESS images?
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'Vcon') | isempty(xCon(i).Vcon) | ...
        ~exist(fullfile(swd,xCon(i).Vcon),'file')
        
        %-Bomb out (nicely) if can't write to results directory
        %---------------------------------------------------------------
        if ~wOK, spm('alert*',{	'Can''t write to the results directory:',...
		['        ',swd],' ','=> can''t compute new contrasts'},...
		mfilename,sqrt(-1));
		spm('Pointer','Arrow')
		error('can''t write contrast image')
        end


	switch(xCon(i).STAT)
	case {'T','P'} %-Implement contrast as sum of scaled beta images
        %---------------------------------------------------------------
            fprintf('\t%-32s: %-10s%20s',sprintf('contrast image %2d',i),...
                                      '(spm_add)','...initialising') %-#

	    Q     = find(abs(xCon(i).c) > 0);
	    V     = xSDM.Vbeta(Q);
	    for j = 1:length(Q)
	        V(j).pinfo(1,:) = V(j).pinfo(1,:)*xCon(i).c(Q(j));
	    end
	    
	    %-Prepare handle for contrast image
	    %-----------------------------------------------------------
	    xCon(i).Vcon = struct(...
	        'fname',  fullfile(swd,sprintf('con_%04d.img',i)),...
                'dim',    [DIM',16],...
                'mat',    M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM contrast - %d: %s',i,xCon(i).name));

            %-Write image
	    %-----------------------------------------------------------
            fprintf('%s%20s',sprintf('\b')*ones(1,20),'...computing')%-#
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);
            xCon(i).Vcon.pinfo(1,1) = spm_add(V,xCon(i).Vcon);
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);
            
            fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
                '...written %s',spm_str_manip(xCon(i).Vcon.fname,'t')))%-#

	case 'F'  %-Implement ESS as sum of squared weighted beta images
        %---------------------------------------------------------------
            fprintf('\t%-32s: %30s',sprintf('ESS image %2d',i),...
                                                     '...computing') %-#

            %-Residual (in parameter space) forming mtx
	    %-----------------------------------------------------------
            h       = spm_FcUtil('Hsqr',xCon(i),xX.xKXs);

	    %-Prepare handle for ESS image
	    %-----------------------------------------------------------
	    xCon(i).Vcon = struct(...
	        'fname',  fullfile(swd,sprintf('ess_%04d.img',i)),...
                'dim',    [DIM',16],...
                'mat',    M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM ESS - contrast %d: %s',i,xCon(i).name));

            %-Write image
	    %-----------------------------------------------------------
            fprintf('%s',sprintf('\b')*ones(1,30))                   %-#
            xCon(i).Vcon  = spm_create_image(xCon(i).Vcon);
            xCon(i).Vcon  = spm_resss(xSDM.Vbeta,xCon(i).Vcon,h);
            xCon(i).Vcon  = spm_create_image(xCon(i).Vcon);


	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(i).STAT,'"'])

	end % (switch(xCon...)

    else

	%-Already got contrast/ESS image - remap it w/ full pathname
        %---------------------------------------------------------------
	xCon(i).Vcon = spm_vol(fullfile(swd,xCon(i).Vcon));

    end % (if ~isfield...)

    spm_progress_bar('Set',100*ii/(length(I)))             	    %-#


    %-Write inference SPM/PPM image(s)
    %===================================================================
    if ~isfield(xCon(i),'Vspm') | isempty(xCon(i).Vspm) | ...
        ~exist(fullfile(swd,xCon(i).Vspm),'file')
	
        %-Bomb out (nicely) if can't write to results directory
        %---------------------------------------------------------------
        if ~wOK, spm('alert*',{	'Can''t write to the results directory:',...
		['        ',swd],' ','=> can''t compute new contrasts'},...
		mfilename,sqrt(-1));
		spm('Pointer','Arrow')
		error('can''t write SPM image')
        end
        fprintf('\t%-32s: %30s',sprintf('spm{%c} image %2d',xCon(i).STAT,i),...
                                                    '...computing')  %-#

	switch(xCon(i).STAT)
	case 'T'                                  %-Compute SPM{t} image
        %---------------------------------------------------------------
	cB  = spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	l   = spm_sample_vol(xSDM.VResMS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	VcB = xCon(i).c'*xX.Bcov*xCon(i).c; 
	Z   = cB./sqrt(l*VcB);

        str = sprintf('[%.2g]',xX.erdf);


	case 'P'                                  %-Compute PPM{P} image
        %---------------------------------------------------------------
	c     = xCon(i).c;
	cB    = spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	VcB   = c'*PPM.Cby*c;
	for j = 1:length(PPM.l)

		% hyperparameter and Taylor approximation
		%-------------------------------------------------------
		l   = spm_sample_vol(xSDM.VHp(j),XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
		VcB = VcB + (c'*PPM.dC{j}*c)*(l - PPM.l(j));
	end

	% get Bayesian threshold
        %---------------------------------------------------------------
	spm('Pointer','Arrow')
	str   = 'contrast threshold (default: prior s.d)'
 	g     = spm_input(str,'+1','e',sprintf('%0.2f',sqrt(c'*PPM.Cb*c)));

	% posterior probability cB > g
        %---------------------------------------------------------------
	Z     = 1 - spm_Ncdf(g,cB,VcB);
        str   = sprintf('[%.2f]',g);
	xCon(i).name  = [xCon(i).name ' ' str];
	xCon(i).eidf  = g;


	case 'F'                                  %-Compute SPM{F} image
        %---------------------------------------------------------------
	MVM  = spm_sample_vol(xCon(i).Vcon,XYZ(1,:),XYZ(2,:),XYZ(3,:),0)/trMV;
	RVR  = spm_sample_vol(xSDM.VResMS, XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	Z    = MVM./RVR;

        str = sprintf('[%.2g,%.2g]',xCon(i).eidf,xX.erdf);

	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(i).STAT,'"'])
	end


        %-Write full statistic image
        %---------------------------------------------------------------
        fprintf('%s%30s',sprintf('\b')*ones(1,30),'...writing')      %-#
        xCon(i).Vspm = struct(...
	    'fname',  fullfile(swd,sprintf('spm%c_%04d.img',xCon(i).STAT,i)),...
	    'dim',    [DIM',16],...
	    'mat',    M,...
	    'pinfo',  [1,0,0]',...
	    'descrip',sprintf('SPM{%c_%s} - contrast %d: %s',...
	                       xCon(i).STAT,str,i,xCon(i).name));

        tmp = zeros(DIM');
	tmp(cumprod([1,DIM(1:2)'])*XYZ -sum(cumprod(DIM(1:2)'))) = Z;

	xCon(i).Vspm = spm_write_vol(xCon(i).Vspm,tmp);

	clear tmp Z
        fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(i).Vspm.fname,'t')))  %-#

    else

	%-Already got statistic image - remap it w/ full pathname
        %---------------------------------------------------------------
	xCon(i).Vspm = spm_vol(fullfile(swd,xCon(i).Vspm));

    end % (if ~isfield...)

end % (for ii = 1:length(I))
spm_progress_bar('Clear')                                            %-#



%-Compute (unfiltered) SPM pointlist for requested masked conjunction
%=======================================================================
fprintf('\t%-32s: %30s','SPM computation','...initialising')         %-#

%-Check conjunctions - Must be same STAT w/ same df
%-----------------------------------------------------------------------
if (length(Ic) > 1) & (any(diff(double(cat(1,xCon(Ic).STAT)))) | ...
		       any(abs(diff(cat(1,xCon(Ic).eidf))) > 1) )
	error('illegal conjunction: can only conjoin SPMs of same STAT & df')
end

%-Compute mask first
%-----------------------------------------------------------------------
if ~isempty(Im)
	fprintf('%s%30s',sprintf('\b')*ones(1,30),'...masking'), end %-#

for i = Im
	Mask = spm_sample_vol(xCon(i).Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	um   = spm_u(pm,[xCon(i).eidf,xX.erdf],xCon(i).STAT);
	if Ex
		Q = Mask <= um;
	else
		Q = Mask >  um;
	end
	XYZ  = XYZ(:,Q);
	QQ   = QQ(Q);
	if isempty(Q), break, end
end
if ~isempty(Im) & ~isempty(Q)
	Msk  = XYZ(1,:)+(XYZ(2,:)-1)*DIM(1)+(XYZ(3,:)-1)*DIM(1)*DIM(2);
else
	Msk  = 0;  %-Statistic images are zero-masked
end


%-Compute conjunction SPM (as minimum SPM) within mask...
%-----------------------------------------------------------------------
if isempty(XYZ)
    fprintf('\n')                                                    %-#
    warning(sprintf('No voxels survive masking at p=%4.3g',pm))
    Z     = [];
else
    fprintf('%s%30s',sprintf('\b')*ones(1,30),'...computing')        %-#
    Z     = Inf;
    for i = Ic
        Z = min(Z,spm_sample_vol(xCon(i).Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0));
    end
end

%-Degrees of Fredom and STAT string describing marginal distribution
%-----------------------------------------------------------------------
edf   = [xCon(Ic(1)).eidf xX.erdf];
if length(Ic) > 1
	str = sprintf('^{%d}',length(Ic));
else
	str = '';
end

switch xCon(Ic(1)).STAT
case 'T'
	STATstr = sprintf('%c%s_{%.1f}','T',str,edf(2));
case 'F'
	STATstr = sprintf('%c%s_{%.1f,%.1f}','F',str,edf(1),edf(2));
case 'P'
	STATstr = sprintf('%s^{%0.2f}','PPM',edf(1));
end


fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#

%-Save mappped statistic image(s) to put in VOL
%-----------------------------------------------------------------------
VspmSv = cat(1,xCon(Ic).Vspm);

%-Save contrast structure (if wOK), with relative pathnames to image files
%=======================================================================
if wOK
    for i = I;
	if ~isempty(xCon(i).Vcon)
        	xCon(i).Vcon = spm_str_manip(xCon(i).Vcon.fname,'t');
	end
	if ~isempty(xCon(i).Vspm)
        	xCon(i).Vspm = spm_str_manip(xCon(i).Vspm.fname,'t');
	end
    end
    save(fullfile(swd,'xCon.mat'),'xCon')
    fprintf('\t%-32s: %30s\n','contrast structure','...saved to xCon.mat')%-#
end

%-Various parameters...
%-----------------------------------------------------------------------
STAT = xCon(Ic(1)).STAT;
n    = length(Ic);
u    = -Inf;
k    = 0;

%-Return to previous directory & clean up interface
%-----------------------------------------------------------------------
spm('Pointer','Arrow')


%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================

%-Height threshold - classical inference
%-----------------------------------------------------------------------
if STAT ~= 'P'


    %-Get height threshold
    %-------------------------------------------------------------------
    str   = 'FWE|FDR|none';
    d = spm_input('p value adjustment to control','+1','b',str,[],3);

    switch d

	case 'FWE' % family-wise false positive rate
        %---------------------------------------------------------------
	u  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
	u  = spm_uc(u,edf,STAT,xSDM.R,n,xSDM.S);

	case 'FDR' % False discovery rate
	%---------------------------------------------------------------	
	u  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
	u  = spm_uc_FDR(u,edf,STAT,n,VspmSv,0);

	otherwise  %-NB: no adjustment
	% p for conjunctions is p of the conjunction SPM
   	%---------------------------------------------------------------
	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1);
	if u <= 1; u = spm_u(u^(1/n),edf,STAT); end

    end

%-Height threshold - Bayesian inference
%-----------------------------------------------------------------------
elseif STAT == 'P'

	u  = spm_input(['p value threshold for PPM'],'+0','r',.95,1);

end % (if STAT)

%-Calculate height threshold filtering
%-------------------------------------------------------------------
Q      = find(Z > u);

%-Apply height threshold
%-------------------------------------------------------------------
Z      = Z(:,Q);
XYZ    = XYZ(:,Q);
QQ     = QQ(Q);
if isempty(Q)
	warning(sprintf('No voxels survive height threshold u=%0.2g',u)),
end


%-Extent threshold (disallowed for conjunctions)
%-----------------------------------------------------------------------
if ~isempty(XYZ) & length(Ic) == 1

    %-Get extent threshold [default = 0]
    %-------------------------------------------------------------------
    k     = spm_input('& extent threshold {voxels}','+1','r',0,1,[0,Inf]);

    %-Calculate extent threshold filtering
    %-------------------------------------------------------------------
    A     = spm_clusters(XYZ);
    Q     = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= k; Q = [Q j]; end
    end

    % ...eliminate voxels
    %-------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    QQ    = QQ(Q);

    if isempty(Q)
	warning(sprintf('No voxels survive extent threshold k=%0.2g',k)), end

else

    k = 0;

end % (if ~isempty(XYZ))


%=======================================================================
% - E N D
%=======================================================================

%-Assemble output structures of unfiltered data
%=======================================================================
SPM    = struct('swd',		swd,...
		'title',	titlestr,...
		'Z',		Z,...
		'n',		n,...
		'STAT',		STAT,...
		'df',		edf,...
		'STATstr',	STATstr,...
		'Ic',		Ic,...
		'Im',		Im,...
		'pm',		pm,...
		'Ex',		Ex,...
		'u',		u,...
		'k',		k,...
		'XYZ',		XYZ,...
		'XYZmm',	M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
		'QQ',		QQ);

VOL    = struct('S',		S,...
		'R',		R,...
		'FWHM',		xSDM.FWHM,...
		'M',		M,...
		'iM',		iM,...
		'VOX',		sqrt(sum(M(1:3,1:3).^2))',...
		'DIM',		DIM,...
		'Vspm',		VspmSv,...
		'Msk',		Msk);

% RESELS per voxel (density) if it exists
%-----------------------------------------------------------------------
if isfield(xSDM,'VRVP'), VOL.VRVP = spm_vol(xSDM.VRVP); end

