function [SPM,VOL,xX,xCon,xSDM] = spm_getSPM
% computs a specified and thresholded following parameter estimation
% FORMAT [SPM,VOL,xX,xCon,xSDM] = spm_getSPM;
%
% SPM    - structure containing SPM, distribution & filtering details
% .swd   - SPM working directory - directory containing current SPM.mat
% .title - title for comparison (string)
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
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
% 
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xCon   - Contrast definitions structure array
%        - (see also spm_FcUtil.m for structure, rules & handling)
% .name  - Contrast name
% .STAT  - Statistic indicator character ('T' or 'F')
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
% .Vcon  - Name of contrast (for 'T's) or ESS (for 'F's) image
% .Vspm  - Name of SPM image
%
% xSDM   - structure containing contents of SPM.mat file
%          ( see spm_spm.m for contents...
%          ( ...except that xX & XYZ are removed, being returned separately
%          ( ...and xSDM.Vbeta & xSDM.VResMS are re-mmapped
%
%_______________________________________________________________________
%
% **** files written
% **** reference to contrast manager
% **** orthogonalisation order
%
%
% spm_getSPM prompts for an SPM and applies thresholds {u & k}
% to a point list of voxel values (specified with their locations {XYZ})
% This allows the SPM be displayed and characterized in terms of regionally 
% significant effects by subsequent routines.
% 
% Either single contrasts can be examined or conjunctions of different
% contrasts.  In the latter case a new SPM is created, that reflects
% the miminum of all specified effects.  A conjunction is therefore
% the conjoint expression of two or more effects.  The contrasts are
% successively enforced to be orthogonal so the order of non-orthogonal
% contrast specification is important.
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
% see spm_results_ui.m for further details
%
%_______________________________________________________________________
% %W% Andrew Holmes, Karl Friston & Jean-Baptiste Poline %E%
SCCSid   = '%I%';

%-GUI setup
%-----------------------------------------------------------------------
SPMid  = spm('SFnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',mfilename)

%-Select SPM.mat & go to that directory
%-----------------------------------------------------------------------
swd    = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat'),'H');

%-Preliminaries...
%=======================================================================

%-Load and check/canonicalise SPM.mat
%-----------------------------------------------------------------------
xSDM   = load(fullfile(swd,'SPM.mat'));	%-Contents of SPM.mat (in a struct)

if ~isfield(xSDM,'SPMid')
	%-SPM.mat pre SPM99
 	error('Incompatible SPM.mat - old SPM results format!?')
elseif ~isfield(xSDM,'M')
	%-SPM.mat from SPM99b (which saved mmapped handles) **
	xSDM.M      = xSDM.Vbeta(1).mat;
	xSDM.DIM    = xSDM.Vbeta(1).dim(1:3)';
	xSDM.VM     = 'mask.img';
	xSDM.Vbeta  = {xSDM.Vbeta.fname}';
	xSDM.VResMS = xSDM.VResMS.fname;
end


%-Map beta and ResMS images, canonicalising relative pathnames
% (SPM result images stored with relative pathnames, for robustness.)
%-----------------------------------------------------------------------
xSDM.Vbeta  = ...
	spm_vol([repmat([swd,filesep],length(xSDM.Vbeta),1),char(xSDM.Vbeta)]);
xSDM.VResMS = spm_vol(fullfile(swd,xSDM.VResMS));

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
QQ     = zeros(1,S);
if exist(fullfile(swd,'Yidx.mat'),'file') & exist(fullfile(swd,'Y.mad'),'file')
	load(fullfile(swd,'Yidx.mat'))
	QQ(Yidx) = 1:length(Yidx);
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

%-Canonicalise SPM99b format xCon (which saved mmapped handles) **
%-----------------------------------------------------------------------
for i=1:length(xCon)
	if isstruct(xCon(i).Vcon), xCon(i).Vcon=xCon(i).Vcon.fname; end
	if isstruct(xCon(i).Vspm), xCon(i).Vspm=xCon(i).Vspm.fname; end
end

%-See if can write to current directory (by trying to resave xCon.mat)
%-----------------------------------------------------------------------
wOK = 1;
try
	save(fullfile(swd,'xCon.mat'),'xCon')
catch
	wOK = 0;
	str = {	'Can''t write to the results directory:',...
        	'(problem saving xCon.mat)',...
		['        ',swd],...
		' ','-> results restricted to contrasts already computed'};
	spm('alert!',str,mfilename,1);
end


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================

%-Get contrasts (if multivariate there is only one structure)
%-----------------------------------------------------------------------
nVar    = size(xSDM.VY,2);
if nVar == 1
	[Ic,xCon] = spm_conman(xX,xCon,'T|F',Inf,...
	'	Select contrasts...',' for conjunction',wOK);
else
	Ic = 1;
end

%-Enforce orthogonality of multiple contrasts for conjunction
% (Orthogonality within subspace spanned by contrasts)
%-----------------------------------------------------------------------
% (Don't want to ask orthogonalisation order if not needed.)
if length(Ic) > 1 & ~spm_FcUtil('|_?',xCon(Ic), xX.xKXs)

    %-Get orthogonalisation order from user
    Ic = spm_input('orthogonlization order','+1','p',Ic,Ic);
    
    %-Successively orthogonalise
    %-------------------------------------------------------------------
    i = 1; while(i < length(Ic)), i = i + 1;
	%-NB: This loop is peculiarly controlled to account for the
	%     possibility that Ic may shrink if some contrasts diasppear
	%     on orthogonalisation (i.e. if there are colinearities)
    
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
	str  = [sprintf('contrasts {%d',Ic(1)),...
		sprintf(',%d',Ic(2:end)),'}'];
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


%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
spm('Pointer','Watch')
spm_progress_bar('Init',100,'computing...')                          %-#
nPar   = size(xX.X,2);
I      = [Ic,Im];
for ii = 1:length(I)

    i  = I(ii);

    %-Canonicalise contrast structure with required fields
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'eidf') | isempty(xCon(i).eidf)
	[trMV,trMVMV] = spm_SpUtil('trMV',...
				spm_FcUtil('X1o',xCon(i),xX.xKXs),xX.V);
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

	case 'T'       %-Implement contrast as sum of scaled beta images
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

    spm_progress_bar('Set',100*(2*ii-1)/(2*length([Ic,Im])+2))       %-#

    %-Write statistic image(s)
    %-------------------------------------------------------------------
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
	Z   = spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0)./...
	(sqrt(spm_sample_vol(xSDM.VResMS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0)*...
	                                (xCon(i).c'*xX.Bcov*xCon(i).c) ));
        str = sprintf('[%.2g]',xX.erdf);

	case 'F'                                  %-Compute SPM{F} image
        %---------------------------------------------------------------
	if isempty(trMV)
	    trMV = spm_SpUtil('trMV',spm_FcUtil('X1o',xCon(i),xX.xKXs),xX.V);
        end
	Z =(spm_sample_vol(xCon(i).Vcon,XYZ(1,:),XYZ(2,:),XYZ(3,:),0)/trMV)./...
	   (spm_sample_vol(xSDM.VResMS, XYZ(1,:),XYZ(2,:),XYZ(3,:),0));

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

	xCon(i).Vspm       = spm_write_vol(xCon(i).Vspm,tmp);

	clear tmp Z
        fprintf('%s%30s\n',sprintf('\b')*ones(1,30),sprintf(...
            '...written %s',spm_str_manip(xCon(i).Vspm.fname,'t')))  %-#

    else

	%-Already got statistic image - remap it w/ full pathname
        %---------------------------------------------------------------
	xCon(i).Vspm = spm_vol(fullfile(swd,xCon(i).Vspm));

    end % (if ~isfield...)

    spm_progress_bar('Set',100*(2*ii-0)/(2*length([Ic,Im])+2))       %-#

end % (for ii = 1:length(I))



%-Compute (unfiltered) SPM pointlist for requested masked conjunction
%=======================================================================
fprintf('\t%-32s: %30s','SPM computation','...initialising')         %-#

%-Check conjunctions - Must be same STAT w/ same df
%-----------------------------------------------------------------------
if (length(Ic) > 1) & (any(diff(double(cat(1,xCon(Ic).STAT)))) | ...
		       any(abs(diff(cat(1,xCon(Ic).eidf))) > 1e-6) )
	error('illegal conjunction: can only conjoin SPMs of same STAT & df')
end

%-Compute mask first
%-----------------------------------------------------------------------
if ~isempty(Im), fprintf('%s%30s',sprintf('\b')*ones(1,30),'...masking'), end %-#
for i = Im

	Mask  = spm_sample_vol(xCon(i).Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
	um    = spm_u(pm,[xCon(i).eidf,xX.erdf],xCon(i).STAT);
	if Ex
		Q     =  Mask <= um;
	else
		Q     =  Mask >  um;
	end
	XYZ   = XYZ(:,Q);
	QQ    = QQ(Q);
	if isempty(Q), break, end
end

spm_progress_bar('Set',100*((2*length([Ic,Im])+1)/(2*length([Ic,Im])+2)))%-#

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
elseif nVar > 1
	str = sprintf('^{%d-variate}',nVar);
else
	str = '';
end

switch xCon(Ic(1)).STAT
case 'T'
	STATstr = sprintf('%c%s_{%.4g}','T',str,edf(2));
case 'F'
	STATstr = sprintf('%c%s_{[%.4g,%.4g]}','F',str,edf(1),edf(2));
end


fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')               %-#
spm_progress_bar('Set',100)                                          %-#


%-Save contrast structure (if wOK), with relative pathnames to image files
%=======================================================================
if wOK
    for i = [Ic,Im];
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
spm_progress_bar('Clear')                                            %-#
spm('Pointer','Arrow')


%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================

%-Height threshold
%-----------------------------------------------------------------------
if ~isempty(XYZ)

    %-Get height threshold
    %-------------------------------------------------------------------
    if spm_input('corrected height threshold','+1','y/n',[1,0],2)
	u  = spm_input('corrected p value','+0','r',0.05,1,[0,1]);
	u  = spm_uc(u,edf,STAT,xSDM.R,n);
    else
	%-NB: Uncorrected p for conjunctions is p of the conjunction SPM
	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1);
	if u <= 1; u = spm_u(u^(1/n),edf,STAT); end
    end

    %-Calculate height threshold filtering
    %-------------------------------------------------------------------
    Q     = find(Z > u);

    %-Apply height threshold
    %-------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    QQ    = QQ(Q);

    if isempty(Q)
        warning(sprintf('No voxels survive height threshold u=%0.2g',u)), end

end % (if ~isempty(XYZ))


%-Extent threshold (only for allowed cases)
%-----------------------------------------------------------------------
if ~isempty(XYZ) & length(Ic) == 1 & STAT == 'T'

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
		'DIM',		DIM);

% RESELS per voxel (density) if it exists
%-----------------------------------------------------------------------
if isfield(xSDM,'VRVP'), VOL.VRVP = spm_vol(xSDM.VRVP); end

