function [SPM,VOL,xX,xSDM] = spm_getSPM
% computs a specified and thresholded following parameter estimation
% FORMAT [SPM,VOL,xX,xSDM] = spm_getSPM;
%
% SPM    - structure containing SPM, distribution & filtering detals
% .swd   - SPM working directory - directory containing current SPM.mat
% .title - title for comparison (string)
% .c     - contrast(s) - in cell array
% .Z     - minimum of n Statistics {filtered on u and k}
% .n     - number of conjoint tests        
% .STAT  - distribution {Z, T, X or F}     
% .df    - degrees of freedom [df{interest}, df{residual}]
% .u     - height threshold
% .k     - extent threshold {resels}
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
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
% 
%
% xX     - Design Matrix structure
%        - (see spm_spm.m for structure)
%
% xSDM   - structure containing contents of SPM.mat file
%          ( see spm_spm.m for contents...
%          ( ...except that xX & XYZ are removed, being returned separately
%
%_______________________________________________________________________
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
% number of resels comprising its clusters {k}. The height threshold is
% specified as above in terms of an [un]corrected p value or
% statistic.  Clusters can also be thresholded on the basis of their
% spatial extent.  (SPM99 expresses all cluster volumes in terms of
% resels.  These are effective resolution elements and comprise several
% voxels, depending on the smoothness estimator.) If you want to see
% all voxels simply enter 0.  In this instance the 'set-level'
% inference can be considered an 'omnibus test' based on the number of
% clusters that obtain.
% 
% see spm_list.m and spm_P.m for further details
%
%_______________________________________________________________________
% %W% Karl Friston, Andrew Holmes %E%


%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

%-Select SPM.mat
%-----------------------------------------------------------------------
swd  = spm_str_manip(spm_get(1,'SPM.mat','Select SPM.mat for analysis'),'H');


%-Preliminaries...
%=======================================================================

%-Get Stats data from SPM.mat
%-----------------------------------------------------------------------
xSDM = load(fullfile(swd,'SPM.mat'));	%-Contents of SPM.mat (in a struct)
xX   = xSDM.xX;				%-Design definition structure
XYZ  = xSDM.XYZ;			%-XYZ coordinates
S    = xSDM.S;				%-search Volume {voxels}
R    = xSDM.R;				%-search Volume {resels}
xSDM = rmfield(xSDM,{'xX','XYZ'});	%-Remove (large) duplicate fields

%-Compute mm<->voxel matrices
%-----------------------------------------------------------------------
 M  = xSDM.Vbeta(1).mat;
iM  = inv(M);
dim = xSDM.Vbeta(1).dim(1:3);

%-Build index from XYZ into corresponding Y.mad locations
%-----------------------------------------------------------------------
tmp = fullfile(swd,'Yidx.mat');
if exist(tmp,'file'), load(tmp), else, Yidx = []; end
QQ         = zeros(1,S);
QQ(Yidx)   = 1:length(Yidx);


%-Canonicalise relative pathnames in xSDM by prepending swd
% (SPM result images kept with relative pathnames so SPM is robust to
% (user moving files around, but use full pathnames here so code is robust
% (to changing directory whilst using the results section
%-----------------------------------------------------------------------
nbeta = length(xSDM.Vbeta);
for i=1:nbeta, xSDM.Vbeta(i).fname = fullfile(swd,xSDM.Vbeta(i).fname); end
xSDM.VResMS.fname = fullfile(swd,xSDM.VResMS.fname);


%-Load contrast definitions (if available)
%-----------------------------------------------------------------------
tmp = fullfile(swd,'xCon.mat');
if exist(tmp,'file'), load(tmp), else, xCon = []; end


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================

%-Get contrasts...
%-----------------------------------------------------------------------
[Ic,xCon] = spm_conman(xX,xCon,'T|F',Inf,...
	'Select contrasts...',' for conjunction',1);


%-Enforce orthogonality of multiple contrasts for conjunction
% (Orthogonality within subspace spanned by contrasts)
%-----------------------------------------------------------------------
tol = 1e-10;					%-Tolerance

Xc = xX.xKXs.X * cat(2,xCon(Ic).c);		%-NB: This amalgamates columns of
						% F-contrasts, perhaps causing
						% orthogonalisation when only an
						% F-contrast has non-ortho. cols

if length(Ic)>1 & any(any(triu(Xc'*Xc,1)>tol))	%-(Probably) need 2 orthogonalise

    %-Get orthogonalisation order from user
    Ic = spm_input('orthogonlization order','+1','p',Ic,Ic);
    
    %-Successively orthogonalise...
    i=1; while(i<length(Ic)), i = i +1;
	%-NB: This loop is peculiarly controlled to account for the
	%     possibility that Ic may shrink if some contrasts diasppear
	%     on orthogonalisation (i.e. if there are colineaities)
    
	%-Orthogonalise (subspace spanned by) contrast i wirit previous
	%---------------------------------------------------------------
	Xc = xX.xKXs.X*cat(2,xCon(Ic(i    )).c);
	XC = xX.xKXs.X*cat(2,xCon(Ic(1:i-1)).c);
	Xc = Xc - XC*pinv(XC)*Xc;
	c  = xX.pKX*Xc;
    
	%-See if this orthogonalised contrast has already been entered
	% or is colinear with a previous one. Define a new contrast if
	% neither is the case.
	%---------------------------------------------------------------
	d  = ones(1,length(xCon));
	for j = 1:length(xCon)
	    if xCon(Ic(i)).STAT==xCon(j).STAT & size(c,2)==size(xCon(j).c,2)
		d(j) = max(max(abs(c-xCon(j).c)));
	    end
	end 
	if any(max(abs(c))<tol)
	    %-A contrast column was colinear with a previous one
	    % (NB: Might be a single constraint of an F-contrast)
	    c(:,max(abs(c))<tol) = [];
	end
	if isempty(c)
	    %-Contrast was completely colinear with a previous one - drop it
	    Ic(i)  = [];
	    i      = i -1;
	elseif any(d<tol)
	    %-Contrast unchanged or already defined - note index
	    Ic(i)  = min(find(d<tol));
	else
	    %-Define orthogonalised contrast as new contrast
	    j = length(xCon)+1;
    
	    xCon(j).name = [xCon(Ic(i)).name,'  (orthogonalized wirit {',...
		    sprintf('%d,', Ic(1:i-2)),sprintf('%d})',Ic(i-1  ))];
	    xCon(j).c    = c;
	    xCon(j).STAT = xCon(Ic(i)).STAT;
	    Ic(i)        = j;
	end
    end % while...
end % if length(Ic)...


%-Get contrasts for masking
%-----------------------------------------------------------------------
if spm_input('mask with other contrast(s)','+1','y/n',[1,0],2)
	[Im,xCon] = spm_conman(xX,xCon,'T&F',-Inf,...
		'Select contrasts for masking...',' for masking',1);

	%-Threshold for mask (uncorrected p-value)
	pm = spm_input('uncorrected mask p-value','+1','r',0.05,1,[0,1]);
else
	Im = [];
	pm = [];
end


%-Preliminary save of the contrast structure (to see if we can write)
%-----------------------------------------------------------------------
try
	save(fullfile(swd,'xCon.mat'),'xCon')
catch
	str = {	'Problem saving xCon.mat:','',...
		'SPM needs write access to results directory!'};
	msgbox(str,sprintf('%s%s: %s...',spm('ver'),...
		spm('GetUser',' (%s)'),mfilename),'warn','modal')
	error(['can''t write to results directory: ',swd])
end

%-Get title for comparison
%-----------------------------------------------------------------------
if length(Ic)==1
	str = xCon(Ic).name;
else
	str = [	sprintf('contrasts {%d',Ic(1)),...
		sprintf(',%d',Ic(2:end)),'}'];
end
if length(Im)==1
	str = sprintf('%s (masked by %s)',str,xCon(Im).name);
elseif ~isempty(Im)
	str = [	sprintf('%s (masked by {%d',str,Im(1)),...
		sprintf(',%d',Im(2:end)),'})'];
end
title = spm_input('title for comparison','+1','s',str);


spm('Pointer','Watch')


%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
nPar = size(xX.X,2);
I    = [Ic,Im];

for i=I

    %-Canonicalise contrast structure with required fields
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'X1o') | isempty(xCon(i).X1o)
	xCon(i).KX1o = spm_SpUtil('cTestSp',xX.xKXs,xCon(i).c);
    end
    if ~isfield(xCon(i),'trMV') | isempty(xCon(i).trMV)
	[xCon(i).trMV,xCon(i).trMVMV] = ...
	    spm_SpUtil('trMV',xCon(i).KX1o,xX.V);
    end
    if ~isfield(xCon(i),'eidf') | isempty(xCon(i).eidf)
        xCon(i).eidf = xCon(i).trMV^2 / xCon(i).trMVMV;
    end


    %-Write contrast/ESS images?
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'Vcon') | isempty(xCon(i).Vcon) | ...
        ~exist(fullfile(swd,xCon(i).Vcon.fname),'file')
	
	switch(xCon(i).STAT)
	case 'T'       %-Implement contrast as sum of scaled beta images
        %---------------------------------------------------------------
	    Q = find(abs(xCon(i).c)>0);
	    V = xSDM.Vbeta(Q);
	    for j=1:length(Q)
	        V(j).pinfo(1,:)=V(j).pinfo(1,:)*xCon(i).c(Q(j));
	    end
	    
	    %-Prepare handle for contrast image
	    xCon(i).Vcon = struct(...
	        'fname',  fullfile(swd,sprintf('con_%04d.img',i)),...
                'dim',    [dim,16],...
                'mat',    M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM contrast image - %d: %s',i,xCon(i).name));

            %-Write image
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);
            xCon(i).Vcon.pinfo(1,1) = spm_add(V,xCon(i).Vcon);
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);

	case 'F'  %-Implement ESS as sum of squared weighted beta images
        %---------------------------------------------------------------
            %-Residual (in parameter space) forming mtx
            [u,s,v] = svd(spm_SpUtil('BetaRC',xX.xKXs,xCon(i).c));
            RB      = u*sqrt(s)*v';

	    %-Prepare handle for ESS image
	    xCon(i).Vcon = struct(...
	        'fname',  fullfile(swd,sprintf('ess_%04d.img',i)),...
                'dim',    [dim,16],...
                'mat',    M,...
                'pinfo',  [1,0,0]',...
                'descrip',sprintf('SPM ESS - contrast %d: %s',i,xCon(i).name));

            %-Write image
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);
            xCon(i).Vcon = spm_resss(xSDM.Vbeta,xCon(i).Vcon,RB);
            xCon(i).Vcon            = spm_create_image(xCon(i).Vcon);

	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(i).STAT,'"'])
	end

    else
        %-Canonicalise relative pathnames of contrast/ESS image
        xCon(i).Vcon.fname = fullfile(swd,xCon(i).Vcon.fname);
    end % (if ~isfield...)


    %-Write statistic image(s)
    %-------------------------------------------------------------------
    if ~isfield(xCon(i),'Vspm') | isempty(xCon(i).Vspm) | ...
        ~exist(fullfile(swd,xCon(i).Vspm.fname),'file')
	
	switch(xCon(i).STAT)
	case 'T'                                  %-Compute SPM{t} image
        %---------------------------------------------------------------
	Z =     spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0) ./...
	  (sqrt(spm_sample_vol(xSDM.VResMS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0) * ...
	                                (xCon(i).c'*xX.Bcov*xCon(i).c)    ) );

	case 'F'                                  %-Compute SPM{F} image
        %---------------------------------------------------------------
	Z = (spm_sample_vol(xCon(i).Vcon, XYZ(1,:),XYZ(2,:),XYZ(3,:),0)...
	                                / xCon(i).trMV                 ) ./ ...
	    (spm_sample_vol(xSDM.VResMS,  XYZ(1,:),XYZ(2,:),XYZ(3,:),0));

	otherwise
        %---------------------------------------------------------------
	    error(['unknown STAT "',xCon(i).STAT,'"'])
	end

        %-Write full statistic image
        %---------------------------------------------------------------
        xCon(i).Vspm = struct(...
	    'fname',  fullfile(swd,sprintf('spm_%04d_%c.img',i,xCon(i).STAT)),...
	    'dim',    [dim,16],...
	    'mat',    M,...
	    'pinfo',  [1,0,0]',...
	    'descrip',sprintf('SPM{%c} - contrast %d: %s',...
	                                        xCon(i).STAT,i,xCon(i).name));

        tmp = zeros(dim);
	tmp(cumprod([1,dim(1:2)])*XYZ -sum(cumprod(dim(1:2)))) = Z;

	xCon(i).Vspm       = spm_write_vol(xCon(i).Vspm,tmp);
        xCon(i).Vspm.fname = spm_str_manip(xCon(i).Vspm.fname,'t');

	clear tmp Z

    end % (if ~isfield...)

    xCon(i).Vcon.fname = spm_str_manip(xCon(i).Vcon.fname,'t');

end % (for i=I)


%-Save contrast structure
%=======================================================================
save(fullfile(swd,'xCon.mat'),'xCon')


%-Compute (unfiltered) SPM pointlist for requested masked conjunction
%=======================================================================

%-Check conjunctions - Must be same STAT w/ same df
%-----------------------------------------------------------------------
if (length(Ic)>1) & (any(diff(double(cat(1,xCon(Ic).STAT)))) | ...
					 any(diff(cat(1,xCon(Ic).eidf))>tol) )
	error('illegal conjunction: can only conjoin SPMs of same STAT & df')
end

%-Compute mask first...
%-----------------------------------------------------------------------
for i=Im
	V = xCon(i).Vspm; V.fname=fullfile(swd,V.fname);
	Q = spm_sample_vol(V,XYZ(1,:),XYZ(2,:),XYZ(3,:),0) > ...
		spm_u(pm,[xCon(i).eidf,xX.erdf],xCon(i).STAT);
	XYZ   = XYZ(:,Q);
	QQ    = QQ(Q);
	if isempty(Q), break, end
end


%-Compute conjunction SPM (as minimum SPM) within mask...
%-----------------------------------------------------------------------
if isempty(XYZ)
	warning(sprintf('No voxels survive masking at p=%4.3g',pm))
	Z = [];
else
	Z = Inf;
	for i=Ic
		V = xCon(i).Vspm; V.fname=fullfile(swd,V.fname);
		Z = min(Z,spm_sample_vol(V,XYZ(1,:),XYZ(2,:),XYZ(3,:),0));	
	end
end

%-Various parameters...
%-----------------------------------------------------------------------
n    = length(Ic);
STAT = xCon(Ic(1)).STAT;
edf  = [xCon(Ic(1)).eidf, xX.erdf];
u    = -Inf;
k    = 0;

%-Unfiltered data
% uSPM   = struct('Z',		Z,...
% 		'XYZ',		XYZ,...
% 		'QQ',		QQ);

%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================
spm('Pointer','Arrow')

if ~isempty(XYZ)                                      %-Height threshold
%-----------------------------------------------------------------------
    %-Get height threshold
    %-------------------------------------------------------------------
    if spm_input('corrected height threshold','+1','y/n',[1,0],2)
	u  = spm_input('corrected p value','+0','r',0.05,1,[0,1]);
	u  = spm_U(u,edf,STAT,xSDM.R,n);
    else
	%-NB: Uncorrected p for conjunctions is p of each component comparison
	u  = spm_input(['threshold {',STAT,' or p value}'],'+0','r',0.001,1);
	if u <= 1; u = spm_u(u,edf,STAT); end
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


if ~isempty(XYZ)                                      %-Height threshold
%-----------------------------------------------------------------------
    %-Get extent threshold [default = 0]
    %-------------------------------------------------------------------
    k     = spm_input('& extent threshold {resels}','+1','r',0,1,[0,Inf]);

    %-Calculate extent threshold filtering
    %-------------------------------------------------------------------
    A     = spm_clusters(XYZ);
    Q     = [];
    for i = 1:max(A)
        j = find(A==i);
        if length(j)*R(end)/S >= k; Q = [Q j]; end
    end

    % ...eliminate voxels
    %-------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    QQ    = QQ(Q);

    if isempty(Q)
	warning(sprintf('No voxels survive extent threshold k=%0.2g',k)), end

end % (if ~isempty(XYZ))


%=======================================================================
% - E N D
%=======================================================================

%-Assemble output structures of unfiltered data
%=======================================================================
SPM    = struct('swd',		swd,...
		'title',	title,...
		'c',		{{xCon(Ic).c}},...
		'Z',		Z,...
		'n',		n,...
		'STAT',		STAT,...
		'df',		edf,...
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
		'DIM',		xSDM.Vbeta(1).dim(1:3)');
