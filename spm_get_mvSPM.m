function [SPM,VOL,xX,xCon,xSDM] = spm_get_mvSPM(swd,xSDM)
% second-level mulitvariate inference following parameter estimation
% FORMAT [SPM,VOL,xX,xCon,xSDM] = spm_get_mvSPM(swd,xSDM);
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
% .X0    - Reduced design matrix (spans design space under Ho)
% .iX0   - Indicies of design matrix columns to form the reduced design matrix.
% .X1o   - Remaining design space (orthogonal to X0).
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
%___________________________________________________________________________
% %W% Andrew Holmes, Karl Friston & Jean-Baptiste Poline %E%
SCCSid   = '2.23';

%-GUI setup
%---------------------------------------------------------------------------
SPMid  = spm('SFnBanner',mfilename,SCCSid);
spm_help('!ContextHelp',mfilename)


%-Get Stats data from SPM.mat
%---------------------------------------------------------------------------
xX     = xSDM.xX;			%-Design definition structure
XYZ    = xSDM.XYZ;			%-XYZ coordinates
S      = xSDM.S;			%-search Volume {voxels}
R      = xSDM.R;			%-search Volume {resels}
Sess   = xSDM.Sess;			% Session structure
xSDM   = rmfield(xSDM,{'xX','XYZ'});	%-Remove (large) duplicate fields

%-Get/Compute mm<->voxel matrices & image dimensions from SPM.mat
%---------------------------------------------------------------------------
M      = xSDM.M;
iM     = inv(M);
DIM    = xSDM.DIM;

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



%-Build index from XYZ into corresponding Y.mad locations
%---------------------------------------------------------------------------
QQ     = zeros(1,S);
if exist(fullfile(swd,'Yidx.mat'),'file') & exist(fullfile(swd,'Y.mad'),'file')
	load(fullfile(swd,'Yidx.mat'))
	QQ(Yidx) = 1:length(Yidx);
end

%-restrict analysis to 1st level p < UFp
%-----------------------------------------------------------------------
Q      = find(QQ);
XYZ    = XYZ(:,Q);
QQ     = QQ(Q);


%-2nd-level specification
%===========================================================================


%-get parameter estimates (now the response variable at this level)
%---------------------------------------------------------------------------
Ts     = spm_input('trials to test','+1','m',Sess{1}.name);
q      = length(Sess);
p      = length(Sess{1}.ind{Ts});
str    = sprintf('[1:%i]',p);
Cs     = spm_input('coefficents to enter','+1','n',str,[1 Inf],p);
for  i = 1:q
	I{i} = Sess{i}.col(Sess{i}.ind{Ts}(Cs));
end
p      = length(Cs);					% p-variate repsonse

%-get second level design matrix
%---------------------------------------------------------------------------
str    = sprintf('ones(%i,1)',q);
X      = spm_input('2nd-level design matrix','+1','n',str,[q Inf]);

if spm_input('any confounds?','+1','y/n',[1 0]);

	G  = spm_input('confounds','+1','r',str,[q Inf]);
	X  = X - G*(pinv(G)*X);
else
	G  = [];
end




%-compute Chi-squared feild SPM{X}
%===========================================================================

% degrees of freedom
%---------------------------------------------------------------------------
h     = rank(X);					% condition
r     = q - h - rank(G);				% residuals
erdf  = p*h;						% Chi-squared
Z     = zeros(1,length(XYZ));
iX    = pinv(X);
iG    = pinv(G);

spm_progress_bar('Init',100,'computing...')                          %-#

for i = 1:length(XYZ)

	% re-organize parameter estimates into n-variate variable response Y
	%-------------------------------------------------------------------
	Y     = zeros(q,p);
	for j = 1:q
		for k = 1:p
			Y(j,k) = spm_sample_vol(xSDM.Vbeta(I{j}(k)), ...
				 XYZ(1,i),XYZ(2,i),XYZ(3,i),0);	
		end
	end

	% Remove confounds
	%-------------------------------------------------------------------
	if length(G)
		Y = Y - G*(iG*Y);
	end

	% ManCova
	%-------------------------------------------------------------------
	Res   = Y - X*(iX*Y);				% residuals

	% test for the alternative hypothesis
	%-------------------------------------------------------------------
	Z(i)  = (r - ((p - h + 1)/2))*log(det(Y'*Y)/det(Res'*Res));

	% progress
	%-------------------------------------------------------------------
	if ~rem(i,100)
		spm_progress_bar('Set',100*i/length(XYZ))
	end
end



%-Create/Get title string for comparison
%-----------------------------------------------------------------------
titlestr  = sprintf('2nd-level %i-variate inference',p);

%-Generate STAT string describing marginal distribution
%-----------------------------------------------------------------------
STATstr   = sprintf('%c%s_{%i}','X','',erdf);


%-Various parameters...
%-----------------------------------------------------------------------
Ic   = [];
Im   = [];
pm   = [];
Ex   = [];
n    = 1;
STAT = 'X';
edf  = [h, erdf];
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
