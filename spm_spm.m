function spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,TH,Dnames,Fnames,SIGMA,RT, VMask)
% Statistical analysis with the General linear model
% FORMAT spm_spm(V,H,C,B,G,CONTRAST,ORIGIN,TH,Dnames,Fnames,SIGMA,RT);
%
% V   - {12 x q} matrix of identifiers of memory mapped data {q scans}
% H   - {q  x h} condition subpartition of the design matrix {h conditions}
% C   - {q  x c} covariate subpartition of the design matrix {c covariates}
% B   - {q  x n} block     subpartition of the design matrix {n subjects}
% G   - {q  x g} confound  subpartition of the design matrix {g covariates}
%
% CONTRAST - matrix of contrasts, one per row, with p elements.
% ORIGIN   - the voxel correpsonding to [0 0 0] in mm
% TH       - thresholds for each image defining voxels of interest
% Dnames   - string matrix of Dnames for effects in the design matrix
% Fnames   - string matrix of Filenames corresponding to observations
% SIGMA    - Gaussian parameter of K for correlated observations
% RT       - Repeat time for EPI fMRI (generally interscan interval)
% VMask	   - File where voxels within the brain are > 0;
%_______________________________________________________________________
%
% spm_spm is the heart of the SPM package and implements the general
% linear model in terms of a design matrix (composed of H C B
% and G) and the data (V).  The significance of the effects in H and
% C can be tested with the SPM{F}. Significant compounds of the estimated
% parameters are assessed with a quotient that has the t distribution
% under the null hyypothesis.  The resulting SPM{t} is transformed to
% the Unit Gaussian distribution [SPM{Z}] and characterized by further
% analysis using the theory of Gaussian Fields (see spm_projections.m
% for more details)
%
% The outputs of this routine are a series of .mat files containing
% paramter estimates, adjusted values, SPM{Z} etc that are written to
% CWD (see spm_defults.m).  IMPORTANT: Existing results are overwritten
% without prompting
%
% Voxels are retained for further analysis if the F ratio for that
% voxel is significant (p < UFp uncorrected) and all the voxels have a
% reasonably high activity [the threhsold is specified as a fraction
% (usually 0.8) of the whole brain mean].
%
%   SPMF.mat contains a 1 x N vector of F values reflecting the omnibus
% significance of effects [of interest] at each of the N 'significant'
% voxels.  'Significance' is defined by the p-value of the F threshold
% (p < UFp).
%
%   XYZ.mat contains a 3 x N matrix of the x,y and z location of the
% voxels in SPMF in mm (usually referring the the standard anatomical
% space (Talairach and Tournoux 1988)} (0,0,0) corresponds to the
% centre of the voxel specified by ORIGIN in the *.hdr of the original
% and related data.
%
%   BETA.mat contains a p x N matrix of the p parameter estimates at
% each of the N voxels.  These parameters include all effects
% specified by the design matrix.
%
%   XA.mat contains a q x N matrix of adjusted activity values having 
% removed the effects of no interest at each of the N voxels for all q
% scans.
%
%   SPMt.mat contains a c x N matrix of the c SPM{Z} defined by the c  
% contrasts supplied for all N voxels at locations XYZ.
%
%   SPM.mat contains a collection of matrices that pertain to the 
% analysis; including the partitions of the design matrix [H C B G], the
% number of voxels analyzed (S), image and voxel dimensions [V],
% smoothness estimates of the SPM{Z} [W], threshold for SPM{F}
% [UF] and the contrasts used [CONTRAST].  See below for a complete
% listing.
%
%   RES.mat contains 1 x N vector of residual sum of squares for the voxels
% in XYZ.mat
%
% Output to the results window includes maximum intensity projections
% of the SPM{F}, the design matrix and a series of pages for the SPM{Z}
% (see 'Results' in the help application).
%
% Variables saved in SPM.mat
%-----------------------------------------------------------------------
% H	-	condition partition of design matrix
% C	-	covariate partition of design matrix
% B	-	block     partition of design matrix
% G	-	confound  partition of design matrix
% S	-	Lebegue measure or volume {voxels}
% UF	-	Threshold for F ratio of variances
% V(1)	-	x image size {voxels}
% V(2)	-	y image size {voxels}
% V(3)	-	z image size {voxels}
% V(4)	-	x voxel size {mm}
% V(5)	-	y voxel size {mm}
% V(6)	-	z voxel size {mm}
% V(7)	-	z origin {voxels}
% V(8)	-       y origin {voxels}
% V(9)	-       z origin {voxels}
% W     -       Smoothness {Guassian parameter - voxels}
% df    -       degrees of freedom due to error
% Fdf   -       degrees of freedom for the F ratio [Fdf(2) = df]
% TH    -       vector of thresholds used to eliminate extracranial voxels
% SIGMA -       Gaussian parameter of K for correlated observations
% RT    -       Repeat time for EPI fMRI (generally interscan interval)
% BCOV  -       cov{BETA} = RES.BCOV = covariances of parameter estimates
% Dnames   -    Sting matrix of parameters in the design matrix
% Fnames   -    string matrix of Filenames corresponding to observations
% CONTRAST -    row vectors of contrasts
%
% Results matrices in .mat files (at voxels satisfying P{F > f} < UFp)
%-----------------------------------------------------------------------
% XA 	-	adjusted data  		{i.e. confounds removed}
% BETA 	-	parameter estimates	{pertaining to the design matrix}
% XYZ	-	location 		{x, y and z in mm [Talairach]}
% SPMF	-	SPM{F}
% SPMt	-	SPM{Z}
% RES 	-	residual SSQ
%
%_______________________________________________________________________
% %W% Jean-Baptiste Poline, Andrew Holmes, Karl Friston %E%

% ANALYSIS PROPER
%=======================================================================

%-Delete files from previous analyses, if they exist
%-----------------------------------------------------------------------
% spm_unlink XA.mat BETA.mat XYZ.mat SPMF.mat SPMt.mat RES.mat 


% temporal convolution of the design matrix - dispersion = SIGMA
% the block partition is deliberately omitted here (B is used to associate
% scans and subjects in subsequent routines)
%-----------------------------------------------------------------------
q	= size([H C B G],1);
bdim	= size([H C B G],2);
K     = spm_sptop(SIGMA,q);
if ~isempty(H),	H     = K*H;	end;
if ~isempty(C),	C     = K*C;	end;
if ~isempty(G),	G     = K*G;	end;

%-Critical value for F comparison at probability threshold (UFp)
%---------------------------------------------------------------------------
DESMTX = [H C B G];
Fdf    	= spm_AnCova([H C],[B G],SIGMA);
UFp	= 1;
% UF	= spm_invFcdf(1 - UFp,Fdf);
df     	= Fdf(2);
TH	= TH(:)';

%---------------
dx = V(1).dim(1);
dy = V(1).dim(2);
dz = V(1).dim(3);

% find a proper block size for the main loop
% minimum size = one line;
% maximum size = one plane;
% MaxMem is the maximum amount of data that will be processed at a time 
%-----------------------------------------------------------------------
MaxMem	= 2^20;	% In bytes
blksz	= MaxMem/8/q;
if dy < 2, error('dy < 2'); end;		% At least 2 lines
nl 	= max(min(round(blksz/dx), dy), 1); 	% nl : max # of lines
clear 	blksz;
clines	= 1:nl:dy; blines = diff([clines dy+1]);
bchsz	= length(clines);			% bchsz : number of blocks


%-- Intialise the name of the new mask : current mask & conditions on voxels
%---------------------------------------------------------------------------
VNewMsk = V(1); i=max(find(VMask.fname=='/')); if isempty(i) i=0; end;
VNewMsk.fname = [CWD '/' 'spm_msk_' VMask.fname(i+1:length(VMask.fname))];
eval(['spm_unlink ' VNewMsk.fname]);
VNewMsk.pinfo = [1 0 0]';
VNewMsk.dim(4) = spm_type('uint8');
VNewMsk.descrip = 'some info on this analysis ?';
VNewMsk = spm_create_image(VNewMsk);


%-- Intialise names of the betas 
%---------------------------------------------------------------------------
Vbeta(1:bdim) = deal(struct(...
		'fname',[],...
		'dim',[V(1).dim(1:3) spm_type('float')],...
		'mat',V(1).mat,'pinfo',[1 0 0],...
		'descrip','spm betas'));
for i=1:bdim
   Vbeta(i).fname = [CWD '/' strrep(sprintf('beta%3.0d', i),' ','0') '.img']; 
   eval(['spm_unlink ' Vbeta(i).fname]);
   Vbeta(i) = spm_create_image(Vbeta(i));
end;


%-- Intialise names of the residual sum of squares 
%---------------------------------------------------------------------------
VResSS	= struct(	'fname','ResSS','dim',[V(1).dim(1:3) 64],...
			'mat',V(1).mat,'pinfo',[1 0 0],...
			'descrip','spm Residuals Sum of Squares');
VResSS	= spm_create_image(VResSS);


%-- Intialise other variables used in the loop 
%---------------------------------------------------------------------------
AbPm	= zeros(1,dx*dy);		% above plane (mask)
vox_ind = [];				% voxels indices in plane

tmp_msk	= [];				% temporary mask 
VMaskFl	= isempty(VMask);

AbVox	= struct('res',[],'msz',[],'ofs',[],'ind',[]);   % msz : mask size
CrVox	= struct('res',[],'msz',[],'ofs',[],'ind',[]);   % ofs : mask offset
							 % ind : indices

z0	= zeros(dx*dy,1);
xA	= [1:dx]'*ones(1,V(1).dim(2));
yA	= ones(dx,1)*[1:V(1).dim(2)];
yA	= yA(:); xA	= xA(:);

S      = 0;                                     	% Volume analyzed
sx_res = 0;              
sy_res = 0;             
sz_res = 0; 
nx     = 0;
ny     = 0;
nz     = 0;
i_res  = round(linspace(1,q,min([q 64])))';		% RSSQ for smoothness

%-Cycle over groups of lines to avoid memory problems
%-----------------------------------------------------------------------
spm_progress_bar('Init',100,'AnCova',' ');

for pl_i = 1:dz	% loop over planes (2 or 3 dimensional data)

  	zA 	= pl_i + z0;
	CrBl	= [];	% current betas in plane
	CrResSS	= [];	% current betas in plane

	for bch = 1:bchsz	% loop over groups of lines 
				% bch : index of lines;

	   cl	= clines(bch); 	 	% cl = indice of line in plane
	   bl	= blines(bch);  	% bl = number of lines
	   CrLm	= ones(1,dx*bl);	% current lines (mask)

	   %----------  construct lists of voxels -----------------------

	   I = ((cl-1)*dx+1):((cl+bl-1)*dx);	% lines cl:cl+bl-1 
	   x = xA(I); y = yA(I); z = zA(I);	% x y z of blk of lines

	   %---------- get mask from disk for these lines ---------------
	   if ~isempty(VMask)
	   	CrLm		= spm_sample_vol(VMask, x,y,z,0);	
	   	CrLm		= (CrLm > 0)';
	   end
	   u = find(CrLm); mu = ones(size(u)); clear X;

	   %---------- get the data in mask -----------------------------
	   if ~isempty(u)
	      for j = 1:q
		d      	= spm_sample_vol(V(j),x(u),y(u),z(u),0)';
		mu      = mu & finite(d) & d>TH(j);	% condition FTH
		X(j,:) 	= d;
	      end
	   end
	   %---------- construct the current mask plane  ----------------
	   %-- removes voxels that didn't fulfill condition FTH

	   CrLm(u(find(~mu))) = zeros(1,sum(~mu));
	   CrVox(1).ofs = I(1) - 1;	
	   CrVox(1).ind = u(find(mu));
	   %disp(['sum(mu) ' num2str(sum(mu)) 'length(u)' num2str(length(u))]);

	   %---------- proceed with General Linear Model ----------------
	   if sum(mu) 
	   
		%-- get voxels above threshold
		X = X(:,find(mu));

		%-- volume analysed S
		S     = S + sum(mu);

		%-- Convolve over scans
		X     = K*X;

		%-AnCova; employing pinv to allow for non-unique designs	
		%---------------------------------------------------------------
		[Fdf F BETA T ResSS] = spm_AnCova([H C],[B G],SIGMA,X,CONTRAST);

		%-- save current betas in mem as we go along
		%-- if there is not enought memory, could be saved as .mad and 
		%-- get back at the end of the plane
		CrBl 	= [CrBl BETA];
		CrResSS = [CrResSS ResSS];

		% Smoothness estimation - Normalize residuals for i_res
		%---------------------------------------------------------------
		CrVox(1).res 	= X(i_res,:) - DESMTX(i_res,:)*BETA;
		CrVox(1).res	= ...
			CrVox(1).res ./ (ones(size(i_res))*sqrt(ResSS/df));

		clear X;		
		
	   	%---------- Compute spatial derivatives ...
		%-----------------------------------------------------------
		% The code avoids looping over the voxels (uses arrays)
		% and optimize memory by working on the masks.
		% It works on all the voxels contained in the mask and i_res.
		% It constructs a list of indices that correspond
		% to the voxels in the original mask AND the 
		% displaced mask. It then finds the location of these
		% in the original list of indices. If necessary, it 
		% constructs the voxel list corresponding to the displaced 
		% mask (eg in y-dim when not at the begining of the plane).
		%-----------------------------------------------------------
		%-- z dim
		%

		j_jk 	= cumsum(CrLm);			% valid for x y z
		
		if  pl_i ~= 1 %-- A plane above

		   jk	= find(AbPm(I) & CrLm);
		   if jk 	% there are some voxels above and below....
			
			k_jk	= cumsum(AbPm(I)); 
			sz_res 	= sz_res + sum(sum( ...
				  (CrVox(1).res(:,j_jk(jk)) ...
				  - AbVox(bch).res(:,k_jk(jk))).^2 ));
			nz 	= nz + length(jk);

		   end % if jk 	% there are some voxels above and below

		end % if  pl_i ~= 1


		%--------------------------------------------------------- 
		%-- y dim
		%
		if bch == 1		% first bunch : no previous line

		   if bl > 1
		   	jk = find( [CrLm(dx+1:dx*bl) zeros(1,dx)] & CrLm );
		   	if jk % some voxels...

			   k_jk = cumsum( CrLm(dx+1:dx*bl) );
			   sy_res= sy_res + sum(sum( ...
				(CrVox(1).res(:,k_jk(jk) + j_jk(dx)) - ...
			 	 CrVox(1).res(:,j_jk(jk)) ).^2  ));
				 % sum(CrLm(1:dx)) = j_jk(dx) 
		      	   ny 	= ny + length(jk);
		   	end % if jk % some voxels...
		  
		   end % bl > 1

		else % if bch == 1 -> a previous line exists

		   %-- Append previous line and CrLm 
		   %-- minus its last line in tmp_msk

		   tmp_msk = [AbPm( ((cl-2)*dx+1):(cl-1)*dx ) ...
				CrLm(1:dx*(bl-1)) ];

		   %-- get tmp_msk corresponding voxels tmp_vox
		   tmp_vox = ...
			[AbVox(bch-1).res(:,find(AbVox(bch-1).ind>dx*(nl-1))) ...
			 CrVox(1).res(:, find(CrVox(1).ind <= dx*(bl-1) )) ];
	
		   jk	= find( CrLm & tmp_msk);
		   if jk % some voxels...
		      k_jk	= cumsum(tmp_msk);
		      sy_res 	= sy_res + sum(sum( ...
		   		  (CrVox(1).res(:,j_jk(jk)) - ...
				  tmp_vox(:,k_jk(jk))).^2  ));
		      ny 	= ny + length(jk);
		   end % if jk % some voxels...

		end % if bch == 1

		%--------------------------------------------------------- 
		%-- x dim
		%
		%-- 1. shift the mask to the left (arbitrary)
		%-- 2. add 0 at the end of the lines

		tmp_msk = [CrLm(2:length(CrLm)) 0];
		tmp_msk(dx:dx:bl*dx) = zeros(1,bl);
		
		jk	= find(tmp_msk & CrLm);
		if jk 	% there are some voxels on the left ....

		   sx_res = sx_res + sum(sum( ...
		   	(CrVox(1).res(:,j_jk(jk)+1) - ...
			 CrVox(1).res(:,j_jk(jk))).^2  ));
		   	% j_jk(jk)+1 = position of voxels on the right
			% of tmp_msk & CrLm
		   nx 	 = nx + length(jk);

		end % if jk 	% there are some voxels on the left ....|

	   %--------------------------------------------------------- 
	   %---------- compute partial derivatives ...              |
	   %--------------------------------------------------------- 

	   end % if sum(mu) 

	   %---------- roll ...
	   %-- AbVox(bch) is overwritten 

	   AbVox(bch).ind = CrVox(1).ind;
	   AbVox(bch).ofs = CrVox(1).ofs;
	   AbVox(bch).res = CrVox(1).res;

	   AbPm(I) = CrLm;
	   %---------- roll ...|

	end %  bch = 1:bchsz

	%------------------------ last bunch of lines =>  all plane in AbVox
	%-- write  mask betas ResSS in file

	spm_write_plane(VNewMsk, reshape(AbPm,dx,dy), pl_i);

	%-- first find the indices in AbVox
	vox_ind = [];
	for i_tmp = 1:bchsz
	   vox_ind = [vox_ind AbVox(i_tmp).ind+AbVox(i_tmp).ofs];
	end
	
	for i_beta =1:bdim
	   bet_tmp	= zeros(1,dx*dy);	
	   if ~isempty(vox_ind)
	  	bet_tmp(vox_ind) = CrBl(i_beta,:); 
	   end;
	   spm_write_plane(Vbeta(i_beta), reshape(bet_tmp,dx,dy),pl_i);
	end % for i_beta =1:bdim
		
	%-- idem for ResSS
	bet_tmp	= zeros(1,dx*dy);	
	bet_tmp(vox_ind) = CrResSS;
	spm_write_plane(VResSS,reshape(bet_tmp,dx,dy),pl_i);		   

	   
	spm_progress_bar('Set',100*(bch + bchsz*(pl_i-1))/(bchsz*dz));

disp(['plane --------------------------------- ' num2str(pl_i)]);

end % for pl_i = 1:dz : loop over planes
%-----------------------------------------------------------------------

spm_progress_bar('Clear');


%-Smoothness estimates %-----------------------------------------------------------------------
Lc2z   = spm_lambda(df);
if dz == 1,
	if any(~[nx ny], error(['W: nx ny ' num2str([nx ny])]); end;
	L_res  = [sx_res/nx sy_res/ny]*(df - 2)/(df - 1)/df;
else
	if any(~[nx ny nz], error(['W: nx ny nz' num2str([nx ny nz])]); end;
	L_res  = [sx_res/nx sy_res/ny sz_res/nz]*(df - 2)/(df - 1)/df;
end
W      = (2*Lc2z*L_res).^(-1/2)
W = W(:,[1:2]); end; %--- 2 dimensional data
FWHM   = sqrt(8*log(2))*W.*sqrt(sum(V(1).mat(1:3,1:3).^2));


%-Save design matrix, and other key variables; S UF CONTRAST W V and df
%-----------------------------------------------------------------------
% V      = [V(1:6,1); ORIGIN(:)];
% [Fdf,F,BETA,T,RES,BCOV] = spm_AnCova([H C],[B G],SIGMA);
%save SPM H C B G S UF V W CONTRAST df Fdf TH Dnames Fnames SIGMA RT BCOV
%save SPM H C B G S V W CONTRAST df Fdf TH Dnames Fnames SIGMA RT BCOV


%-Display and print SPM{F}, Design matrix and textual information
%=======================================================================
FWHM   = sqrt(8*log(2))*W.*sqrt(sum(V(1).mat(1:3,1:3).^2));
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph); spm_clf(Fgraph)
if exist('SPMF.mat')
	load XYZ
	load SPMF
	axes('Position',[-0.05 0.5 0.8 0.4]);
	spm_mip(sqrt(SPMF),XYZ,V)
	str = sprintf('SPM{F} p < %0.2f, df: %0.1f,%0.1f',UFp,Fdf);
	title(str,'FontSize',16)
end
text(240,220,sprintf('Search volume: %d voxels',S))
text(240,240,sprintf('Image size: %d %d %d voxels',V(1:3)))
text(240,260,sprintf('Voxel size  %0.1f %0.1f %0.1f mm',V(4:6)))
text(240,280,sprintf('Resolution {FWHM} %0.1f %0.1f %0.1f mm',FWHM))



%-Print out contrasts
%-----------------------------------------------------------------------
axes('Position',[0.1 0.1 0.8 0.4],'XLim',[0,1],'YLim',[0,1]); axis off
line([0 1],[1 1],'LineWidth',3);
line([0 1],[0.92,0.92],'LineWidth',3);
line([0 1],[0.85 0.85],'LineWidth',1);

text(0,0.96,'Results directory:')
text(0.23,0.96,pwd,'FontSize',12,'Fontweight','Bold')
text(0,0.88,'Contrasts','FontSize',12,'Fontweight','Bold')

x0 = 0.25; y0 = 0.82; dx = 0.8/size(CONTRAST,1); dy = 0.04;
for j = 1:size([H C],2)
	text(0,y0 - (j - 1)*dy,Dnames(j,:),'FontSize',10)
end

line([0 1],[1,1]*(y0 - size([H C],2)*dy),'LineWidth',1);
for i = 1:size(CONTRAST,1)
	text(x0 + dx*(i - 1),0.88,int2str(i),'FontSize',10)
	for j = 1:size([H C],2)
		str = sprintf('%-6.3g',CONTRAST(i,j));
		text(x0 + dx*(i - 1),y0 - (j - 1)*dy,str,'FontSize',10)
	end
end

spm_print

%-Display, characterize and print SPM{Z}
%-----------------------------------------------------------------------
if exist('SPMt.mat')
	load SPMt
	U     = spm_invNcdf(1 - 0.01);
	% [P,EN,Em,En] = spm_P(1,W,U,0,S);
	% K     = round(En);
	K = 0;
	for i = 1:size(CONTRAST,1)
	    spm_projections(SPMt(i,:),XYZ,U,K,V,W,S,DESMTX,CONTRAST(i,:),df);
	    spm_print
	end
end




