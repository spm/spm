function spm_spm(VY,X,Xnames,xM,K,c,UFp,varargin)
%
% FORMAT spm_spm(VY,X,Xnames,xM,K,c,UFp,<extra parameters for SPM.mat>)
% FORMAT spm_spm(VY,X,Xnames,TH,K,c,UFp,<extra parameters for SPM.mat>)
%
% First contrast in xCon is used for F-thresholding
%
% K should be a sparse smoothing matrix, which should *not* smooth over
% block boundaries (creates havoc with contrast validity!)
%
% Raw (unsmoothed) data saved in squashed *.mad format (which takes a
% while, seeing as there's no F-thresholding any more, and that the
% *.mad file format involves computation of scale and offsets for each
% voxel!)
%
% Y.mad file saved only for plotting
%
% Christensen reference for F-contrasts
%
%-** Need to make sure results section can obtain relevant info from SPMcfg file
%_______________________________________________________________________
% %W% Jean-Baptiste Poline, Andrew Holmes, Karl Friston %E%
SCCSid   = '%I%';

%-Parameters
%-----------------------------------------------------------------------
Def_UFp  = 0.0001;	%-Default F-threshold for Y.mad pointlist filtering
maxMem   = 2^20;	%-Max data chunking size, in bytes
maxRes4S = 64;		%-Maximum #res images for smoothness estimation

%-Condition arguments
%-----------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end
if nargin<4, xM  = zeros(size(X,1),1); end
if nargin<5, K   = speye(size(X,1)); end
if nargin<6, c   = []; end
if nargin<7, UFp = []; end
if isempty(UFp), UFp = Def_UFp; end

%-If xM is not a structure then assumme it's a vector of thresholds
if ~isstruct(xM), xM = struct(	'T',	[],...
				'TH',	xM,...
				'I',	0,...
				'VM',	{[]},...
				'xs',	struct('Masking','analysis threshold'));
end


%-Say hello
%-----------------------------------------------------------------------
SPMid  = spm('FnBanner',mfilename,SCCSid);
Finter = spm('FigName','Stats: estimation'); spm('Pointer','Watch')


%-Delete files from previous analyses...
%-----------------------------------------------------------------------
if exist('./SPM.mat','file')==2
	warning(sprintf(['Existing results in this directory will be ',...
		'overwritten\n\t (pwd = %s) '],pwd))
end
spm_unlink SPM.mat XYZ.mat Y.mad
spm_unlink mask.img mask.hdr ResSS.img ResSS.hdr
[files,null] = spm_list_files(pwd,'beta_????.*');
for i=1:size(files,1), spm_unlink(deblank(files(i,:))), end


%=======================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%=======================================================================

%-Check Y images have same dimensions, orientation & voxel size
%-----------------------------------------------------------------------
if any(any(diff(cat(1,VY.dim),1,1),1)&[1,1,1,0]) %NB: Bombs for single image
	error('images do not all have the same dimensions'), end
if any(any(any(diff(cat(3,VY.mat),1,3),3)))
	error('images do not all have same orientation & voxel size'), end


%-Condition temporal convolution matrix
%-----------------------------------------------------------------------
if any(size(K)~=size(X,1)), error('K not a temporal smoothing matrix'), end
%-** Check K doesn't smooth the block (B) partition?


%-Design parmameters
%-----------------------------------------------------------------------
[nScan,nbeta] = size(X);		%-#scans & #parameters
xXd           = spm_glm('KCC',X,K);	%-Design matrix psudoinverses & traces


%-F-contrast for Y.mad pointlist filtering (only)
%-----------------------------------------------------------------------
%-**if no F-contrast specified, do raw F for all effects (i.e. vs. raw data SS)
% if ~UFp (UFp==0) => don't write Y.mad
%     UFp == 1     => no F-filtering
% o/w              => F-filtering of data at upper tail probability p = UFp

if isempty(c) & UFp>0 & UFp<1		%-Want to F-threshold but no contrast!
	UFp=1;
	warning('no F-con specified: no F-filtering')
end

if UFp>0 & UFp<1			%-We're going to F-filter for Y.mad file
	if ~spm_DesUtil('allCon',X,c), error('Invalid F-contrast'), end
	R     = c'*pinv(c*pinv(X'*X)*c')*c;	%-See Christensen
	Fdf   = pi;				%-****
	trROV = Fdf;				%-****
	%-Compute F-con mtx, df, trROV

	%-Initialise XYZ location vector
	%---------------------------------------------------------------
	XYZ = [];
elseif UFp == 1
	UF = 0;
elseif UFp == 0
	
end



%-Image dimensions
%-----------------------------------------------------------------------
xdim = VY(1).dim(1); ydim = VY(1).dim(2); zdim = VY(1).dim(3);
YNaNrep = spm_type(VY(1).dim(4),'nanrep');


%-Intialise the name of the new mask : current mask & conditions on voxels
%-----------------------------------------------------------------------
VM = struct(		'fname',	'mask',...
			'dim',		[VY(1).dim(1:3),spm_type('uint8')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:resultant analysis mask');
VM = spm_create_image(VM);


%-Intialise beta image files
%-----------------------------------------------------------------------
Vbeta(1:nbeta) = deal(struct(...
			'fname',	[],...
			'dim',		[VY(1).dim(1:3) spm_type('float')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	''));
for i=1:nbeta
	Vbeta(i).fname = sprintf('beta_%04d.img',i);
	Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,Xnames{i});
	spm_unlink(Vbeta(i).fname)
	Vbeta(i) = spm_create_image(Vbeta(i));
end


%-Intialise residual sum of squares image file
%-----------------------------------------------------------------------
VResSS = struct(	'fname',	'ResSS',...
			'dim',		[VY(1).dim(1:3) spm_type('double')],...
			'mat',		VY(1).mat,...
			'pinfo',	[1 0 0]',...
			'descrip',	'spm_spm:Residual Sum of Squares');
VResSS	= spm_create_image(VResSS);




%=======================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%=======================================================================

%-Find a suitable block size for the main loop, which proceeds a bunch
% of lines at a time (minimum = one line; maximum = one plane)
% (maxMem is the maximum amount of data that will be processed at a time)
%-----------------------------------------------------------------------
blksz	= maxMem/8/nScan;			%-block size (in bytes)
if ydim<2, error('ydim < 2'), end		%-need at least 2 lines
nl 	= max(min(round(blksz/xdim),ydim),1); 	%-max # lines / block
clines	= 1:nl:ydim;				%-bunch start line #'s
blines  = diff([clines ydim+1]);		%-#lines per bunch
nbch    = length(clines);			%-#bunches


%-Intialise other variables used in the loop 
%=======================================================================
AbPm	= zeros(1,xdim*ydim);				%-above plane (mask)
vox_ind = [];						%-voxels indices in plane

tmp_msk	= [];						%-temporary mask 

AbVox	= struct('res',[],'ofs',[],'ind',[]);		%-voxels above
CrVox	= struct('res',[],'ofs',[],'ind',[]);		%-current voxlels
							% res : residuals
   							% ofs : mask offset
							% ind : indices

xords	= [1:xdim]'*ones(1,ydim); xords = xords(:)';	%-plane X coordinates
yords	= ones(xdim,1)*[1:ydim];  yords = yords(:)';	%-plane Y coordinates


%-Smoothness estimation variables
%-----------------------------------------------------------------------
S      = 0;                                     %-Volume analyzed (in voxels)
sx_res = 0; sy_res = 0; sz_res = 0;		%-sum((dr./d{x,y,z}).^2)
nx     = 0; ny     = 0; nz     = 0;		%-# {x,y,z} partial derivs
%-Indexes of residual images to sample for smoothness estimation
i_res  = round(linspace(1,nScan,min(nScan,maxRes4S)))';


%-Cycle over groups of lines within planes to avoid memory problems
%=======================================================================
spm_progress_bar('Init',100,'model estimation','');

for pl_i = 1:zdim			%-loop over planes (2D or 3D data)

    zords   = pl_i*ones(xdim*ydim,1)';	%-plane Z coordinates
    CrBl    = [];			%-current plane betas
    CrResSS = [];			%-current plane ResSS
    
    for bch = 1:nbch			%-loop over bunches of lines

	fprintf('\r\tPlane %3d/%-3d : plank %3d/%-3d%20s',pl_i,zdim,bch,nbch,' ')

	cl   = clines(bch); 	 	%-line index of first line of bunch
	bl   = blines(bch);  		%-number of lines for this bunch

	%-construct list of voxels in this bunch of lines
	I    = ((cl-1)*xdim+1):((cl+bl-1)*xdim);	%-lines cl:cl+bl-1
	xyz  = [xords(I); yords(I); zords(I)];		%-coords of voxels in bch


	%-Get data & construct analysis mask for this bunch of lines
	%===============================================================
	fprintf('%s%20s',sprintf('\b')*ones(1,20),'...read & mask data')%-#
	CrLm = logical(ones(1,xdim*bl));		%-current lines mask

	%-Compute explicit mask
	% (note that these may not have same orientations as Y images)
	%---------------------------------------------------------------
	for i = 1:length(xM.VM)
		tM   = inv(xM.VM(i).mat)*VY(1).mat;	%-Reorientation matrix
		tmp  = tM * [xyz;ones(1,size(xyz,2))];	%-Coords in mask image
		%-Load mask image within current mask & update mask
		CrLm(CrLm) = spm_sample_vol(xM.VM(i),...
				tmp(1,CrLm),tmp(1,CrLm),tmp(3,CrLm),0) > 0;
	end
	
	%-Get the data in mask, compute threshold & implicit masks
	%---------------------------------------------------------------
	Y = zeros(nScan,xdim*bl);
	for j = 1:nScan
		if ~any(CrLm), break, end		%-Break if empty mask
		Y(j,CrLm) = spm_sample_vol(VY(j),...	%-Load data in mask
				xyz(1,CrLm),xyz(2,CrLm),xyz(3,CrLm),0);
		CrLm(CrLm) = Y(j,CrLm)>xM.TH(j);	%-Threshold (& NaN) mask
		if xM.I & ~YNaNrep & xM.TH(j)<0		%-Use implicit 0 mask
			CrLm(CrLm) = abs(Y(j,CrLm))>eps;
		end
	end

	%-Apply mask
	%---------------------------------------------------------------
	Y = Y(:,CrLm);				%-Data matrix within mask
	S = S + sum(CrLm);			%-Volume analysed
	CrVox.ofs = I(1) - 1;			%-Lines offset
	CrVox.ind = find(CrLm);			%-Voxel indexes (within bunch)



	%-Proceed with General Linear Model & smoothness estimation
	%===============================================================
	if any(CrLm)
	
		%-Save raw data in squashed *.mad file format (see spm_append)
		% (This is done so that the data is available for plotting later)
		% (Idea is that results are self contained, but this will be a  )
		% (big file, so might consider using the raw data images, and   )
		% (assumme/hope they don't move or get changed!                 )
		fprintf('%s%20s',sprintf('\b')*ones(1,20),'...saving data')%-#
		%-**spm_append(Y,'Y.mad',2)

		%-save xyz co-ordinates in matrix for saving later
		XYZ = [XYZ, xyz];

		fprintf('%s%20s',sprintf('\b')*ones(1,20),'...par est.')%-#

		%-Temporal smoothing (done in place to save memory)
		%-------------------------------------------------------
		Y = K*Y;			%-(so now Y should read KY :-)

		%-General linear model: least squares estimation
		% (Using pinv to allow for non-unique designs             )
		% (Including temporal convolution of design matrix with K )
		% (Note that Y has already been temporally smoothed       )
		%-------------------------------------------------------
		beta  = xXd.PKX * Y;			%-Parameter estimates
		res   = Y - xXd.KX*beta;		%-Residuals
		ResSS = sum(res.^2);			%-Residual sum of squares
		%-**[Fdf F betA T ressS] = spm_AnCova(K*X,[],K,Y);
		clear Y

		%-Save betas for current plane in memory as we go along
		% (if there is not enough memory, could save directly as float )
		% (analyze images, or *.mad and retreived when plane complete  )
		CrBl 	= [CrBl,beta];
		CrResSS = [CrResSS,ResSS];

		%-Subsample (i_res) normalised residuals for moothness estimation
		%-------------------------------------------------------
		CrVox.res=res(i_res,:)./(ones(size(i_res))*sqrt(ResSS/xXd.trRV));
		clear res				%-clear res


	   	%-Smoothness estimation: compute spatial derivatives...
		%=======================================================
		%-The code avoids looping over the voxels (using arrays)
		% and optimizes memory usage by working on the masks.
		%-It works on all the voxels contained in the mask and i_res.
		%-It constructs a list of indices that correspond
		% to the voxels in the original mask AND the 
		% displaced mask. It then finds the location of these
		% in the original list of indices. If necessary, it 
		% constructs the voxel list corresponding to the displaced 
		% mask (eg in y-dim when not at the begining of the plane).
		fprintf('%s%20s',sprintf('\b')*ones(1,20),'...smoothness est')%-#


		%-Construct utility vector to map from mask (CrLm) to voxel
		% index within the mask. (For voxels with CrLm true, j_jk
		% is the index of the voxel within the mask)
		j_jk = cumsum(CrLm);			% (valid for x y z)


		%- z dim
		%-------------------------------------------------------
		if pl_i>1				%-There's a plane below
		    %-Indices of inmask voxels with inmask z-1 neighbours
		    jk = find(CrLm & AbPm(I));
		    if ~isempty(jk)
			%-Compute indices of inmask z-1 voxels
			k_jk   = cumsum(AbPm(I)); 
			sz_res = sz_res + sum(sum( ...
					(CrVox.res(:,j_jk(jk)) ...
					- AbVox(bch).res(:,k_jk(jk))).^2 ));
			nz     = nz + length(jk);
		   end % (if ~isempty(jk))
		end % (if pl_i~=1)


		%- y dim
		%-------------------------------------------------------
		if bch==1		%-first bunch : no previous line
		    if bl > 1
			%-Indices of inmask voxels with inmask y+1 neighbours
		   	jk = find(CrLm & [CrLm(xdim+1:xdim*bl),zeros(1,xdim)]);
		   	if ~isempty(jk)
			    %-Compute indices of inmask y+1 voxels
			    k_jk   = cumsum(CrLm(xdim+1:xdim*bl))+j_jk(xdim);
			    % NB: sum(CrLm(1:xdim))==j_jk(xdim)
			    %-Compute partial derivs, sum ^2's over resids & vox
			    sy_res = sy_res + sum(sum( ...
				(CrVox.res(:,k_jk(jk)) - ...
			 	 CrVox.res(:,j_jk(jk)) ).^2  ));
		      	   ny      = ny + length(jk);
		   	end % (if ~isempty(jk))
		    end % (if bl > 1)


		else % (if bch==1)	%-a previous line exists
		    %-Make tmp_mask as mask shifted y-1, by prepending previous
		    % line to CrLm minus its last line
		    tmp_msk = [AbPm( ((cl-2)*xdim+1):(cl-1)*xdim ),...
				CrLm(1:xdim*(bl-1)) ];
		    %-get residuals for y-1 shifted space (of tmp_vox)
		    tmp_vox = ...
		      [AbVox(bch-1).res(:,find(AbVox(bch-1).ind>xdim*(nl-1))),...
		       CrVox.res(:,find(CrVox.ind <= xdim*(bl-1))) ];
	
		    %-Inmask voxels with inmask (incl. prev line) y-1 neighbours
		    jk	= find(CrLm & tmp_msk);
		    if ~isempty(jk)
			%-Compute indices of inmask y-1 voxels
			k_jk	= cumsum(tmp_msk);
			sy_res 	= sy_res + sum(sum( ...
					(CrVox.res(:,j_jk(jk)) - ...
					tmp_vox(:,k_jk(jk))).^2  ));
			ny 	= ny + length(jk);
		    end % (if ~isempty(jk))

		end % (if bch==1)


		%- x dim
		%-------------------------------------------------------
		%-Shift the mask to the left, add 0 at line ends
		tmp_msk = [CrLm(2:length(CrLm)),0];
		tmp_msk(xdim:xdim:bl*xdim) = 0;
		
		%-Indices inmask of voxels with inmask x+1 neighbours
		jk = find(CrLm & tmp_msk);
		if ~isempty(jk)
		   sx_res = sx_res + sum(sum( ...
		   	(CrVox.res(:,j_jk(jk)+1) - ...
			 CrVox.res(:,j_jk(jk))).^2  ));
		   % NB: j_jk(jk)+1 = position of voxels on right of j_jk(jk)
		   nx     = nx + length(jk);
		end % (if ~isempty(jk))

	    end % (if any(CrLm))

	    %-Roll... (AbVox(bch) is overwritten)
	    %-----------------------------------------------------------
	    AbVox(bch).ind = CrVox.ind;		%-Voxel indexes (within bunch)
	    AbVox(bch).ofs = CrVox.ofs;		%-Bunch voxel offset (in plane)
	    AbVox(bch).res = CrVox.res;		%-Sample of residuals
	    AbPm(I)        = CrLm;		%-"above plane" mask

    end % (for bch = 1:nbch)
    

    %-Plane complete, write out plane data to image files
    %===================================================================
    fprintf('%s%20s',sprintf('\b')*ones(1,20),'...saving plane')%-#

    %-Mask image
    %-------------------------------------------------------------------
    %-AbPm (& AbVox) now contains a complete plane mask
    spm_write_plane(VM, reshape(AbPm,xdim,ydim), pl_i);

    %-Construct voxel indices for AbPm
    Q   = find(AbPm);
    tmp = zeros(xdim,ydim);


    %-Write beta images
    %-------------------------------------------------------------------
    for i = 1:nbeta
	tmp(Q) = CrBl(i,:);
	spm_write_plane(Vbeta(i),tmp,pl_i);
    end % (for i = 1:nbeta)
	    
    %-Write ResSS image
    %-------------------------------------------------------------------
    tmp(Q) = CrResSS;
    spm_write_plane(VResSS,tmp,pl_i);		   

	   
    %-Report progress
    fprintf('%s%20s',sprintf('\b')*ones(1,20),' ')%-#
    spm_progress_bar('Set',100*(bch + nbch*(pl_i-1))/(nbch*zdim));

end % (for pl_i = 1:zdim)
sprintf('\n')


%=======================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%=======================================================================

%-Smoothness estimates (For Gaussianised SPM{t} : SPM{G-t})
%-----------------------------------------------------------------------
df = xXd.df;
Lc2z   = spm_lambda(df);
if zdim == 1,
	if any(~[nx,ny]), error(sprintf('W: nx=%d, ny=%d',nx,ny)), end;
	Lambda = diag([sx_res/nx, sy_res/ny, Inf]*(df-2)/(df-1)/xXd.trRV);
else
	if any(~[nx,ny,nz]),error(sprintf('W: nx=%d, ny=%d, nz=%d',nx,ny,nz)),end
	Lambda = diag([sx_res/nx, sy_res/ny, sz_res/nz]*(df-2)/(df-1)/xXd.trRV);
end
W      = (2*Lc2z*Lambda).^(-1/2);
%-**FWHM   = sqrt(8*log(2))*W.*sqrt(sum(VY(1).mat(1:3,1:3).^2));


%-Save remaining results files and analysis parameters
%=======================================================================

%-Save coordinates of within mask voxels (for Y.mad pointlist use)
%-----------------------------------------------------------------------
save XYZ XYZ

%-Save analysis parameters in SPM.mat file
%-----------------------------------------------------------------------
SPMvars = {	'SPMid','VY','X','Xnames','xM','K',...
		'xXd','S','Lambda','W'};
if nargin>5
	%-Extra arguments were passed for saving in the SPM.mat file
	for i=7:nargin
		SPMvars = [SPMvars,{inputname(i)}];
		eval([inputname(i),' = varargin{i-6};'])
	end
end
save('SPM',SPMvars{:})

%-End: Cleanup GUI
%=======================================================================
spm_progress_bar('Clear')
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('\n\n')




keyboard
return





%=======================================================================
% - O L D   D I S P L A Y   S T U F F
%=======================================================================


%-Display and print SPM{F}, Design matrix and textual information
%=======================================================================
FWHM   = sqrt(8*log(2))*W.*sqrt(sum(VY(1).mat(1:3,1:3).^2));
Fgraph = spm_figure('FindWin','Graphics');
figure(Fgraph); spm_clf(Fgraph)
if exist('SPMF.mat')
	load XYZ
	load SPMF
	axes('Position',[-0.05 0.5 0.8 0.4]);
	spm_mip(sqrt(SPMF),XYZ,VY)
	str = sprintf('SPM{F} p < %0.2f, df: %0.1f,%0.1f',UFp,Fdf);
	title(str,'FontSize',16)
end
text(240,220,sprintf('Search volume: %d voxels',S))
text(240,240,sprintf('Image size: %d %d %d voxels',VY(1:3)))
text(240,260,sprintf('Voxel size  %0.1f %0.1f %0.1f mm',VY(4:6)))
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

x0 = 0.25; y0 = 0.82; xdim = 0.8/size(CONTRAST,1); ydim = 0.04;
for j = 1:size([H C],2)
	text(0,y0 - (j - 1)*ydim,Dnames(j,:),'FontSize',10)
end

line([0 1],[1,1]*(y0 - size([H C],2)*ydim),'LineWidth',1);
for i = 1:size(CONTRAST,1)
	text(x0 + xdim*(i - 1),0.88,int2str(i),'FontSize',10)
	for j = 1:size([H C],2)
		str = sprintf('%-6.3g',CONTRAST(i,j));
		text(x0 + xdim*(i - 1),y0 - (j - 1)*ydim,str,'FontSize',10)
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
	    spm_projections(SPMt(i,:),XYZ,U,K,VY,W,S,X,CONTRAST(i,:),df);
	    spm_print
	end
end




