function spm_segment(PF,PG,opts)
% Segment an MR image into Gray, White, CSF & other.
%
% --- The Prompts Explained ---
%
% 'select MRI(s) for subject '
% If more than one volume is specified (eg T1 & T2), then they must be
% in register (same position, size, voxel dims etc..).
%
% 'select Template(s) '
% If the images have been spatially normalised, then there is no need to
% select any images. Otherwise, select one or more template images which
% will be used for affine normalisation of the images. Note that the
% affine transform is only determined from the first image specified
% for segmentation. 
%
%_______________________________________________________________________
%
%                      The algorithm is three step:
%
% 1) Determine the affine transform which best matches the image with a
%    template image. If the name of more than one image is passed, then
%    the first image is used in this step. This step is not performed if
%    no template images are specified.
%    Note: the templates and probability images are Right-is-Right
%    orientation. The affine normalisation assumes that the input images
%    are the same - although the algorithm works well enough with input
%    images in the Left-is-Right orientation.
%
% 2) Perform Cluster Analysis with a modified Mixture Model and a-priori
%    information about the likelihoods of each voxel being one of a
%    number of different tissue types. If more than one image is passed,
%    then they they are all assumed to be in register, and the voxel
%    values are fitted to multi-normal distributions.
%
% 3) Write the segmented image. The names of these images have
%    "_seg1", "_seg2" & "_seg3" appended to the name of the
%    first image passed.
%
%_______________________________________________________________________
%
% Ref:
% Hartigan, J. A. 1975. Clustering Algorithms, pp. 113-129.
% John Wiley & Sons, Inc., New York.
%
%_______________________________________________________________________
%
% The template image, and a-priori likelihood images are modified
% versions of those kindly supplied by Alan Evans, MNI, Canada
% (ICBM, NIH P-20 project, Principal Investigator John Mazziotta).
%_______________________________________________________________________
% %W% (c) John Ashburner %E%

% Programmers notes
%
% FORMAT spm_segment(PF,PG,opts)
% PF   - name(s) of image(s) to segment (must have same dimensions).
% PG   - name(s) of template image(s) for realignment.
%      - or a 4x4 transformation matrix which maps from the image to
%        the set of templates.
% opts - options string.
%        - 't' - write images called *_seg_tmp* rather than *_seg*
%        - 'f' - fix number of voxels in each cluster

debug = 0;

if nargin<3
	% set to ' ' rather than '' to get rid of annoying warnings
	opts = ' ';
	if nargin<2
		PG = '';
	end
end

global SWD
DIR1   = [SWD '/templates/'];
DIR2   = [SWD '/apriori/'];

if (nargin==0)
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Segmentation')
	spm_help('!ContextHelp','spm_segment.m');

	n     = spm_input('number of subjects',1,'e',1);
	if (n < 1)
		spm_figure('Clear','Interactive');
		return;
	end


	for i = 1:n
		PF = spm_get(Inf,'.img',...
			['Select MRI(s) for subject ' num2str(i)]);
		eval(['PF' num2str(i) ' = PF;']);
	end

	PG = '';
	tmp = spm_input('Are they spatially normalised?', 2, 'y/n');

	if (tmp == 'n')
		% Get template
		%-----------------------------------------------------------------------
		templates = str2mat([DIR1 'T1.img'], [DIR1 'T2.img']);

		% Get modality of target
		respt = spm_input('Modality of first image?',3,'m',...
			'modality - T1 MRI|modality - T2 MRI|--other--',...
			[1 2 0],1);
		if (respt > 0)
			PG = deblank(templates(respt,:));
		else
			ok = 0;
			while (~ok)
				PG = spm_get(Inf,'.img',['Select Template(s) for normalisation'],...
					'', DIR1);
				if (size(PG,1)>0)
					dims = zeros(size(PG,1),9);
					for i=1:size(PG,1)
						[dim vox dummy dummy dummy origin dummy]...
							 = spm_hread(deblank(PG(i,:)));
						dims(i,:) = [dim vox origin];
					end
					if size(dims,1) == 1 | ~any(diff(dims))
						ok = 1;
					end
				else
					ok = 0;
				end
			end
		end
	end

	for i = 1:n
		eval(['PF = PF' num2str(i) ';']);
		set(spm_figure('FindWin','Interactive'),'Name',...
			['Segmenting ' num2str(i) '..'],'Pointer','Watch');
		drawnow;
		if (size(PF,1)~=0) spm_segment(PF,PG,opts); end
	end
	spm_figure('Clear','Interactive');
	return;
end

%_______________________________________________________________________
%_______________________________________________________________________

%- A-Priori likelihood images.
PB    = str2mat([DIR2 'gray.img'],[DIR2 'white.img']);

niter     = 48;
nc        = [1 1 3]; % Number of clusters for each probability image

% Starting parameters for affine normalisation.
% - assume input image uses neurological convention

for i=1:size(PF,1)
	M = spm_get_space(PF(i,:));
	MF(:,i) = M(:);
end
MF1 = reshape(MF(:,1),4,4);

if ~isempty(PG) & isstr(PG)
	% Affine normalisation
	MG = spm_get_space(PG(1,:));
	spm_smooth(PF(1,:),'./spm_seg_tmp.img',8);

	% perform affine normalisation at different sampling frequencies
	% with increasing numbers of parameters.
	disp('Determining Affine Mapping');
	global sptl_Ornt
	if prod(size(sptl_Ornt)) == 12
		prms = [sptl_Ornt ones(1,size(PG,1))]';
	else
		prms = [0 0 0  0 0 0  1 1 1  0 0 0 ones(1,size(PG,1))]';
	end
	spm_chi2_plot('Init','Affine Registration','Convergence');
	prms = spm_affsub3('affine3', PG, './spm_seg_tmp.img', 1, 8, prms);
	prms = spm_affsub3('affine3', PG, './spm_seg_tmp.img', 1, 6, prms);
	spm_chi2_plot('Clear');
	spm_unlink ./spm_seg_tmp.img ./spm_seg_tmp.hdr ./spm_seg_tmp.mat
	MM = spm_matrix(prms);

elseif all(size(PG) == [4 4])
	MM = PG;
else
	MM = spm_matrix([0 0 0 0 0 0 1 1 1 0 0 0]');
end


for i=1:size(PF,1)
	VF(:,i) = spm_map(PF(i,:));
end

for j=1:size(PB,1)
	VB(:,j)  = spm_map(PB(j,:));
	tmp      = spm_get_space(PB(j,:));
	Mat(:,j) = tmp(:);
end

% Cluster Analysis
%-----------------------------------------------------------------------
% Voxels to sample during the cluster analysis
sxy    = 1/4;
sz     = 6;

bb1    = [ [-78 78]' [-112 76]' [-50 85]'];

faces = MM*MF1*[
1	VF(2,1)/2	VF(3,1)/2	1
VF(1,1)	VF(2,1)/2	VF(3,1)/2	1
VF(1,1)/2	1	VF(3,1)/2	1
VF(1,1)/2	VF(2,1)	VF(3,1)/2	1
VF(1,1)/2	VF(2,1)/2	1	1
VF(1,1)/2	VF(2,1)/2	VF(3,1)	1]';

faces = faces(1:3,:)';
bb    = [max([min(faces); bb1(1,:)]); min([max(faces) ; bb1(2,:)])];
sb    = diff(bb);
d     = round(sb(1:2)*sxy);


lkp=[];
for i=1:size(nc,2)
	lkp = [lkp ones(1,nc(i))*i];
end

n  = size(lkp,2);
m  = size(PF,1);
nb = size(VB,2);

mn = zeros(m,n);	% Means
cv = zeros(m*m,n);	% (Co)variances
mg = zeros(1,n);	% Number of voxels/cluster

sumbp = zeros(1,n);

osumpr = -Inf;

spm_chi2_plot('Init','Segmenting','Log-likelihood','Iteration #');

disp('Segmenting'); 
for iter = 1:niter
	sumpr= 0;
	mom0 = zeros(1,n)+eps;
	mom1 = zeros(m,n);
	mom2 = zeros(m*m,n)+eps;

	for p = bb(1,3):sz:bb(2,3)

		clear dat bp pr

		B     = spm_matrix([-bb(1,1)*sxy -bb(1,2)*sxy -p 0 0 0 sxy sxy 1]);

		% Ignore voxels of value zero - since we don't know if it is because they
		% are truly zero - or if it is because they are zeros due to some kind of
		% image editing or from outside the original FOV
		msk = zeros(d(1)*d(2),1);

		for i=1:m
			% Sample the image data using
			% nearest neighbour to reduce smoothing effects.
			M = inv(B*MM*reshape(MF(:,i),4,4));
			dt = spm_slice_vol(VF(:,i),M, d(1:2),1);
			msk = msk | dt(:)==0;

			% Most cases will not have -ve values - but a few (e.g., PET will have)
			% - so get rid of them.
			dat(:,i) = abs(dt(:));
		end
		msk        = find(msk);

		if prod(size(msk)) ~= size(dat,1)
			dat(msk,:) = [];
	
			k = size(dat,1);
			% A priori probability data for eg WM GM CSF scalp etc..
			bp = zeros(d(1)*d(2),nb+1);
			for j=1:nb
				M = inv(B*reshape(Mat(:,j),4,4));
				tmp = spm_slice_vol(VB(:,j),M, d(1:2),1)/nc(j);
				bp(:,j) = tmp(:);
			end
			bp(msk,:) = [];
			% Other tissue
			bp(:,nb+1) = ...
				abs(ones(k,1) -...
				bp(:,1:nb)*nc(1:nb)')/nc(nb+1);
	
			if (iter==1)
				% Initial probability estimates based upon
				% a priori knowledge
				for i=1:n
					pr(:,i) = bp(:,lkp(i));
					sumbp(i) = sumbp(i) + sum(bp(:,lkp(i)));
				end
			else
				% Compute PDFs for each cluster
				for i=1:n
					c = reshape(cv(:,i),m,m);
					amp = 1/sqrt((2*pi)^m * det(c)) *mg(i);
					dst = (dat-ones(k,1)*mn(:,i)')/sqrtm(c);
					dst = dst.*dst;
					if (size(dst,2)>1)
						dst = sum(dst')';
					end
					pr(:,i)=amp*exp(-0.5*dst).*bp(:,lkp(i));
				end
			end

			% For computing log likelihood
			sumpr = sumpr + sum(log(sum(pr')));

			% PDFs must sum to 1
			sp = sum(pr')';
			msk = find(~sp);
			sp(msk) = ones(size(msk));

			for i=1:n
				pr(:,i) = pr(:,i)./sp;
			end

			% Compute new n, mean & var for each cluster - step 1
			for i=1:n
				mom0(:,i) = mom0(:,i) + sum(pr(:,i));
				mom1(:,i) = mom1(:,i) + ...
					sum((pr(:,i)*ones(1,m)).*dat)';
				mom2(:,i) = mom2(:,i) + reshape(...
					((pr(:,i)*ones(1,m)).*dat)'*dat,m*m,1);
			end
		end
	end

	if iter>1, spm_chi2_plot('Set',sumpr); end

	% Compute new n, mean & var for each cluster - step 2
	for i=1:n
		mg(:,i) = mom0(:,i)/sumbp(i);
		if any(opts == 'f') & i<=nb , mg(:,i) = 1; end;
		mn(:,i) = mom1(:,i)/mom0(:,i);
		tmp = (mom0(:,i).*mn(:,i))*mn(:,i)';

		%stabilise the covariance matrix
		tmp = tmp - eye(size(tmp))*eps*100;

		cv(:,i) = (0.000001+mom2(:,i) - tmp(:))/mom0(:,i);
	end

	if iter==1
		nn = 0;
		for j=1:length(nc)
			% Split the clusters
			for i=2:nc(j)
				cv(:,nn+i) = cv(:,nn+i)*0.8^(1-i);
				mn(:,nn+i) = mn(:,nn+i)*0.8^(1-i);
			end
			nn = nn + nc(j);
		end
	end


	% Background Cluster.
	%    Strictly speaking, since voxels contain absolute values,
	%    the distributions should be modified accordingly. However
	%    in practice, the only allowance made for this is for the
	%    distribution of a background cluster. The mean of this
	%    cluster is assumed to be zero, and in order for the model
	%    to fit properly, the number of voxels contained in this 
	%    cluster is doubled.
	mn(:,n)   = zeros(m,1);
	mg(:,n)   = 2*mg(:,n);
	cv(:,n)   = (0.000001+mom2(:,n))/mom0(:,n);

	% Stopping criterion
	if iter == 2
		sumpr2 = sumpr;
	elseif iter > 3
		if (sumpr-osumpr)/(sumpr-sumpr2) < 0.0003
			break;
		end
	end
	osumpr = sumpr;

end
global loglik
loglik = spm_chi2_plot('Data');
spm_chi2_plot('Clear');



%-----------------------------------------------------------------------

% Create headers, open files etc...
%-----------------------------------------------------------------------
dm     = VF(1:3);
vx     = VF(4:6);
M      = MF1;
orgn   = M\[0 0 0 1]';
orgn   = round(orgn(1:3)');
planes = 1:VF(3);

k = prod(dm(1:2));

if any(opts == 't')
	app  = '_seg_tmp';
	nimg = 2;
else
	app  = '_seg';
	nimg = nb+1;
end

for j=1:nimg
	iname = [spm_str_manip(PF(1,:),'rd') app num2str(j)];
	fp(j) = fopen([iname '.img'],'w');
	if (fp(j) == -1)
		open_error_message([iname '.img']);
		error(['Failed to open ' iname '.img']);
	end
	spm_hwrite([iname '.img'],dm,vx,1/255,2,0,orgn,...
		'Segmented image');
	spm_get_space(iname,M);
end


% Write the segmented images.
%-----------------------------------------------------------------------
spm_progress_bar('Init',dm(3),'Writing Segmented','planes completed');
disp('Writing Segmented');
clear pr bp dat

for pp=1:size(planes,2)
	p = planes(pp);
	B  = spm_matrix([0 0 -p]);
	M2 = B*inv(MM*MF1);

	for i=1:m
		% The image data
		M1 = inv(B*(MF1\reshape(MF(:,i),4,4)));
		tmp = spm_slice_vol(VF(:,i), M1, dm(1:2), 1);

		% Routine is designed for images without -ve values.
		dat(:,i) = abs(tmp(:));
	end

	bp = zeros(size(dat,1),nb+1);
	for j=1:nb
		M = inv(M2*reshape(Mat(:,j),4,4));
		tmp = spm_slice_vol(VB(:,j), M, dm(1:2),1);
		bp(:,j) = tmp(:)/nc(j);
	end

	bp(:,nb+1) = ...
		abs(ones(k,1) -...
		bp(:,1:nb)*nc(1:nb)')/nc(nb+1);

	for i=1:n
		c = reshape(cv(:,i),m,m);
		amp = 1/sqrt((2*pi)^m * det(c))*mg(i);
		dst = (dat-ones(k,1)*mn(:,i)')/sqrtm(c);
		dst = dst.*dst;
		if (size(dst,2)>1)
			dst = sum(dst')';
		end
		pr(:,i) = amp*exp(-0.5*dst) .* bp(:,lkp(i));
	end

	sp = (sum(pr')'+eps)/255;

	for j=1:nimg
		tmp = find(lkp(1:(length(lkp)-1))==j);
		if length(tmp) == 1
			dat = pr(:,tmp);
		else
			dat = sum(pr(:,tmp)')';
		end
		dat = round(dat./sp);
		if fwrite(fp(j),dat,'uchar') ~= prod(size(dat))
			write_error_message([spm_str_manip(PF(1,:),'rd') app num2str(j)]);
			error(['Failed to write ' [spm_str_manip(PF(1,:),'rd') app num2str(j)]]);
		end
	end

	spm_progress_bar('Set',pp);
end
spm_progress_bar('Clear');

for j=1:nimg
	fclose(fp(j));
end

for v=[VF VB]
	spm_unmap(v);
end


% Do the graphics
%=======================================================================
spm_figure('Clear','Graphics');
figure(spm_figure('FindWin','Graphics'));

% Show some text
%-----------------------------------------------------------------------
axes('Position',[0.05 0.8 0.9 0.2],'Visible','off');
text(0.5,0.80, 'Segmentation','FontSize',16,'FontWeight','Bold',...
	'HorizontalAlignment','center');

text(0,0.65, ['Image:  ' spm_str_manip(PF(1,:),'k60d')],...
	'FontSize',14,'FontWeight','Bold');

text(0,0.40, 'Means:','FontSize',12,'FontWeight','Bold');
text(0,0.30, 'Std devs:' ,'FontSize',12,'FontWeight','Bold');
text(0,0.20, 'N vox:','FontSize',12,'FontWeight','Bold');
for j=1:nb
	text((j+0.5)/(nb+2),0.40, num2str(mn(1,j)),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center');
	text((j+0.5)/(nb+2),0.30, num2str(sqrt(cv(1,j))),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center');
	text((j+0.5)/(nb+2),0.20, num2str(mg(1,j)*sumbp(j)),...
		'FontSize',12,'FontWeight','Bold',...
		'HorizontalAlignment','center');
end
if m > 1
	text(0,0.10, ...
	'Note: only means and variances for the first image are shown','FontSize',12);
end

% and display a few images.
%-----------------------------------------------------------------------
V = spm_map(deblank(PF(1,:)));
for j=1:nimg
	iname = [spm_str_manip(PF(1,:),'rd') app num2str(j) '.img'];
	VS(:,j) = spm_map(iname);
end
M1 = spm_get_space(iname);
M2 = spm_get_space(PF(1,:));
for i=1:5
	M = spm_matrix([0 0 i*V(3)/6]);
	img = spm_slice_vol(V,M,V(1:2),0);
	img(1,1) = eps;
	axes('Position',[0.05 0.75*(1-i/5)+0.05 0.9/(nb+2) 0.75/5],'Visible','off');
	imagesc(rot90(img));
	axis('off','image');

	for j=1:nimg
		img = spm_slice_vol(VS(:,j),M2\M1*M,V(1:2),0);
		axes('Position',...
			[0.05+j*0.9/(nb+2) 0.75*(1-i/5)+0.05 0.9/(nb+2) 0.75/5],...
			'Visible','off');
		image(rot90(img*64));
		axis('off','image');
	end
end
spm_unmap(V);
for i=1:size(VS,2)
	spm_unmap(VS(:,i));
end
clear VS;

spm_print;
drawnow;


if any(opts == 't')
	fprintf('Smoothing.\n');
	for i=1:nimg
		iname1 = [spm_str_manip(PF(1,:),'rd') app num2str(i)];
		iname2 = [spm_str_manip(PF(1,:),'rd') '_sseg_tmp' num2str(i)];
		spm_smooth([iname1 '.img'],[iname2 '.img'],8);
		spm_unlink([iname1 '.img'], [iname1 '.hdr'], [iname1 '.mat']);
	end
end
return;




function open_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;

function write_error_message(q)
f=spm_figure('findwin','Graphics'); 
if ~isempty(f), 
	figure(f); 
	spm_figure('Clear','Graphics'); 
	spm_figure('Clear','Interactive'); 
	ax=axes('Visible','off','Parent',f); 
	text(0,0.60,'Error opening:', 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.55,spm_str_manip(q,'k40d'), 'FontSize', 25, 'Interpreter', 'none'); 
	text(0,0.40,'  Please check that you have write permission.', 'FontSize', 16, 'Interpreter', 'none'); 
end
return;

