function spm_realign(arg1,arg2,arg3,arg4,arg5,arg6)
% within mode image realignment routine
%
% --- The Prompts Explained ---
%
% 'number of subjects'
% Enter the number of subjects you wish to realign.
%
% For fMRI, it will ask you the number of sessions for each subject.
% In the coregistration step, the sessions are first realigned to
% each other, by aligning the first scan from each session to the
% first scan of the first session.  Then the images within each session
% are aligned to the first image of the session.
% The parameter estimation is performed this way because it is assumed
% (rightly or not) that there may be systematic differences
% in the images between sessions.
% The adjustment step (correcting for resampling artifacts) is also
% performed completely independantly between each of the fMRI sessions.
%
% 'select scans for subject ..'
% Select the scans you wish to realign. All operations are relative
% to the first image selected.
%
% ......... Note that not all of the following prompts may be used: .........
%
% 'Which option?'
% 	'Coregister only'
% 	Only determine the parameters required to transform each of the
% 	images 2..n to the same space as image 1.
% 	The determined parameters for image XXXX.img are saved in the
%	file XXXX.mat. This is a Matlab file, containing the matrix 'M'.
% 	The location of an image voxel (in mm) can be determined by
% 	computing M(1:3,:)*[xcoord ycoord zcoord 1].
%	Note that if the coregistration is performed more than once on
%	the unresliced data, the starting estimates are obtained from
% 	parameters stored in the '.mat' files.
%	Note that for PET, the coregistration is a two step process.
%	First of all, the images are all realigned to the first in
%	the series.  A mean of these realigned images is created, and
%	a second pass realignment is performed to realign all the
%	images to the mean. Finally, the parameters are corrected
%	for any differences estimated by registering the first image in
%	the series to the mean image.
%
% 	'Reslice Only'
% 	Reslice the specified images according to the contents of the
% 	previously determined parameters. The images are resliced to be
% 	in the same space as the first one selected.  For fMRI, this is
%	the first image of the first session.
%
% 	'Coregister & Reslice'
% 	Combine the above two steps together.
%
% 'Mask brain when registering?'
% The object is to match brains together using a rigid body transformation.
% It is possible that nonrigid movement of structures outside the brain
% (eyeballs, toung, cheeks etc) may interfere with this registration.
% To get around this problem, non-brain areas of the image can be masked
% out during the parameter estimation step.  In order to acheive this, a
% template (of the appropriate modality) is first registered with the images.
% This allows a brain/non-brain mask to be overlayed on the images during
% the registration process.
%
%
% Options for reslicing:
%
% 'Create what?'
% 	'All Images (1..n)'
% 	This reslices all the images - including the first image selected
% 	- which will remain in it's original position.
%
%	'Images 2..n'
% 	Reslices images 2..n only. Useful for if you wish to reslice
% 	(for example) a PET image to fit a structural MRI, without
% 	creating a second identical MRI volume.
%
%	'All Images + Mean Image'
% 	In addition to reslicing the images, it also creates a mean of the
% 	resliced image.
%
%	'Mean Image Only'
% 	Creates the mean image only.
%
% 'Mask the images?'
% To avoid artifactual movement-related variance.
% Because of subject motion, different images are likely to have different
% patterns of zeros from where it was not possible to sample data.
% With masking enabled, the program searches through the whole time series
% looking for voxels which need to be sampled from outside the original
% images. Where this occurs, that voxel is set to zero for the whole set
% of images.
%
% 'Adjust sampling errors?' (fMRI only)
% Adjust the data (fMRI) to remove interpolation errors arising from the
% reslicing of the data.  The adjustment for each fMRI session is performed
% independantly of any other session.
%
% 'Reslice Interpolation Method?'
% 	'Trilinear Interpolation'
% 	Use trilinear interpolation (first-order hold) to sample the images
%       during the writing of realigned images.
%
% 	'Sinc Interpolation'
% 	Use a sinc interpolation to sample the images during the writing
%	of realigned images.
% 	This is slower than bilinear interpolation, but produces better
% 	results. It is especially recommended for fMRI time series.
%
%__________________________________________________________________________
% Refs:
%
% Friston KJ, Ashburner J, Frith CD, Poline J-B, Heather JD & Frackowiak
% RSJ (1995) Spatial registration and normalization of images Hum. Brain
% Map. 2:165-189
%
% Friston KJ, Williams SR, Howard R Frackowiak RSJ and Turner R (1995)
% Movement-related effect in fMRI time-series.  Mag. Res. Med. 35:346-355
%
%__________________________________________________________________________
% %W% Karl Friston & John Ashburner %E%



global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn SWD

if (nargin == 0)
	% User interface.
	%_______________________________________________________________________
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Realignment')
	spm_help('!ContextHelp','spm_realign.m');

	pos = 1;

	n     = spm_input('number of subjects', pos, 'e', 1);
	if (n < 1)
		spm_figure('Clear','Interactive');
		return;
	end

	P = cell(n,1);
	sessions = cell(n,1);

	pos = pos + 1;
	for i = 1:n
		if strcmp(MODALITY,'FMRI'),
			ns = spm_input(['num sessions for subject 'num2str(i)], pos, 'e', 1);
			sess = zeros(1,ns);
			pp = [];
			for s=1:ns
				p = '';
				while size(p,1)<1,
					p = spm_get(Inf,'.img',...
						['scans for subj ' num2str(i) ', sess' num2str(s)]);
				end;
				if s==1,
					pp = p;
				else
					pp=str2mat(pp,p);
				end
				sess(s) =  size(p,1);
			end
			sessions{i} = cumsum(sess);
			P{i} = pp;
		else
			ns = 1;
			P{i} = spm_get(Inf,'.img',...
				['select scans for subject ' num2str(i)]);
			sessions{i} = size(P{i},1);
		end
	end

	FlagsC = ' ';
	FlagsR = ' ';

	if sptl_WhchPtn == 1
		WhchPtn = 3;
	else
		WhchPtn = spm_input('Which option?', pos, 'm',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3],3);
		pos = pos + 1;
	end

	% Co-registration options
	%-----------------------------------------------------------------------
	if (WhchPtn == 1 | WhchPtn == 3)
		Q = '';
		tmp = spm_input('Mask brain when registering?', pos, 'm',...
			'Mask Brain|Dont Mask Brain',...
			[1 0],2);
		if tmp==1
			if (strcmp(MODALITY, 'FMRI')),
				def = 1;
			else
				def = 2;
			end
			DIR1 = [SWD '/coreg/'];
			templates = str2mat([DIR1 'EPI.img'],[DIR1 'PET.img'], ...
				[DIR1 'T1.img'], [DIR1 'T2.img'],[DIR1 'Transm.img']);
			tmp = spm_input('Modality of images?',3,'m',...
				'EPI MR images|PET images|T1 MR images|T2 MR images|Transm images',...
				[1 2 3 4 5],def);
			Q = str2mat(deblank(templates(tmp,:)),[SWD '/apriori/brainmask.img']);
		end
		pos = pos + 1;
	end

	% Reslicing options
	%-----------------------------------------------------------------------
	if (WhchPtn == 2 | WhchPtn == 3),
		p = spm_input('Reslice interpolation method?',pos,'m',...
			'Trilinear Interpolation|Sinc Interpolation',[1 2],2);
		pos = pos + 1;
		if (p == 2) FlagsR = [FlagsR 'S']; end

		if sptl_CrtWht == 1
			p = 3;
		else
			p = spm_input('Create what?',pos,'m',...
				[' All Images (1..n)| Images 2..n|'...
				 ' All Images + Mean Image| Mean Image Only'],...
				[1 2 3 4],3);
			pos = pos + 1;
		end
		if (p == 2) FlagsR = [FlagsR 'n']; end
		if (p == 3) FlagsR = [FlagsR 'i']; end
		if (p == 4) FlagsR = [FlagsR 'Ni']; end
		if (~any(FlagsR == 'N'))
			if sptl_MskOptn == 1
				FlagsR = [FlagsR 'k'];
			else
				if (spm_input('Mask the resliced images?',pos,'y/n') == 'y')
					FlagsR = [FlagsR 'k'];
				end
				pos = pos + 1;
			end
			if (strcmp(MODALITY, 'FMRI'))
				if (sptl_DjstFMRI == 1)
					FlagsR = [FlagsR 'c'];
				elseif sptl_DjstFMRI ~= 0
					if (spm_input(...
						'Adjust sampling errors?'...
						,pos,'y/n') == 'y')
						FlagsR = [FlagsR 'c'];
					end
					pos = pos + 1;
				end
			end
		end
	end

	set(spm_figure('FindWin','Interactive'),'Name','executing',...
		'Pointer','Watch'); drawnow

	for i = 1:n
		fprintf('\n---- Subject %d -----\n', i);
		if WhchPtn==1 | WhchPtn==3,
			realign_images(P{i},Q,sessions{i});
		end
		if WhchPtn==2 | WhchPtn==3,
			reslice_images(P{i},FlagsR,sessions{i})
		end;
	end
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','',...
		'Pointer','Arrow');
	drawnow
	return;

elseif nargin == 1 & strcmp(arg1,'Defaults'),
	edit_defaults;
	return;
elseif nargin >= 1 & strcmp(arg1,'Reslice'),
	P = arg2;
	Flags = 'n';
	Sessions = size(P,1);
	if nargin>=3,
		Flags = arg3;
		if nargin>=4,
			Sessions = arg4;
		end;
	end;
	reslice_images(P,Flags,Sessions);
end;

%_______________________________________________________________________
%_______________________________________________________________________
%_______________________________________________________________________
function reslice_images(P,Flags,sessions)
% FORMAT reslice_images(P,Flags,sessions)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% sessions - the last scan in each of the sessions.  For example,
%            the images in the second session would be
%            P((sessions(2-1)+1):sessions(2),:).
%
% Flags - options flags
%
%         c - adjust the data (fMRI) to remove movement-related components.
%
%         k - mask output images
%             To avoid artifactual movement-related variance the realigned
%             set of images can be internally masked, within the set (i.e.
%             if any image has a zero value at a voxel than all images have
%             zero values at that voxel).  Zero values occur when regions
%             'outside' the image are moved 'inside' the image during
%             realignment.
%
%         i - write mean image
%             The average of all the realigned scans is written to
%             mean*.img.
%
%         S - use sinc interpolation for reslicing (11x11x11).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be
%             necessary to resample it.
%
%         N - don't reslice any of the images - except possibly create a
%             mean image.
%
%             The spatially realigned and adjusted images are written to
%             the orginal subdirectory with the same filename but prefixed
%             with a 'r'.
%             They are all aligned with the first.

P=spm_vol(P);

fprintf('Reslicing images..\n');
Hold = 1;
if (any(Flags == 'S')) Hold = -11; end

spm_progress_bar('Init',P(1).dim(3),'Reslicing','planes completed');


start_vol = 1;
if (any(Flags == 'n') & ~any(Flags == 'i')) start_vol = 2; end

if any(Flags == 'i')
	Integral = zeros(P(1).dim(1)*P(1).dim(2),P(1).dim(3));
end

tiny = 5e-2; % From spm_vol_utils.c
IC0  = eye(6)*eps*100;
open_mode = 'w';

if any(Flags == 'c'),
	Y1 = zeros(prod(P(1).dim(1:2)),prod(size(P)));
	Y2 = zeros(prod(P(1).dim(1:2)),prod(size(P)));
	Y3 = zeros(prod(P(1).dim(1:2)),prod(size(P)));
end;
X  = zeros(prod(P(1).dim(1:2)),prod(size(P)));


x1=repmat((1:P(1).dim(1))',1,P(1).dim(2));x1=x1(:);
x2=repmat( 1:P(1).dim(2)  ,P(1).dim(1),1);x2=x2(:);

for x3 = 1:P(1).dim(3)
	Count = zeros(prod(P(1).dim(1:2)),1);
	for i = 1:prod(size(P))
		M = inv(P(1).mat\P(i).mat);
		y1=M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
		y2=M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
		y3=M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));

		Mask =        (y1 >= (1-tiny) & y1 <= (P(i).dim(1)+tiny));
		Mask = Mask & (y2 >= (1-tiny) & y2 <= (P(i).dim(2)+tiny));
		Mask = Mask & (y3 >= (1-tiny) & y3 <= (P(i).dim(3)+tiny));
		Count = Count + Mask;

		d  = spm_sample_vol(P(i),y1,y2,y3,Hold);
		X(:,i)= d;

		if any(Flags == 'c'),
			Y1(:,i)=y1;
			Y2(:,i)=y2;
			Y3(:,i)=y3;
		end;
	end
	Mask = (Count == prod(size(P)));

	if any(Flags == 'c'),
		ss = 1;
		for s = 1:length(sessions),
			se = sessions(s);
			for xx=find(Mask)',
				% Basis functions appropriate for sinc interpolation
				A = [	cos(2*pi*Y1(xx,ss:se)') sin(2*pi*Y1(xx,ss:se)') ...
			  		cos(2*pi*Y2(xx,ss:se)') sin(2*pi*Y2(xx,ss:se)') ...
			   		cos(2*pi*Y3(xx,ss:se)') sin(2*pi*Y3(xx,ss:se)')];

				Tac = X(xx,ss:se)';
				% centre the signal and covariates
				Tac = Tac - mean(Tac);
				A   = A - repmat(mean(A,1),size(A,1),1);
				Tac = X(xx,ss:se)' - A*((A'*A + IC0)\(A'*Tac));
				X(xx,ss:se) = Tac';
			end
			ss = se+1;
		end;
	end;


	if any(Flags == 'i')
		Integral(:,x3) = sum(X,2)./Count;
	end

	if any(Flags == 'k'),
		X(find(~Mask),:) = 0;
	end

	if ~any(Flags == 'N'),
		start_vol = 1;
		if (any(Flags == 'n')) start_vol = 2; end

		for i = start_vol:prod(size(P)),
			% May need to change this later if more image formats are
			% supported.
			d = X(:,i)/P(i).pinfo(1,1);
			if any(P(i).dim(4) == [2 4 8]), d = round(d); end
			if P(i).dim(4) == 2,
				d(find(d <   0)) =   0;
				d(find(d > 255)) = 255;
			elseif (P(i).dim(4) == 4)
				d(find(d >  32767)) =  32767;
				d(find(d < -32768)) = -32768;
			elseif (P(i).dim(4) == 8)
				d(find(d >  2^31-1)) =  2^31-1;
				d(find(d < -2^31  )) = -2^31;
			end
			p = deblank(P(i).fname);
			q      = max([find(p == '/') 0]);
			q      = [p(1:q) 'r' p((q + 1):length(p))];

			fp = fopen(q, open_mode);
			if (fp ~= -1),
				if fwrite(fp,d,spm_type(P(i).dim(4))) ~= prod(size(d)),
					fclose(fp);
					write_error_message(q);
					error(['Error writing ' q '. Check your disk space.']);
				end;
				fclose(fp);
			else
				open_error_message(q);
				error(['Error opening ' q '. Check that you have write permission.']);
			end
		end
	end
	open_mode = 'a';
	spm_progress_bar('Set',x3);
end

% Write headers and matrixes
%------------------------------------------------------------------
if ~any(Flags == 'N'),
	vox    = sqrt(sum(P(1).mat(1:3,1:3).^2));
	origin = round(P(1).mat\[0 0 0 1]');
	origin = origin(1:3);
	dim    = P(1).dim(1:3);
	for i = start_vol:prod(size(P)),
		dtype  = P(i).dim(4);
		scale  = P(i).pinfo(1,1);
		p = deblank(P(i).fname);
		q = max([find(p == '/') 0]);
		q = [p(1:q) 'r' p((q + 1):length(p))];
		spm_hwrite(q,dim,vox,scale,dtype,0,origin,'spm - realigned');
		spm_get_space(q,P(1).mat);
	end
end

if any(Flags == 'i')
	% Write integral image (16 bit signed)
	%-----------------------------------------------------------
	mx = max(max(Integral));
	scale  = mx/32767;
	p      = P(1).fname;
	q      = max([find(p == '/') 0]);
	q      = [p(1:q) 'mean' p((q + 1):length(p))];
	fp = fopen(q,'w');
	if (fp ~= -1)
		vox    = sqrt(sum(P(1).mat(1:3,1:3).^2));
		origin = round(P(1).mat\[0 0 0 1]');
		origin = origin(1:3);
		dim    = P(1).dim(1:3);
		for j = 1:P(1).dim(3)
			d = round(Integral(:,j)/scale);
			d(find(d < -32768)) = -32768;
			if fwrite(fp,d,spm_type(4)) ~= prod(size(d))
				fclose(fp);
				write_error_message(q);
				error(['Error writing ' q '. Check your disk space.']);
			end
		end
		spm_hwrite(q,dim,vox,scale,4,0,origin,'spm - mean image');
		spm_get_space(q,P(1).mat);
	else
		open_error_message(q);
		error(['Error opening ' q '. Check that you have write permission.']);
	end
end

spm_figure('Clear','Interactive');
fprintf('Done\n');
return;
%_______________________________________________________________________


%_______________________________________________________________________
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
%_______________________________________________________________________

%_______________________________________________________________________
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
%_______________________________________________________________________



%_______________________________________________________________________
function realign_images(P,Q,sessions)
% FORMAT realign_images(P,Q,sessions)
%
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% sessions - the last scan in each of the sessions.  For example,
%            the images in the second session would be
%            P((sessions(2-1)+1):sessions(2),:).
%
% Q     - an optional matrix od filenames.  These are used for masking
%         out regions which are (roughly) considered to be outside the
%         brain.  The last filename is an image containing values
%         between zero and one, where each value is a weight.  The
%         other images are template images.  The way that this works is
%         that the first image in the series is matched (using an affine
%         transformation) to a linear combinantion of the template images.
%         This affine mapping can then be used to overlay the weight image
%         over the first image of the series.
%


global MODALITY
P=spm_vol(P);
Flags = ' ';
if strcmp(MODALITY, 'PET'), Flags = [Flags 'p']; end;
if nargin<2, Q = ''; end
if ~isempty(Q),
	tmp = clock;
	ref = [pwd '/spm_realign_tmp_' sprintf('%.2d%.2d%.2d%.2d%.4d',tmp(3),tmp(4),tmp(5),round(tmp(6)),round(rem(tmp(6),1)*10000)) '.img'];
	spm_smooth(P(1).fname,ref,8);
	fprintf('Initial registration to templates..\n');
	params  =spm_affsub3('affine2',Q(1:(end-1),:),ref,1,8);
	params = spm_affsub3('affine2',Q(1:(end-1),:),ref,1,6,params);
	tmp = spm_str_manip(ref,'rd');
	spm_unlink(ref,[tmp '.hdr'],[tmp '.mat']);
	MW = inv(spm_matrix(params(1:12)));
	PW = spm_vol(Q(end,:));
	PW.mat = inv(spm_matrix(params(1:12)))*PW.mat;
else
	PW=[];
end

if (length(sessions)==1),
	fprintf('Registering images..\n');
	P = realign_series(P,PW,Flags);
else
	fprintf('Registering together the first image of each session..\n');
	tmp = [1 sessions(1:(end-1))-1];
	Ptmp = realign_series(P(tmp),PW,Flags);
	ss=1;
	for s=1:length(sessions)
		es = sessions(s);
		M = Ptmp(s).mat*inv(P(ss).mat);
		for i=ss:es
			P(i).mat = M*P(i).mat;
		end;
		ss = es+1;
	end;

	ss=1;
	for s=1:length(sessions)
		es = sessions(s);
		fprintf('Registering together images from session %d..\n', s);
		P(ss:es) = realign_series(P(ss:es),PW,Flags);
		ss = es+1;
	end;
end


plot_parameters(P);

% Save Realignment Parameters
%---------------------------------------------------------------------------
fprintf('Saving parameters..\n');
for i=1:prod(size(P)),
	spm_get_space(P(i).fname, P(i).mat);
end;
return;
%_______________________________________________________________________


%_______________________________________________________________________
function P = realign_series(P,Wt,Flags)
% Realign a time series of 3D images to the first of the series.
% FORMAT P = spm_realign_series(P,Wt,Flags)
% P  - a vector of volumes (see spm_vol)
% Wt - an optional masking volume
%-----------------------------------------------------------------------
% P(i).mat is modified to reflect the modified position of the image i.
% The scaling (and offset) parameters are also set to contain the
% optimum scaling required to match the images.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin < 2, Wt = []; end;
if nargin < 3, Flags = ' '; end;

Hold1 = -5;	% Interpolation method
fwhm  = 6;	% Smoothing kernel size
separation = 5;	% Spacing (in millimeters) between samples
register_to_mean = 0;

% because PET images are noisier:
if any(Flags == 'p'), register_to_mean = 1; fwhm = 8; end;

if any(Flags == 'd'), separation = 3; end;
if any(Flags == 's'), fwhm = 10; end;

lkp = [1 2 3 4 5 6];
if P(1).dim(3) < 3, lkp = [1 2 6]; end;

% Points to sample in reference image
%-----------------------------------------------------------------------
skip = sqrt(sum(P(1).mat(1:3,1:3).^2)).^(-1)*separation;
d    = P(1).dim(1:3);
[x1,x2,x3]=meshgrid(1:skip(1):d(1),1:skip(2):d(2),1:skip(3):d(3));
x1 = x1(:);
x2 = x2(:);
x3 = x3(:);


% Possibly mask an area of the sample volume.
%-----------------------------------------------------------------------
if ~isempty(Wt),
	[y1,y2,y3]=coords([0 0 0  0 0 0],P(1).mat,Wt.mat,x1,x2,x3);
	wt = spm_sample_vol(Wt,y1,y2,y3,1);
	msk = find(wt>0.01);
	x1=x1(msk);
	x2=x2(msk);
	x3=x3(msk);
	wt=wt(msk);
else
	wt = [];
end
n=prod(size(x1));


% Compute rate of change of chi2 w.r.t changes in parameters (matrix A)
% Also compute A'A for the whole volume, so that A'A for the masked
% volume can be computed more efficiently. 
%-----------------------------------------------------------------------
V=smooth_vol(P(1),fwhm);
[b,dG1,dG2,dG3]=spm_sample_vol(V,x1,x2,x3,Hold1);
clear V
A0 = make_A(P(1).mat,x1,x2,x3,dG1,dG2,dG3,wt,lkp);
if ~isempty(wt), b = b.*wt; end

if register_to_mean,
	count = ones(n,1);
	ave   = b;
	grad1 = dG1;
	grad2 = dG2;
	grad3 = dG3;
end;

spm_progress_bar('Init',length(P)-1,'Registering Images');
% Loop over images
%-----------------------------------------------------------------------
for i=2:length(P)
	V=smooth_vol(P(i),fwhm);
	ss = Inf;
	countdown = -1;
	Hold = 1;	% Begin with bi-linear interpolation.
	for iter=1:64
		[y1,y2,y3] = coords([0 0 0  0 0 0],P(1).mat,P(i).mat,x1,x2,x3);
		msk        = find((y1>=1 & y1<=d(1) & y2>=1 & y2<=d(2) & y3>=1 & y3<=d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_sample_vol(V, y1(msk),y2(msk),y3(msk),Hold);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = [A0(msk,:) F];
		Alpha      = spm_atranspa(A);
		Beta       = A'*b(msk);
		soln       = Alpha\Beta;

		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = soln(1:(end-1));
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum((F*soln(end)-b(msk)).^2)/length(msk);
		if (pss-ss)/pss < 1e-7	% Stopped converging.
			if Hold == 1
				% Switch to a better (slower) interpolation
				% to finish with
				Hold = Hold1;
			elseif countdown == -1
				% Do three final iterations to finish off with
				countdown = 3;
			end
		end
		if countdown ~= -1
			if countdown==0, break; end;
			countdown = countdown -1;
		end

% figure(2);
% img1 = spm_slice_vol(P(i),inv(spm_matrix([0 0 -25])*inv(P(1).mat)*P(i).mat),P(1).dim(1:2),1);
% img2 = spm_slice_vol(P(1),inv(spm_matrix([0 0 -25])),P(1).dim(1:2),1);
% img3 = spm_slice_vol(Wt,inv(spm_matrix([0 0 -25])*inv(P(1).mat)*Wt.mat),P(1).dim(1:2),1);
% subplot(2,2,1);imagesc(rot90(img2));axis image off;
% subplot(2,2,2);imagesc(rot90(img1));axis image off;
% subplot(2,2,3);imagesc(rot90(img3));axis image off;
% drawnow;
	end

	if register_to_mean,
		% Generate mean and derivatives of mean
		count(msk) = count(msk) + 1;
		[G,dG1,dG2,dG3] = spm_sample_vol(V,y1(msk),y2(msk),y3(msk),Hold1);
		ave(msk)   = ave(msk)   + G.*soln(end);
		grad1(msk) = grad1(msk) + dG1.*soln(end);
		grad2(msk) = grad2(msk) + dG2.*soln(end);
		grad3(msk) = grad3(msk) + dG3.*soln(end);
	end;
	spm_progress_bar('Set',i-1);
end;
spm_progress_bar('Clear');

if ~register_to_mean, return; end;
%_______________________________________________________________________
M=P(1).mat;
A0 = make_A(M,x1,x2,x3,grad1./count,grad2./count,grad3./count,wt,lkp);
if ~isempty(wt), b = (ave./count).*wt;
else, b = (ave./count); end

clear ave grad1 grad2 grad3

% Loop over images
%-----------------------------------------------------------------------
spm_progress_bar('Init',length(P),'Registering Images to Mean');
for i=1:length(P)
	V=smooth_vol(P(i),fwhm);
	ss = Inf;
	countdown = -1;
	for iter=1:64
		[y1,y2,y3] = coords([0 0 0  0 0 0],M,P(i).mat,x1,x2,x3);
		msk        = find((y1>1 & y1<d(1) & y2>1 & y2<d(2) & y3>1 & y3<d(3)));
		if length(msk)<32, error_message(P(i)); end;

		F          = spm_sample_vol(V, y1(msk),y2(msk),y3(msk),Hold1);
		if ~isempty(wt), F = F.*wt(msk); end;

		A          = [A0(msk,:) F];
		Alpha      = spm_atranspa(A);
		Beta       = A'*b(msk);
		soln       = Alpha\Beta;
		p          = [0 0 0  0 0 0  1 1 1  0 0 0];
		p(lkp)     = soln(1:(end-1));
		P(i).mat   = inv(spm_matrix(p))*P(i).mat;

		pss        = ss;
		ss         = sum((F*soln(end)-b(msk)).^2)/length(msk);
		if (pss-ss)/pss < 1e-7 & countdown == -1 % Stopped converging.
			% Do three final iterations to finish off with
			countdown = 3;
		end
		if countdown ~= -1
			if countdown==0, break; end;
			countdown = countdown -1;
		end
	end
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');


% Since we are supposed to be aligning everything to the first
% image, then we had better do so
%-----------------------------------------------------------------------
M = M/P(1).mat;
for i=1:length(P)
	P(i).mat   = M*P(i).mat;
end

return;
%_______________________________________________________________________


%_______________________________________________________________________
function [y1,y2,y3]=coords(p,M1,M2,x1,x2,x3)
% Rigid body transformation of a set of coordinates.
M  = (inv(M2)*inv(spm_matrix(p(1:6)))*M1);
y1 = M(1,1)*x1 + M(1,2)*x2 + M(1,3)*x3 + M(1,4);
y2 = M(2,1)*x1 + M(2,2)*x2 + M(2,3)*x3 + M(2,4);
y3 = M(3,1)*x1 + M(3,2)*x2 + M(3,3)*x3 + M(3,4);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function V=smooth_vol(P,fwhm)
% Convolve the volume in memory.
s = sqrt(sum(P.mat(1:3,1:3).^2)).^(-1)*(fwhm/sqrt(8*log(2)));
x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2));
y  = exp(-(y).^2/(2*(s(2)).^2));
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

V=zeros(P.dim(1:3));
spm_conv_vol(P,V,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function A = make_A(M,x1,x2,x3,dG1,dG2,dG3,wt,lkp)
p0 = [0 0 0  0 0 0  1 1 1  0 0 0];
At=zeros(length(lkp),prod(size(x1)));
for i=1:length(lkp)
	pt         = p0;
	pt(lkp(i)) = pt(i)+1e-6;
	[y1,y2,y3] = coords(pt,M,M,x1,x2,x3);
	tmp        = sum([y1-x1 y2-x2 y3-x3].*[dG1 dG2 dG3],2)/(-1e-6);
	if ~isempty(wt), A(:,i) = tmp.*wt;
	else, A(:,i) = tmp; end
end
return;
%_______________________________________________________________________

%_______________________________________________________________________
function error_message(P)
f=spm_figure('findwin','Graphics');
if ~isempty(f),
	figure(f);
	spm_figure('Clear','Graphics');
	spm_figure('Clear','Interactive');
	ax=axes('Visible','off','Parent',f);
	text(0,0.60,'There is not enough overlap in the', 'FontSize', 25, 'Interpreter', 'none');
	text(0,0.55,'    images to obtain a solution.', 'FontSize', 25, 'Interpreter', 'none');
	text(0,0.40,'  Please check that your header information is OK.','FontSize', 16, 'Interpreter', 'none');
	text(0,0.25, ['Offending image: "' P.fname '".'],'FontSize', 12, 'Interpreter', 'none');
end

fprintf('There is not enough overlap of the images to obtain a solution\n');
error(['The offending image is "' P.fname '".']);
return;
%_______________________________________________________________________



%_______________________________________________________________________
function plot_parameters(P)
fg=spm_figure('FindWin','Graphics');
if ~isempty(fg)

	Params = zeros(prod(size(P)),12);
	for i=1:prod(size(P))
		Params(i,:) = spm_imatrix(P(1).mat\P(i).mat);
	end

	% display results
	% translation and rotation over time series
	%-------------------------------------------------------------------
	spm_figure('Clear','Graphics');
	ax=axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
	set(get(ax,'Title'),'String','Image realignment','FontSize',16,'FontWeight','Bold','Visible','on');
	x     =  0.1;
	y     =  0.9;
	for i = 1:min([prod(size(P)) 12])
		text(x,y,[sprintf('%-4.0f',i) P(i).fname],'FontSize',10,'Interpreter','none','Parent',ax);
		y = y - 0.08;
	end
	if prod(size(P)) > 12
		text(x,y,'................ etc','FontSize',10,'Parent',ax); end

	ax=axes('Position',[0.1 0.35 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
	plot(Params(:,1:3),'Parent',ax)
	s = ['x translation';'y translation';'z translation'];
	text([2 2 2], Params(2, 1:3), s, 'Fontsize',10,'Parent',ax)
	set(get(ax,'Title'),'String','translation','FontSize',16,'FontWeight','Bold');
	set(get(ax,'Xlabel'),'String','image');
	set(get(ax,'Ylabel'),'String','mm');


	ax=axes('Position',[0.1 0.05 0.8 0.2],'Parent',fg,'XGrid','on','YGrid','on');
	plot(Params(:,4:6)*180/pi,'Parent',ax)
	s = ['pitch';'roll ';'yaw  '];
	text([2 2 2], Params(2, 4:6)*180/pi, s, 'Fontsize',10,'Parent',ax)
	set(get(ax,'Title'),'String','rotation','FontSize',16,'FontWeight','Bold');
	set(get(ax,'Xlabel'),'String','image');
	set(get(ax,'Ylabel'),'String','degrees');

	% print realigment parameters
	spm_print
end
return;
%_______________________________________________________________________


%_______________________________________________________________________
function edit_defaults
tmp = 1;
if sptl_WhchPtn == 1, tmp = 2; end;
sptl_WhchPtn = spm_input(...
	['Coregistration and reslicing?'],...
	2, 'm',...
	['Allow separate coregistration and reslicing|'...
	 'Combine coregistration and reslicing'], [-1 1], tmp);


tmp = 2;
if sptl_CrtWht == 1,
	tmp = 1;
end
sptl_CrtWht   = spm_input(['Images to create?'], 3, 'm',...
	'All Images + Mean Image|Full options', [1 -1], tmp);


tmp = 3;
if sptl_DjstFMRI == 1,
	tmp = 1;
elseif sptl_DjstFMRI == 0
	tmp = 2;
end
sptl_DjstFMRI = spm_input(['fMRI adjustment for interpolation?'],4,'m',...
		'   Always adjust |    Never adjust|Optional adjust',...
		[1 0 -1], tmp);


tmp = 2;
if sptl_MskOptn == 1,
	tmp = 1;
end
sptl_MskOptn  = spm_input(['Option to mask images?'], 5, 'm',...
		'  Always mask|Optional mask', [1 -1], tmp);
return;
