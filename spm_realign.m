function spm_realign(P,Flags)
% within mode image realignment routine
%
% --- The Prompts Explained ---
%
% 'number of subjects'
% Enter the number of subjects you wish to realign.
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
%
% 	'Reslice Only'
% 	Reslice the specified images according to the contents of the
% 	previously determined parameters. The images are resliced to be
% 	in the same space as the first one selected.
%
% 	'Coregister & Reslice'
% 	Combine the above two steps together.
%
% 'Interpolation Method?'
% 	'Bilinear Interpolation'
% 	Use bilinear interpolation (first-order hold) to sample the images
%       when determining parameters/reslicing.
%
% 	'Sinc Interpolation'
% 	Use a sinc interpolation to sample the images.
% 	This is slower than bilinear interpolation, but produces better
% 	results. It is especially recommended for fMRI time series.
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
% reslicing of the data.
%
%____________________________________________________________________________
% TO OBTAIN SIMILAR RESULTS TO SPM95, select:
% 'Which option?'        	'Coregister & Reslice'
% 'Interpolation Method?'	'Bilinear Interpolation'       (for PET)
% 'Interpolation Method?'	'Sinc Interpolation'           (for fMRI)
%
% 'Create what?'        	'All Images + Mean Image'
% 'Mask the images?'       	'YES'
% 'Adjust sampling errors?'	'YES' (for fMRI)
%____________________________________________________________________________
%
%
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
% %W% Karl Friston - modified by John Ashburner %E%

% programmers notes
%
% FORMAT spm_realign(P, Flags)
% P     - matrix of filenames {one string per row}
%         All operations are performed relative to the first image.
%         ie. Coregistration is to the first image, and resampling
%         of images is into the space of the first image.
%
% Flags - options flags
%         e - compute realignment parameters
%             Realigns a time-series of scans whose filenames are
%             in the input matrix (P).
%             The realignment uses least squares solutions and partial
%             derivatives of the image process with respect to the spatial 
%             transformations constituting a rigid-body transformation.
%             These movement parameter estimates are based on a first order
%             Taylor expansion in terms of the partial derivatives measured 
%             using smoothed scans.
%             Subsequent scans are realigned with the first.
%
%         r - reslice images
%             The spatially realigned and adjusted images are written to
%             the orginal subdirectory with the same filename but prefixed
%             with a 'r'.
%             They are all aligned with the first.
%
%         c - adjust the data (fMRI) to remove movement-related components
%             The adjustment procedure is based on a autoregression-moving 
%             average-like model of the effect of position on signal and 
%             explicitly includes a spin excitation history effect.
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
%         s - use sinc interpolation for parameter estimates (7x7x7).
%
%         S - use sinc interpolation for reslicing (9x9x9).
%
%         n - don't reslice the first image
%             The first image is not actually moved, so it may not be
%             necessary to resample it.
%
%         N - don't reslice any of the images - except possibly create a
%             mean image.


global MODALITY sptl_WhchPtn sptl_DjstFMRI sptl_CrtWht sptl_MskOptn

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

	pos = pos + 1;
	for i = 1:n
		P = spm_get(Inf,'.img',...
			['select scans for subject ' num2str(i)]);
		eval(['P' num2str(i) ' = P;']);
	end

	Flags = '';

	if sptl_WhchPtn == 1
		p = 3;
	else
		p = spm_input('Which option?', pos, 'm',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3],3);
		pos = pos + 1;
	end

	if (p == 1 | p == 3) Flags = [Flags 'e']; end
	if (p == 2 | p == 3) Flags = [Flags 'r']; end
	p = spm_input('Interpolation Method?',pos,'m',...
		'Bilinear Interpolation|Sinc Interpolation',[1 2],2);
	pos = pos + 1;
	if (p == 2) Flags = [Flags 'sS']; end

	if (any(Flags == 'r'))
		if sptl_CrtWht == 1
			p = 3;
		else
			p = spm_input('Create what?',pos,'m',...
				[' All Images (1..n)| Images 2..n|'...
				 ' All Images + Mean Image| Mean Image Only'],...
				[1 2 3 4],3);
			pos = pos + 1;
		end
		if (p == 2) Flags = [Flags 'n']; end
		if (p == 3) Flags = [Flags 'i']; end
		if (p == 4) Flags = [Flags 'Ni']; end
		if (~any(Flags == 'N'))
			if sptl_MskOptn == 1
				Flags = [Flags 'k'];
			else
				if (spm_input('Mask the images?'       ,pos,'y/n') == 'y')
					Flags = [Flags 'k'];
				end
				pos = pos + 1;
			end
			if (strcmp(MODALITY, 'FMRI'))
				if (sptl_DjstFMRI == 1)
					Flags = [Flags 'c'];
				elseif sptl_DjstFMRI ~= 0
					if (spm_input(...
						'Adjust sampling errors?'...
						,pos,'y/n') == 'y')
						Flags = [Flags 'c'];
					end
					pos = pos + 1;
				end
			end
		end
	end

	set(spm_figure('FindWin','Interactive'),'Name','executing',...
		'Pointer','Watch'); drawnow
	for i = 1:n
		eval(['P = P' num2str(i) ';'])
		eval('spm_realign(P,Flags);','disp(''Realignment Bombed Out'');');
	end
	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','',...
		'Pointer','Arrow');
	drawnow
	return;

elseif nargin == 1 & strcmp(P,'Defaults')
	% Defaults.
	%_______________________________________________________________________

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
end

% Do the work..
%_______________________________________________________________________

nreg   = size(P,1);
Params = zeros(12,nreg); % initial estimate/default values
Params(7:9,nreg) = 1; % good old Matlab 5

% Computation of Realignment Parameters - uses P & Params
%---------------------------------------------------------------------------
if any(Flags == 'e')

	Hold = 1;
	if (any(Flags == 's')) Hold = -7; end


	dims = spm_hread(P(1,:));

	spm_progress_bar('Init',nreg);

	% Images smoothed prior to realignment to:
	% a) Make the Taylor approximations more valid for larger
	%    displacements so that less iterations are required.
	% b) Less samples are needed to encompass more of the lower
	%    frequency signal.
	% c) It masks the errors due to approximations to sinc
	%    interpolation in the resampling.
	%    More work should be done at looking at the FTs of
	%    the images with respect to FTs of the errors due to
	%    the resampling in order to optimize the smoothness
	%    parameter.
	smo = 8;

	global CWD
	ref = [CWD '/spm_realign_tmpG.img'];
	mov = [CWD '/spm_realign_tmpF.img'];

	spm_smooth(P(1,:),ref,smo);

	for i=2:nreg
		spm_smooth(P(i,:),mov,smo);
		if dims(3)>3,
			% 3D registration
			params=spm_affsub3('rigid1',ref,mov,1,8);
			params=spm_affsub3('rigid1',ref,mov,Hold,6,params);
			%params=spm_affsub3('rigid1',P(1,:)  , P(i,:)  ,Hold,6,params);
		else,
			% 2D registration
			params=spm_affsub3('2d1',ref,mov,1,8);
			params=spm_affsub3('2d1',ref,mov,Hold,4,params);
			%params=spm_affsub3('2d1',P(1,:)  , P(i,:)  ,Hold,4,params);
		end;
		Params(:,i) = params(1:12);
		spm_unlink(mov,[spm_str_manip(mov,'s') '.hdr'],[spm_str_manip(mov,'s') '.mat']);
		spm_progress_bar('Set',i);
	end
	spm_unlink(ref,[spm_str_manip(ref,'s') '.hdr'],[spm_str_manip(ref,'s') '.mat']);
	spm_progress_bar('Clear');

	Params = Params';

	% Saving of Realignment Parameters
	%---------------------------------------------------------------------------
	Matrixes = zeros(16,size(P,1));
	for i=2:size(P,1)
		M = spm_get_space(deblank(P(i,:)));
		Matrixes(:,i) = M(:);
	end

	for i=2:size(P,1)
		M = reshape(Matrixes(:,i),4,4);
		spm_get_space(deblank(P(i,:)), spm_matrix(Params(i,:))*M);
	end


	fg=spm_figure('FindWin','Graphics');
	if ~isempty(fg)
		% display results
		% translation and rotation over time series
		%-------------------------------------------------------------------
		spm_figure('Clear','Graphics');
		ax=axes('Position',[0.1 0.65 0.8 0.2],'Parent',fg,'Visible','off');
		set(get(ax,'Title'),'String','Image realignment','FontSize',16,'FontWeight','Bold','Visible','on');
		x     =  0.1;
		y     =  0.9;
		for i = 1:min([size(P,1) 12])
			text(x,y,[sprintf('%-4.0f',i) P(i,:)],'FontSize',10,'Interpreter','none','Parent',ax);
			y = y - 0.08;
		end
		if size(P,1) > 12
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

		% print realigment parameters and figure 2 attributes
		spm_print
	end
end




% Application of Realignment Parameters - uses P
%---------------------------------------------------------------------------
if any(Flags == 'r')

	Hold = 1;
	if (any(Flags == 'S')) Hold = -9; end

	% Get properties of image to realign to
	%-------------------------------------------------------------------
	p = spm_str_manip(P(1,:), 'd');
	M = spm_get_space(p);
	[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(p);

	% Get properties of all the images
	%-------------------------------------------------------------------
	Headers = zeros(size(P,1),(3+1+1));
	Matrixes = zeros(16,size(P,1));
	V     = zeros(12,size(P,1));
	for i = 1:size(P,1)
		p  = deblank(P(i,:));

		M1  = M\spm_get_space(p);
		Matrixes(:,i) = M1(:);

		[dim dummy scale typ dummy dummy] = spm_hread(p);
		Headers(i,:) = [dim scale typ];

		V(:,i) = spm_map(spm_str_manip(P(i,:),'d'));
	end


	spm_progress_bar('Init',DIM(3),'Reslicing','planes completed');

	Xm = kron(ones(DIM(2),1),[1:DIM(1)]');
	Ym = kron([1:DIM(2)]',ones(DIM(1),1));
	open_mode = 'w';
	if (any(Flags == 'i') | any(Flags == 'c'))
		Integral = zeros(prod(DIM(1:2)),DIM(3));
	end

	start_vol = 1;
	if (any(Flags == 'n') & ~any(Flags == 'i')) start_vol = 2; end

	for j = 1:DIM(3)
		B     = spm_matrix([0 0 -j 0 0 0 1 1 1]);
		if (any(Flags == 'c'))
			X = zeros(DIM(1)*DIM(2),size(P,1));
		end

		% get masks for this plane
		%-----------------------------------------------------------
		if (any(Flags == 'i') | any(Flags == 'k') |...
			any(Flags == 'c'))
			Count = zeros(prod(DIM(1:2)),1);
			for i = 1:size(P,1)
				M0  = reshape(Matrixes(:,i),4,4);
				M1 = inv(B*M0);

				tiny = 5e-2; % From spm_vol_utils.c

				% Find range of slice
				% -----------------------------------------------------------
				tmp = M1(1,1)*Xm + M1(1,2)*Ym + M1(1,4);
				Mask= (tmp >= (1-tiny) & tmp <= (Headers(i,1)+tiny));
				tmp = M1(2,1)*Xm + M1(2,2)*Ym + M1(2,4);
				Mask= Mask & (tmp >= (1-tiny) & tmp <= (Headers(i,2)+tiny));
				tmp = M1(3,1)*Xm + M1(3,2)*Ym + M1(3,4);
				Mask= Mask & (tmp >= (1-tiny) & tmp <= (Headers(i,3)+tiny));

				Count = Count + Mask;
			end
			Mask = (Count == size(P,1));
		end

		% get values for this plane
		%-----------------------------------------------------------
		for i = start_vol:size(P,1)
			M0  = reshape(Matrixes(:,i),4,4);
			M1 = inv(B*M0);

			d  = spm_slice_vol(V(:,i),M1,DIM(1:2),Hold);

			if (any(Flags == 'i') | any(Flags == 'c'))
				Integral(:,j) = Integral(:,j) + d(:);
			end

			if (~any(Flags == 'c'))
				% don't need to load the whole time series.
				if (i > 1 | ~any(Flags == 'n'))
					d = d/Headers(i,4);
					if any(Flags == 'k'), d = d(:).*Mask; end
					if any(Headers(i,5) == [2 4 8]), d = round(d); end
					if Headers(i,5) == 2,
						d(find(d <   0)) =   0;
						d(find(d > 255)) = 255;
					elseif (Headers(i,5) == 4)
						d(find(d >  32767)) =  32767;
						d(find(d < -32768)) = -32768;
					elseif (Headers(i,5) == 8)
						d(find(d >  2^31-1)) =  2^31-1;
						d(find(d < -2^31  )) = -2^31;
					end
				end

				if (~any(Flags == 'N'))
					p  = spm_str_manip(P(i,:), 'd');
					q  = max([find(p == '/') 0]);
					q  = [p(1:q) 'r' p((q + 1):length(p))];
					fp = fopen(q, open_mode);
					if (fp ~= -1)
						if fwrite(fp,d,spm_type(Headers(i,5))) ~= prod(size(d))
							fclose(fp);
							write_error_message(q);
							error(['Error writing ' q '. Check your disk space.']);
						end
						fclose(fp);
					else
						open_error_message(q);
						error(['Error opening ' q '. Check that you have write permission.']);
					end
				end
			else
				X(:,i) = d(:);
			end
		end

		if (any(Flags == 'i') | any(Flags == 'c'))
			Integral(:,j) = Integral(:,j)./(Count + eps);
		end


		if (any(Flags == 'c'))
			% Adjust for motion
			%---------------------------------------------------
			Mask = (Count == size(P,1));

			X = X'; % Transpose to reduce paging

% Covary out perioric function of movement to correct interpolation
% errors.  This needs to be done with a different design matrix for each
% voxel - so it is a little on the slow side.
%-----------------------------------------------------------------------

dXmotion = Matrixes(1,:)';
dYmotion = Matrixes(2,:)';
dZmotion = Matrixes(3,:)';

% The regularization could do with some work to optimize it.
%-----------------------------------------------------------------------
IC0 = eye(6)*0.000001;

for y=1:DIM(2)
	Xmotiony  = (Matrixes(5,:)*y + Matrixes( 9,:)*j + Matrixes(13,:))';
	Ymotiony  = (Matrixes(6,:)*y + Matrixes(10,:)*j + Matrixes(14,:))'; % - y;
	Zmotiony  = (Matrixes(7,:)*y + Matrixes(11,:)*j + Matrixes(15,:))'; % - j;

	for x = 1:DIM(1)
		xy = x+(y-1)*DIM(1);
		if (Mask(xy))

			% Motion relative to the first image in the series
			Xmo = Xmotiony + dXmotion*x; % - x;
			Ymo = Ymotiony + dYmotion*x;
			Zmo = Zmotiony + dZmotion*x;

			% Basis functions appropriate for sinc interpolation
			A = [cos(2*pi*Xmo) sin(2*pi*Xmo) ...
			     cos(2*pi*Ymo) sin(2*pi*Ymo) ...
			     cos(2*pi*Zmo) sin(2*pi*Zmo)  ];

			% centre the signal and covariates
			Tac = X(:,xy) - Integral(xy,j);
			A   = A       - ones(size(A,1),1) * mean(A);

			X(:,xy) = X(:,xy) - A*((A'*A + IC0)\(A'*Tac));
		end
	end
end
			X = X';
			start_vol = 1;
			if (any(Flags == 'n')) start_vol = 2; end

			for i = start_vol:size(P,1)
				d = X(:,i)/Headers(i,4);
				if any(Flags == 'k'), d = d.*Mask; end
				if any(Headers(i,5) == [2 4 8]), d = round(d); end
				if Headers(i,5) == 2,
					d(find(d <   0)) =   0;
					d(find(d > 255)) = 255;
				elseif (Headers(i,5) == 4)
					d(find(d >  32767)) =  32767;
					d(find(d < -32768)) = -32768;
				elseif (Headers(i,5) == 8)
					d(find(d >  2^31-1)) =  2^31-1;
					d(find(d < -2^31  )) = -2^31;
				end
				p = deblank(P(i,:));
				q      = max([find(p == '/') 0]);
				q      = [p(1:q) 'r' p((q + 1):length(p))];

				fp = fopen(q, open_mode);
				if (fp ~= -1)
					if fwrite(fp,d,spm_type(Headers(i,5))) ~= prod(size(d))
						fclose(fp);
						write_error_message(q);
						error(['Error writing ' q '. Check your disk space.']);
					end
					fclose(fp);
				else
					open_error_message(q);
					error(['Error opening ' q '. Check that you have write permission.']);
				end

			end
		end
		open_mode = 'a';
		spm_progress_bar('Set',j);
	end

	% Write headers and matrixes
	%------------------------------------------------------------------
	if (~any(Flags == 'N'))
		for i = start_vol:size(P,1)
			p = deblank(P(i,:));
			q = max([find(p == '/') 0]);
			q = [p(1:q) 'r' p((q + 1):length(p))];
			spm_hwrite(q,DIM,VOX,Headers(i,4),Headers(i,5),0,ORIGIN,'spm - realigned');
			spm_get_space(q,M);
		end
	end

	if any(Flags == 'i')
		% Write integral image (16 bit signed)
		%-----------------------------------------------------------
		mx = max(max(Integral));
		SCALE  = mx/32767;
		p      = deblank(P(1,:));
		q      = max([find(p == '/') 0]);
		q      = [p(1:q) 'mean' p((q + 1):length(p))];
		fp = fopen(q,'w');
		if (fp ~= -1)
			for j = 1:DIM(3)
				d = round(Integral(:,j)/SCALE);
				d(find(d < -32768)) = -32768;
				if fwrite(fp,d,spm_type(4)) ~= prod(size(d))
					fclose(fp);
					write_error_message(q);
					error(['Error writing ' q '. Check your disk space.']);
				end
			end
			spm_hwrite(q,DIM,VOX,SCALE,4,0,ORIGIN,'spm - mean image');
			spm_get_space(q,M);
		else
			open_error_message(q);
			error(['Error opening ' q '. Check that you have write permission.']);
		end
	end

	for i = 1:size(P,1); spm_unmap(V(:,i)); end
end


spm_figure('Clear','Interactive');
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
