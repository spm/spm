function spm_coregister(PGF, PFF, PGG, PFG, others,flags)
% Between and within mode image coregistration.
%
%____________________________________________________________________________
%
% The TARGET image is the image to which the OBJECT image is realigned.
% If there are any OTHER images, then the same transformations are applied to
% these images as are applied to the OBJECT image.
%
% eg 1) to realign a structural MR image to a sequence of PET images:
%  TARGET: meanPET1.img
%  OBJECT: MRI.img
%   OTHER: -
%
% eg 2) to realign a sequence of PET images to a structural MR image:
%  TARGET: MRI.img
%  OBJECT: meanPET1.img
%   OTHER: PET1.img PET2.img PET3.img etc...
%____________________________________________________________________________
%
% The program has two modes of operation:
% 1) If the modalities of the target image(s) and the object image(s) are
%    the same, then the program performs within mode coregistration by
%    minimising the sum of squares difference between the target and object.
%
% 2) If the modalities differ, then the following is performed:
%    i)   Affine normalisation of object to a template of the same modality,
%         and affine normalisation of the target to a template of the same
%         modality.  Only the parameters which describe rigid body
%         transformations are allowed to differ between these normalisations.
%         This produces a rough coregistration of the images.
%    ii)  The images are partitioned into gray matter, white matter, csf and
%         (possibly) scalp using spm_segment.m.  The mappings from images to
%         templates derived from the previous step are used to map from the
%         images to a set of a-priori probability images of GM, WM and CSF.
%    iii) These partitions are then registered together simultaneously, using
%         the results of step i as a starting estimate.
%
% Realignment parameters are stored in the ".mat" files of the "object" and
% the "other" images.
%____________________________________________________________________________
% Refs:
%
% Ashburner J & Friston KJ (1997) Multimodal Image Coregistration and
% Partitioning - a Unified Framework. NeuroImage 6:209-217
%
% Ashburner J, Neelin P, Collins DL, Evans AC & Friston KJ (1997)
% Incorporating Prior Knowledge into Image Registration. NeuroImage 6:344-352
%____________________________________________________________________________
% %W% John Ashburner FIL %E%

% programers notes
%
% FORMAT spm_coregister(PGF, PFF, PGG, PFG, others,flags)
% 		PGF    - target image (the image to act as template).
% 		PFF    - object image (the image to reposition).
% 		PGG    - image to affine normalise target image to.
% 		PFG    - image to affine normalise object image to.
% 		others - other images to apply same transformation to.
%		flags  - any flags
%			n - only do first pass rough coregistration.
% FORMAT spm_coregister(PGF, PFF)
% 		PGF    - target image.
% 		PFF    - object image.
%
% This form simply does a graphic display of how well the coregistration has
% worked.

linfun = inline('fprintf(''  %-60s%s'', x,sprintf(''\b'')*ones(1,60))');

do_disp = 0;

global SWD
DIR1 = [fullfile(SWD,'templates'),filesep];

global SWD sptl_WhchPtn
if (nargin == 0)
	% Act as user interface if there are no arguments
	%_______________________________________________________________________

	SPMid = spm('FnBanner',mfilename,'%I%');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Coregister');
	spm_help('!ContextHelp','spm_coregister.m');

	% get number of subjects
	nsubjects = spm_input('number of subjects',1, 'e', 1);
	if (nsubjects < 1)
		spm_figure('Clear','Interactive');
		return;
	end

	if (sptl_WhchPtn ~= 1)
		p = spm_input('Which option?',2,'m',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3],3);
	else
		p = 3;
	end

	if (p == 1 | p == 3)
		flags = ' ';
		templates = str2mat([DIR1 'PET.img'], ...
			[DIR1 'T1.img'], [DIR1 'T2.img'],...
			[DIR1 'PD.img'],[DIR1 'EPI.img'],[DIR1 'Transm.img']);

		% Get modality of target
		%-----------------------------------------------------------------------
		respt = spm_input('Modality of first target image?',3,'m',...
			'target - PET|target - T1 MRI|target - T2 MRI|target - PD MRI|target - EPI|target - Transm',...
			[1 2 3 4 5 6],1);
		PGG = deblank(templates(respt,:));

		% Get modality of object
		%-----------------------------------------------------------------------
		respo = spm_input('Modality of first object image?',4,'m',...
			'object - PET|object - T1 MRI|object - T2 MRI|object - PD MRI|object - EPI|object - Transm',...
			[1 2 3 4 5 6],2);
		PFG = deblank(templates(respo,:));

		if (respo==6 | respt==6),
			% only perform the first step of the registration
			% because transmission/CT images do not segment very
			% well.
			flags = [flags 'n'];
		end
		if (respt == respo)
			n_images = 1;
		else
			n_images = Inf;
		end

		for i = 1:nsubjects
			% select target(s)
			PGF = [];
			while size(PGF,1)<1
				PGF = spm_get(n_images,'.img',...
					['select target image for subject ' num2str(i)]);
			end

			% select object(s)
			PFF = [];
				while size(PFF,1)<1
				PFF = spm_get(n_images,'.img',...
					['select object image for subject ' num2str(i)]);
			end

			% select others
			others = spm_get(Inf,'.img',...
				['select other images for subject ' num2str(i)]);

			eval(['PGF'    num2str(i) ' = PGF;']);
			eval(['PFF'    num2str(i) ' = PFF;']);
			eval(['others' num2str(i) ' = others;']);
		end
	end

	if p==2,
		for i = 1:nsubjects,
			% select target space
			PGF = spm_get(1,'.img',...
					['select image defining space for subject ' num2str(i)]);

			% select images to reslice
			PFF = [];
			PFF = spm_get(Inf,'.img',...
				['select images to reslice ' num2str(i)]);

			eval(['PGF'    num2str(i) ' = PGF;']);
			eval(['PFF'    num2str(i) ' = PFF;']);
			eval(['others' num2str(i) ' = [];']);
		end;
	end;


	% For each subject, recursively call the program to perform the
	% registration.
	%-----------------------------------------------------------------------
	spm('Pointer','Watch')
	for i=1:nsubjects,
		spm('FigName',['Coregister: working on subj ' num2str(i)],Finter,CmdLine);
		fprintf('\rCoregistering Subject %d: ', i);

		eval(['PGF    =    PGF' num2str(i) ';']);
		eval(['PFF    =    PFF' num2str(i) ';']);
		eval(['others = others' num2str(i) ';']);

		if p == 1 | p == 3,
			spm_coregister(PGF, PFF, PGG, PFG, others,flags);
		end;
		if p == 2 | p == 3,
			% Write the coregistered images
			%-----------------------------------------------------------------------
			P = str2mat(PGF(1,:),PFF);
			if prod(size(others))>0,
				P = str2mat(P,others);
			end;
			spm_reslice(P,struct('mask',0,'mean',0,'hold',1,'which',1));
		end;
	end;
	fprintf('\r%60s%s', ' ',sprintf('\b')*ones(1,60));
	spm('FigName','Coregister: done',Finter,CmdLine);
	spm('Pointer');
	return;
elseif nargin == 2

	% Do the graphics
	%=======================================================================

	fig = spm_figure('FindWin','Graphics');
	set(0,'CurrentFigure',fig);
	spm_figure('Clear','Graphics');

	ax = axes('Position',[0.1 0.51 0.8 0.45],'Visible','off','Parent',fig);
	text(0.5,0.90, 'Coregistration','FontSize',16,...
		'FontWeight','Bold','HorizontalAlignment','center','Parent',ax);


	VF = spm_vol(PFF(1,:));
	MF = VF.mat;
	VG = spm_vol(PGF(1,:));
	MG = VG.mat;

	Q = MG\MF;
	text(0,0.85, sprintf('X1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(1,:)),'Parent',ax);
	text(0,0.80, sprintf('Y1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(2,:)),'Parent',ax);
	text(0,0.75, sprintf('Z1 = %0.3f*X %+0.3f*Y %+0.3f*Z %+0.3f',Q(3,:)),'Parent',ax);
	text(0.25,0.65, spm_str_manip(PFF,'k22'),...
		'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold','Parent',ax);
	text(0.75,0.65, spm_str_manip(PGF,'k22'),...
		'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold','Parent',ax);

	%-----------------------------------------------------------------------
	for i=1:3
		M = spm_matrix([0 0 i*VF(1).dim(3)/4]);
		ax=axes('Position',[0.1 0.03+(0.75-0.25*i) 0.4 0.25],'Visible','off','Parent',fig);
		set(fig,'CurrentAxes',ax);
		img1 = spm_slice_vol(VF(1),M,VF(1).dim(1:2),1);
		hold on;
		imagesc(img1');
		contour(img1',3,'r');
		axis('off','image');

		ax=axes('Position',[0.5 0.03+(0.75-0.25*i) 0.4 0.25],'Visible','off','Parent',fig);
		set(fig,'CurrentAxes',ax); 
		img2 = spm_slice_vol(VG(1),MG\MF*M,VF.dim(1:2),1);
		hold on;
		imagesc(img2');
		contour(img1',3,'r');
		axis('off','image');
	end
	drawnow;

	return
end


% Do the work
%_______________________________________________________________________

% Get the space of the images
%-----------------------------------------------------------------------
MG = spm_get_space(deblank(PGF(1,:)));
MF = spm_get_space(deblank(PFF(1,:)));

tic;
if strcmp(PGG,PFG), 	% Same modality

	inameG = [];
	inameF = [];

	% Smooth the images
	%-----------------------------------------------------------------------
	disp(PGF);
	disp(PFF);
	for i=1:size(PFF,1),
		inameG = str2mat(inameG, [spm_str_manip(PGF(i,:),'rd') '_tmpG.img']);
		linfun(['Smoothing "' spm_str_manip(PGF(i,:),'a40') '"..']);
		spm_smooth(deblank(PGF(i,:)),deblank(inameG(i+1,:)),8);

		inameF = str2mat(inameF, [spm_str_manip(PFF(i,:),'rd') '_tmpF.img']);
		linfun(['Smoothing "' spm_str_manip(PFF(i,:),'a40') '"..']);
		spm_smooth(deblank(PFF(i,:)),deblank(inameF(i+1,:)),8);
	end;
	inameG = inameG(2:size(inameG,1),:);
	inameF = inameF(2:size(inameF,1),:);
	VG = spm_vol(inameG);
	VF = spm_vol(inameF);
	for i=1:prod(size(VG)),
		VG(i).pinfo(1:2,:) = VG(i).pinfo(1:2,:)/spm_global(VG(i));
	end;
	for i=1:prod(size(VF)),
		VF(i).pinfo(1:2,:) = VF(i).pinfo(1:2,:)/spm_global(VF(i));
	end;

	% Coregister the images together.
	%-----------------------------------------------------------------------
	spm_chi2_plot('Init','Coregistering','Convergence','Iteration');
	linfun('Coarse Coregistration..');
	params = spm_affsub3('rigid2', VG, VF, 1, 8);
	linfun('Fine Coregistration..');
	params = spm_affsub3('rigid2', VG, VF, 1, 6, params);
	MM     = spm_matrix(params);
	spm_chi2_plot('Clear');

	% Delete temporary files
	%-----------------------------------------------------------------------
	for i=1:size(PFF,1)
		delete_image(inameG(i,:));
		delete_image(inameF(i,:));
	end

	if do_disp==1,
		linfun(' ');
		im1 = spm_vol(PGF(1,:));
		im2 = spm_vol(PFF(1,:));
		M1=im1.mat;
		M2=im2.mat;
		d1=im1.dim(1:3);
		d2=im2.dim(1:3);
		fprintf('--------------------------------------------------------------\n');
		fprintf('Method: 3\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1.fname,im2.fname);
		disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
	end;

else 	% Different modalities

	% Rough coregistration
	%-----------------------------------------------------------------------
	linfun(['Smoothing "' spm_str_manip(PGF(1,:),'a40') '"..']);
	spm_smooth(PGF(1,:),fullfile('.','spm_coreg_tmpG.img'),8);
	linfun(['Smoothing "' spm_str_manip(PFF(1,:),'a40') '"..']);
	spm_smooth(PFF(1,:),fullfile('.','spm_coreg_tmpF.img'),8);

	PPF = str2mat(	fullfile('.','spm_coreg_tmpG.img'),...
			fullfile('.','spm_coreg_tmpF.img')    );
	PPG = str2mat(PGG(1,:), PFG(1,:));
	VF = spm_vol(PPF);
	VG = spm_vol(PPG);

	spm_chi2_plot('Init','Rough Coregistration','Convergence','Iteration');
	% Be careful here with the order of the matrix multiplications.
	global sptl_Ornt;
	if prod(size(sptl_Ornt)) == 12
		params = [sptl_Ornt(1:6) sptl_Ornt(1:6) sptl_Ornt(7:12) 1 1]';
	else
		params = [zeros(1,12) 1 1 1 0 0 0  1 1]';
	end
	linfun('First Pass Coregistration - coarse..');
	params = spm_affsub3('register1', VG, VF, 1, 8 , params);
	linfun('First Pass Coregistration - fine..');
	params = spm_affsub3('register1', VG, VF, 1, 6 , params);
	spm_chi2_plot('Clear');

	delete_image('spm_coreg_tmpG');
	delete_image('spm_coreg_tmpF');

	MGR = spm_matrix(params([1:6 ])');
	MFR = spm_matrix(params([7:12])');
	MTA = spm_matrix([0 0 0 0 0 0 params([13:18])']);
	MM  = MGR\MFR; % equivalent to (MTA*MGR)\(MTA*MFR)

	if do_disp==1,
		linfun(' ');
		im1 = spm_vol(PGF(1,:));
		im2 = spm_vol(PFF(1,:));
		M1  = im1.mat;
		M2  = im2.mat;
		d1  = im1.dim(1:3);
		d2  = im2.dim(1:3);
		fprintf('--------------------------------------------------------------\n');
		fprintf('Method: 1\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1.fname,im2.fname);
		disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
	end

	global QUICK_COREG
	if ~any(QUICK_COREG == 1 | any(flags=='n'))

		% Partition the target image(s) into smoothed segments
		%-----------------------------------------------------------------------
		linfun('Segmenting and Smoothing Reference Image(s)..');
		spm_segment(PGF,MTA*MGR,'ft');

		PPG = [ [spm_str_manip(PGF(1,:),'rd') '_sseg_tmp1.img']
			[spm_str_manip(PGF(1,:),'rd') '_sseg_tmp2.img']];
		VG = spm_vol(PPG);

		% Partition the object image(s) into smoothed segments
		%-----------------------------------------------------------------------
		linfun('Segmenting and Smoothing Object Image(s)..');
		spm_segment(PFF,MTA*MFR,'ft');

		PPF = [ [spm_str_manip(PFF(1,:),'rd') '_sseg_tmp1.img']
			[spm_str_manip(PFF(1,:),'rd') '_sseg_tmp2.img']];
		VF = spm_vol(PPF);

		% Coregister the segments together
		%-----------------------------------------------------------------------
		linfun('Coregistering the image segments..');
		disp('Coregistering segments')
		spm_chi2_plot('Init','Coregistering segments','Convergence','Iteration');
		params = [spm_imatrix(MM) 1 1]';
		params = spm_affsub3('rigid2', VG(1:2), VF(1:2), 1, 6, params);
		MM = spm_matrix(params(1:12));
		spm_chi2_plot('Clear');

		% Delete temporary files
		%-----------------------------------------------------------------------
		for P = str2mat(PPG,PPF)';
			delete_image(P');
		end

		if do_disp==1,
			linfun(' ');
			im1 = spm_vol(PGF(1,:));
			im2 = spm_vol(PFF(1,:));
			M1  = im1.mat;
			M2  = im2.mat;
			d1  = im1.dim(1:3);
			d2  = im2.dim(1:3);
			fprintf('--------------------------------------------------------------\n');
			fprintf('Method: 2\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1.fname,im2.fname);
			disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
		end
	end
end

% Save parameters. The matrixes are all loaded into memory, before
% the transformations are applied - just in case any images have
% been included more than once in the list.
%-----------------------------------------------------------------------
linfun('Saving Parameters..');
if isempty(others),
	Images = PFF;
else,
	Images   = str2mat(PFF,others);
end;

Matrixes = zeros(4,4,size(Images,1));

for i=1:size(Images,1)
	M = spm_get_space(deblank(Images(i,:)));
	Matrixes(:,:,i) = M;
end

for i=1:size(Images,1)
	M = Matrixes(:,:,i);
	spm_get_space(deblank(Images(i,:)), MM*M);
end

spm_coregister(PGF, PFF);
spm_print;

linfun(' ');
spm_figure('Clear','Interactive');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function disp_coreg_params(M1,M2,d1,d2)
v1=sqrt(sum(M1(1:3,1:3).^2));
v2=sqrt(sum(M2(1:3,1:3).^2));
M=M1\M2;
c1=[0 0 0;1 0 0;0 1 0;1 1 0;0 0 1;1 0 1;0 1 1;1 1 1]*diag((d1-1).*v1);
ci=[(c1/diag(v1)+1)';ones(1,size(c1,1))];
ci=M2\M1*ci;
ci=ci(1:3,:)';
c2=(ci-1)*diag(v2);

fprintf('Point     x        y        z    new_x    new_y    new_z\n\n');
for i=1:size(c1),
	fprintf('%.1d %9g%9g%9g%9g%9g%9g\n', i, c1(i,:), c2(i,:));
end
fprintf('--------------------------------------------------------------\n');
return
%_______________________________________________________________________
%_______________________________________________________________________
function delete_image(iname)
iname = spm_str_manip(iname,'sd');
spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat'], [iname '.mnc']);
return;
