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

do_disp = 1;

global SWD
DIR1 = [SWD '/coreg/'];

global SWD sptl_WhchPtn
if (nargin == 0)
	% Act as user interface if there are no arguments
	%_______________________________________________________________________

	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Coregistration');
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
		flags = '';
		templates = str2mat([DIR1 'PET.img'], ...
			[DIR1 'T1.img'], [DIR1 'T2.img'],...
			[DIR1 'EPI.img'],[DIR1 'Transm.img']);

		% Get modality of target
		%-----------------------------------------------------------------------
		respt = spm_input('Modality of first target image?',3,'m',...
			'target - PET|target - T1 MRI|target - T2 MRI|target - EPI|target - Transm',...
			[1 2 3 4 5],1);
		PGG = deblank(templates(respt,:));

		% Get modality of object
		%-----------------------------------------------------------------------
		respo = spm_input('Modality of first object image?',4,'m',...
			'object - PET|object - T1 MRI|object - T2 MRI|object - EPI|object - Transm',...
			[1 2 3 4 5],2);
		PFG = deblank(templates(respo,:));

		if (respo==5 | respt==5),
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

	if (p==2)
		for i = 1:nsubjects
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
		end
	end


	% For each subject, recursively call the program to perform the
	% registration.
	%-----------------------------------------------------------------------
	for i=1:nsubjects
		set(spm_figure('FindWin','Interactive'),...
			'Name',['Coregistering subject ' num2str(i)],'Pointer','Watch');
		drawnow;
		eval(['PGF    =    PGF' num2str(i) ';']);
		eval(['PFF    =    PFF' num2str(i) ';']);
		eval(['others = others' num2str(i) ';']);

		if (p == 1 | p == 3)
			spm_coregister(PGF, PFF, PGG, PFG, others,flags);
		end
		if (p == 2 | p == 3)
			% Write the coregistered images
			%-----------------------------------------------------------------------
			P = str2mat(PGF(1,:),PFF);
			if prod(size(others))>0
				P = str2mat(P,others);
			end
			spm_realign('Reslice',P,'n');
		end
		spm_figure('Clear','Interactive'); drawnow;
	end
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


	VF = spm_map(deblank(PFF(1,:)));
	MF = spm_get_space(deblank(PFF(1,:)));
	VG = spm_map(deblank(PGF(1,:)));
	MG = spm_get_space(deblank(PGF(1,:)));

	Q = MG\MF;
	text(0,0.85, sprintf('X1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(1,:)),'Parent',ax);
	text(0,0.80, sprintf('Y1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(2,:)),'Parent',ax);
	text(0,0.75, sprintf('Z1 = %0.2f*X + %0.2f*Y + %0.2f*Z + %0.2f',Q(3,:)),'Parent',ax);
	text(0.25,0.65, spm_str_manip(PFF,'k22'),...
		'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold','Parent',ax);
	text(0.75,0.65, spm_str_manip(PGF,'k22'),...
		'HorizontalAlignment','center','FontSize',14,'FontWeight','Bold','Parent',ax);

	%-----------------------------------------------------------------------
	for i=1:3
		M = spm_matrix([0 0 i*VF(3)/4]);
		ax=axes('Position',[0.1 0.03+(0.75-0.25*i) 0.4 0.25],'Visible','off','Parent',fig);
		set(fig,'CurrentAxes',ax);
		img1 = spm_slice_vol(VF(:,1),M,VF(1:2,1),1);
		hold on;
		imagesc(img1');
		contour(img1',3,'r');
		axis('off','image');

		ax=axes('Position',[0.5 0.03+(0.75-0.25*i) 0.4 0.25],'Visible','off','Parent',fig);
		set(fig,'CurrentAxes',ax); 
		img2 = spm_slice_vol(VG(:,1),MG\MF*M,VF(1:2,1),1);
		hold on;
		imagesc(img2');
		contour(img1',3,'r');
		axis('off','image');
	end
	spm_unmap_vol(VG);
	spm_unmap_vol(VF);
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
if strcmp(PGG,PFG) 	% Same modality


	inameG = [];
	inameF = [];

	% Smooth the images
	%-----------------------------------------------------------------------
	disp('Smoothing')
	disp(PGF);
	disp(PFF);
	for i=1:size(PFF,1)
		inameG = str2mat(inameG, [spm_str_manip(PGF(i,:),'rd') '_tmpG.img']);
		spm_smooth(deblank(PGF(i,:)),deblank(inameG(i+1,:)),8);

		inameF = str2mat(inameF, [spm_str_manip(PFF(i,:),'rd') '_tmpF.img']);
		spm_smooth(deblank(PFF(i,:)),deblank(inameF(i+1,:)),8);
	end
	inameG = inameG(2:size(inameG,1),:);
	inameF = inameF(2:size(inameF,1),:);

	% Coregister the images together.
	%-----------------------------------------------------------------------
	disp('Coregistering')
	spm_chi2_plot('Init','Coregistering','Convergence','Iteration');
	params = spm_affsub3('rigid2', inameG, inameF, 1, 8);
	params = spm_affsub3('rigid2', inameG, inameF, 1, 6, params);
	MM     = spm_matrix(params);
	spm_chi2_plot('Clear');

	% Delete temporary files
	%-----------------------------------------------------------------------
	for i=1:size(PFF,1)
		iname  = spm_str_manip(inameG(i,:),'sd');
		spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat']);

		iname  = spm_str_manip(inameF(i,:),'sd');
		spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat']);
	end

	if do_disp==1,
		im1 = PGF(1,:);
		im2 = PFF(1,:);
		M1=spm_get_space(im1);
		M2=spm_get_space(im2);
		d1=spm_hread(im1);
		d2=spm_hread(im2);
		fprintf('--------------------------------------------------------------\n');
		fprintf('Method: 3\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1,im2);
		disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
	end

else 	% Different modalities

	% Rough coregistration
	%-----------------------------------------------------------------------
	disp('Smoothing');
	spm_smooth(PGF(1,:),'./spm_coreg_tmpG.img',8);
	spm_smooth(PFF(1,:),'./spm_coreg_tmpF.img',8);

	PPF = str2mat('./spm_coreg_tmpG.img', './spm_coreg_tmpF.img');
	PPG = str2mat(PGG(1,:), PFG(1,:));

	disp('Rough coregistration');
	spm_chi2_plot('Init','Rough Coregistration','Convergence','Iteration');
	% Be careful here with the order of the matrix multiplications.
	global sptl_Ornt;
	if prod(size(sptl_Ornt)) == 12
		tmp = spm_imatrix(inv(spm_matrix(sptl_Ornt(1:12))));
		params = [tmp(1:6) tmp(1:6) tmp(7:12) 1 1]';
	else
		params = [zeros(1,12) 1 1 1 0 0 0  1 1]';
	end
	params = spm_affsub3('register1', PPF, PPG, 1, 8 , params);
	params = spm_affsub3('register1', PPF, PPG, 1, 6 , params);
	spm_chi2_plot('Clear');

	spm_unlink ./spm_coreg_tmpG.img ./spm_coreg_tmpG.hdr ./spm_coreg_tmpG.mat
	spm_unlink ./spm_coreg_tmpF.img ./spm_coreg_tmpF.hdr ./spm_coreg_tmpF.mat

	MM = spm_matrix(params([1:6 13:18]))/spm_matrix(params([7:12 13:18]));

	if do_disp==1,
		im1 = PGF(1,:);
		im2 = PFF(1,:);
		M1=spm_get_space(im1);
		M2=spm_get_space(im2);
		d1=spm_hread(im1);
		d2=spm_hread(im2);
		fprintf('--------------------------------------------------------------\n');
		fprintf('Method: 1\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1,im2);
		disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
	end

	global QUICK_COREG
	if ~any(QUICK_COREG == 1 | any(flags=='n'))

		% Partition the target image(s) into smoothed segments
		%-----------------------------------------------------------------------
		disp('Segmenting and smoothing:')
		disp(PGF);
		spm_segment(PGF,inv(spm_matrix(params([1:6 13:18]))),'ft');

		PPG = [ [spm_str_manip(PGF(1,:),'rd') '_sseg_tmp1.img']
			[spm_str_manip(PGF(1,:),'rd') '_sseg_tmp2.img']];

		% Partition the object image(s) into smoothed segments
		%-----------------------------------------------------------------------
		disp('Segmenting and smoothing:')
		disp(PFF);
		spm_segment(PFF,inv(spm_matrix(params([7:12 13:18]))),'ft');

		PPF = [ [spm_str_manip(PFF(1,:),'rd') '_sseg_tmp1.img']
			[spm_str_manip(PFF(1,:),'rd') '_sseg_tmp2.img']];

		% Coregister the segments together
		%-----------------------------------------------------------------------
		disp('Coregistering segments')
		spm_chi2_plot('Init','Coregistering segments','Convergence','Iteration');
		params = [spm_imatrix(MM) 1 1]';
		params = spm_affsub3('rigid2', PPG(1:2,:), PPF(1:2,:), 1, 6, params);
		MM = spm_matrix(params(1:12));
		spm_chi2_plot('Clear');

		% Delete temporary files
		%-----------------------------------------------------------------------
		for P = str2mat(PPG,PPF)';
			iname2 = spm_str_manip(P','rd');
			spm_unlink([iname2 '.img'], [iname2 '.hdr'], [iname2 '.mat']);
		end

		if do_disp==1,
			im1 = PGF(1,:);
			im2 = PFF(1,:);
			M1=spm_get_space(im1);
			M2=spm_get_space(im2);
			d1=spm_hread(im1);
			d2=spm_hread(im2);
			fprintf('--------------------------------------------------------------\n');
			fprintf('Method: 2\nDate: %s\nPatient Number: ??\nFrom: %s\nTo:   %s\n\n', date,im1,im2);
			disp_coreg_params(M1,MM*M2,d1,d2)
fprintf('time=%g seconds\n',toc);
		end
	end
end

% Save parameters. The matrixes are all loaded into memory, before
% the transformations are applied - just in case any images have
% been included more than once in the list.
%-----------------------------------------------------------------------
Images   = str2mat(PFF,others);
Matrixes = zeros(16,size(Images,1));

for i=1:size(Images,1)
	M = spm_get_space(deblank(Images(i,:)));
	Matrixes(:,i) = M(:);
end

for i=1:size(Images,1)
	M = reshape(Matrixes(:,i),4,4);
	spm_get_space(deblank(Images(i,:)), MM*M);
end

spm_coregister(PGF, PFF);
spm_print;

disp('Done');
spm_figure('Clear','Interactive');
return




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

