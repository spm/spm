function spm_coregister(PGF, PFF, PGG, PFG, others)
% Between and within mode image coregistration.
%
% FORMAT spm_coregister
%
% OR     spm_coregister(PGF, PFF, PGG, PFG, others)
% 		PGF    - target image (the image to act as template).
% 		PFF    - object image (the image to reposition).
% 		PGG    - image to affine normalise target image to.
% 		PFG    - image to affine normalise object image to.
% 		others - other images to apply same transformation to.
%____________________________________________________________________________
%
% Without any arguments, spm_coregister acts as its own user interface.
%
% The program has two modes of operation:
% 1) If the modalities of the target image(s) and the object image(s) are
%    the same, then the program performs within mode coregistration by
%    minimising the sum of squares difference between the PGF and PFF.
%
% 2) If the modalities differ, then the following is performed:
%    i)   Affine normalisation of PGF to PGG, and affine normalisation of
%         PFF to PFG. Only the parameters which describe rigid body
%         transformations are allowed to differ in these normalisations.
%         This produces a rough coregistration.
%    ii)  The images are partitioned into gray matter, white matter, csf and
%         (possibly) scalp using spm_segment.m
%    iii) These partitions are then registered together simultaneously.
%
% Realignment parameters are stored in the ".mat" files of the "object" and
% the "other" images.
%____________________________________________________________________________
% %W% John Ashburner FIL %E%


global SWD sptl_WhchPtn
if (nargin == 0)
	% Act as user interface if there are no arguments
	%_______________________________________________________________________

	spm_figure('Clear','Interactive');
	set(spm_figure('FindWin','Interactive'),'Name','Coregistration');

	% get number of subjects
	nsubjects = spm_input('number of subjects',1);

	if (sptl_WhchPtn ~= 1)
		p = spm_input('Which option?',2,'m',...
			'Coregister only|Reslice Only|Coregister & Reslice',...
			[1 2 3]);
	else
		p = 3;
	end

	if (p == 1 | p == 3)
		templates = str2mat([SWD '/templates/PET.img'], [SWD '/templates/T1.img'], [SWD '/templates/T2.img']);

		% Get modality of target
		%-----------------------------------------------------------------------
		respt = spm_input('Modality of first target image?',3,'m',...
			'target - PET|target - T1 MRI|target - T2 MRI',...
			[1 2 3]);
		PGG = deblank(templates(respt,:));

		% Get modality of object
		%-----------------------------------------------------------------------
		respo = spm_input('Modality of first object image?',4,'m',...
			'object - PET|object - T1 MRI|object - T2 MRI',...
			[1 2 3]);
		PFG = deblank(templates(respo,:));

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
			spm_coregister(PGF, PFF, PGG, PFG, others);
		end
		if (p == 2 | p == 3)
			% Write the coregistered images
			%-----------------------------------------------------------------------
			P = str2mat(PGF(1,:),PFF);
			if prod(size(others))>0
				P = str2mat(P,others);
			end
			spm_realign(P,'rn');
		end
		spm_figure('Clear','Interactive'); drawnow;
	end
	return;
end


% Do the work
%_______________________________________________________________________

% Get the space of the images
%-----------------------------------------------------------------------
MG = spm_get_space(deblank(PGF(1,:)));
MF = spm_get_space(deblank(PFF(1,:)));

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
	params = spm_affsub3('rigid2', inameG, inameF, 1, 8);
	params = spm_affsub3('rigid2', inameG, inameF, 1, 4, params);
	MM     = spm_matrix(params);

	% Delete temporary files
	%-----------------------------------------------------------------------
	for i=1:size(PFF,1)
		iname  = spm_str_manip(inameG(i,:),'sd');
		spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat']);

		iname  = spm_str_manip(inameG(i,:),'sd');
		spm_unlink([iname '.img'], [iname '.hdr'], [iname '.mat']);
	end


else 	% Different modalities

	% Rough coregistration
	%-----------------------------------------------------------------------
	disp('Smoothing');
	spm_smooth(PGF(1,:),'./spm_seg_tmpG.img',8);
	spm_smooth(PFF(1,:),'./spm_seg_tmpF.img',8);

	PPF = str2mat('./spm_seg_tmpG.img', './spm_seg_tmpF.img');
	PPG = str2mat(PGG(1,:), PFG(1,:));

	disp('Rough coregistration');
	global sptl_Ornt;
	params = [sptl_Ornt(1:6) sptl_Ornt(1:6) sptl_Ornt(7:12) 1 1]';
	params = spm_affsub3('register1', PPG, PPF, 1, 12, params);
	params = spm_affsub3('register1', PPG, PPF, 1, 8 , params);
	params = spm_affsub3('register1', PPG, PPF, 1, 4 , params);

	spm_unlink ./spm_seg_tmpG.img ./spm_seg_tmpG.hdr ./spm_seg_tmpG.mat
	spm_unlink ./spm_seg_tmpF.img ./spm_seg_tmpF.hdr ./spm_seg_tmpF.mat

	MM = spm_matrix(params([1:6 13:18]))\spm_matrix(params([7:12 13:18]));

	global QUICK_COREG
	if ~any(QUICK_COREG == 1)

		% Partition the target image(s) into smoothed segments
		%-----------------------------------------------------------------------
		disp('Segmenting and smoothing:')
		disp(PGF);
		spm_segment(PGF,spm_matrix(params([1:6 13:18])),'t');
		PPG = [];
		for i=1:4
			iname1 = [spm_str_manip(PGF(1,:),'rd') '_seg_tmp'  num2str(i)];
			iname2 = [spm_str_manip(PGF(1,:),'rd') '_sseg' num2str(i)];
			spm_smooth([iname1 '.img'],[iname2 '.img'],8);
			spm_unlink([iname1 '.img'], [iname1 '.hdr'], [iname1 '.mat']);
			PPG = str2mat(PPG, [iname2 '.img']);
		end
		PPG = PPG(2:5,:);

		% Partition the object image(s) into smoothed segments
		%-----------------------------------------------------------------------
		disp('Segmenting and smoothing:')
		disp(PFF);
		spm_segment(PFF,spm_matrix(params([7:12 13:18])),'t');

		PPF = [];
		for i=1:4
			iname1 = [spm_str_manip(PFF(1,:),'rd') '_seg_tmp'  num2str(i)];
			iname2 = [spm_str_manip(PFF(1,:),'rd') '_sseg' num2str(i)];
			spm_smooth([iname1 '.img'],[iname2 '.img'],8);
			spm_get_space(deblank([iname2 '.img']), MM*spm_get_space([iname2 '.img']));
			spm_unlink([iname1 '.img'], [iname1 '.hdr'], [iname1 '.mat']);
			PPF = str2mat(PPF, [iname2 '.img']);
		end
		PPF = PPF(2:5,:);

		% Coregister the segments together
		%-----------------------------------------------------------------------
		disp('Coregistering segments')
		params = spm_affsub3('rigid1', PPG(1:3,:), PPF(1:3,:), 1, 8);
		params = spm_affsub3('rigid1', PPG(1:3,:), PPF(1:3,:), 1, 4 , params);

		MM = spm_matrix(params(1:12))*MM;

		% Delete temporary files
		%-----------------------------------------------------------------------
		for i=1:4
			iname2 = [spm_str_manip(PGF(1,:),'rd') '_sseg' num2str(i)];
			spm_unlink([iname2 '.img'], [iname2 '.hdr'], [iname2 '.mat']);

			iname2 = [spm_str_manip(PFF(1,:),'rd') '_sseg' num2str(i)];
			spm_unlink([iname2 '.img'], [iname2 '.hdr'], [iname2 '.mat']);
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

disp('Done');
spm_figure('Clear','Interactive');

