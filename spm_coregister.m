function spm_coregister(PGF, PFF, PGG, PFG, others,flags)
% Between and within mode image coregistration.
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

if nargin == 2, do_disp(PGF, PFF); return; end;

show_results = 0;

% Get the space of the images
%-----------------------------------------------------------------------
MG = spm_get_space(deblank(PGF(1,:)));
MF = spm_get_space(deblank(PFF(1,:)));

if strcmp(PGG,PFG), 	% Same modality
	tic;
	MM = coreg_within_mod(PGF, PFF);

	if show_results==1,
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

	tic;
	[MM,MGR,MFR,MTA] = coreg_step1(PGF, PFF, PGG, PFG);
	if show_results==1,
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
	end;

	global sptl_QckCrg
	if ~any(any(sptl_QckCrg) | any(flags=='n'))
		MM = coreg_steps23(PGF, PFF, PGG, PFG, MM,MGR,MFR,MTA);
		if show_results==1,
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
		end;
	end;
end;

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

linfun('Creating Display..');
do_disp(PGF, PFF);
spm_print;

linfun(' ');
spm_figure('Clear','Interactive');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function MM = coreg_within_mod(PGF, PFF)
% Within modality registration

inameG = [];
inameF = [];

% Smooth the images
%-----------------------------------------------------------------------
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
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [MM,MGR,MFR,MTA] = coreg_step1(PGF, PFF, PGG, PFG)
% Rough between modality coregistration
%-----------------------------------------------------------------------
linfun(['Smoothing "' spm_str_manip(PGF(1,:),'a40') '"..']);
spm_smooth(PGF(1,:),fullfile('.','spm_coreg_tmpG.img'),8);
linfun(['Smoothing "' spm_str_manip(PFF(1,:),'a40') '"..']);
spm_smooth(PFF(1,:),fullfile('.','spm_coreg_tmpF.img'),8);

PPF = str2mat(	fullfile('.','spm_coreg_tmpG.img'),...
		fullfile('.','spm_coreg_tmpF.img')    );
PPG = str2mat(PGG(1,:), PFG(1,:));
VF  = spm_vol(PPF);
VG  = spm_vol(PPG);

spm_chi2_plot('Init','Rough Coregistration','Convergence','Iteration');
% Be careful here with the order of the matrix multiplications.
global sptl_Ornt;
if prod(size(sptl_Ornt)) == 12,
	params = [sptl_Ornt(1:6) sptl_Ornt(1:6) sptl_Ornt(7:12) 1 1]';
else,
	params = [zeros(1,12) 1 1 1 0 0 0  1 1]';
end;
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
return;
%_______________________________________________________________________

%_______________________________________________________________________
function MM = coreg_steps23(PGF, PFF, PGG, PFG, MM,MGR,MFR,MTA)
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
spm_chi2_plot('Init','Coregistering segments','Convergence','Iteration');
params = [spm_imatrix(MM) 1 1]';
params = spm_affsub3('rigid2', VG(1:2), VF(1:2), 1, 6, params);
MM = spm_matrix(params(1:12));
spm_chi2_plot('Clear');

% Delete temporary files
%-----------------------------------------------------------------------
for P = str2mat(PPG,PPF)',
	delete_image(P');
end;
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
%_______________________________________________________________________

%_______________________________________________________________________
function do_disp(PGF, PFF)
% Do the graphical display

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
for i=1:3,
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
end;
drawnow;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function linfun(x)
fprintf('%-60s%s', x,sprintf('\b')*ones(1,60));
return;
%_______________________________________________________________________
