function spm_mireg_ui
% Between modality coregistration using Mutual Information
%
% ____________________________________________________________________________
% 
%  The TARGET image is the image to which the OBJECT image is realigned.
%  If there are any OTHER images, then the same transformations are applied to
%  these images as are applied to the OBJECT image.
% 
%  eg 1) to realign a structural MR image to a sequence of PET images:
%   TARGET: meanPET1.img
%   OBJECT: MRI.img
%    OTHER: -
% 
%  eg 2) to realign a sequence of PET images to a structural MR image:
%   TARGET: MRI.img
%   OBJECT: meanPET1.img
%    OTHER: PET1.img PET2.img PET3.img etc...
% ____________________________________________________________________________
%
% The registration method used here is based on the work described in:
% A Collignon, F Maes, D Delaere, D Vandermeulen, P Suetens & G Marchal
% (1995) "Automated Multi-modality Image Registration Based On
% Information Theory". In the proceedings of Information Processing in
% Medical Imaging (1995).  Y. Bizais et al. (eds.).  Kluwer Academic
% Publishers.
% 
% The mutual Information is essentially given by:
% H  = H/(sum(H(:))+eps);
% s1 = sum(H,1);
% s2 = sum(H,2);
% H  = H.*log2((H+eps)./(s2*s1+eps));
% mi = sum(H(:));
% 
% where H is a 256x256 histogram, and mi is the mutual information.
% As an attempt to improve the convergence properties, the algorithm
% minimises exp(-mi).  This is what gets plotted for each line
% minimisation of the optimisation.
% 
% The optimisation has been taken from "Numerical Recipes in C"
% (1992, 2nd Ed.), by WH Press, SA Teukolsky, WT Vetterling &
% BP Flannery.
% 
% At the end, the voxel-to-voxel affine transformation matrix is
% displayed, along with the histograms for the images in the original
% orientations, and the final orientations.  The registered images are
% displayed at the bottom.
%
% Realignment parameters are stored in the ".mat" files of the "object" and
% the "other" images.
% _______________________________________________________________________
% %W% John Ashburner %E%

SPMid = spm('FnBanner',mfilename,'2.8');
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MI Coregister');
spm_help('!ContextHelp','spm_mireg_ui.m');

% get number of subjects
nsubjects = spm_input('number of subjects',1, 'e', 1);
if nsubjects < 1,
	spm_figure('Clear','Interactive');
	return;
end;

global sptl_WhchPtn
if sptl_WhchPtn ~= 1,
	p = spm_input('Which option?',2,'m',...
	'Coregister only|Reslice Only|Coregister & Reslice',...
	[1 2 3],3);
else,
	p = 3;
end;


if p == 1 | p == 3,
	for i = 1:nsubjects,
		mireg(i) = struct('VG',[],'VF',[],'PO','');
		
		% select target(s)
		PG = spm_get(1,'.img', ['select target image for subject ' num2str(i)]);
		mireg(i).VG = spm_vol(PG);
		
		% select object(s)
		PF = spm_get(1,'.img', ['select object image for subject ' num2str(i)]);
		mireg(i).VF = spm_vol(PF);

		% select others
		PO = spm_get(Inf,'.img', ['select other images for subject ' num2str(i)]);
		if isempty(PO),
			mireg(i).PO = PF;
		else,
			mireg(i).PO = str2mat(PF,PO);
		end;
	end;
end;

if p==2,
	for i = 1:nsubjects,
		mireg(i) = struct('VG',[],'VF',[],'VO',[]);
		% select target space
		PG = spm_get(1,'.img', ['select image defining space for subject ' num2str(i)]);
		mireg(i).VG = spm_vol(PG);

		PO = spm_get(Inf,'.img', ['select images to reslice ' num2str(i)]);
		mireg(i).VO = PO;
	end;
end;

% For each subject, recursively call the program to perform the
% registration.
%-----------------------------------------------------------------------
spm('Pointer','Watch')
for i=1:nsubjects,
	spm('FigName',['MI Coregister: working on subj ' num2str(i)],Finter,CmdLine);
	fprintf('\nCoregistering Subject %d\n', i);

	if p == 1 | p == 3,
		x = spm_mireg(mireg(i).VG, mireg(i).VF);
		M = inv(spm_matrix(x));

		MM = zeros(4,4,size(mireg(i).PO,1));
		for j=1:size(mireg(i).PO,1),
			MM(:,:,j) = spm_get_space(deblank(mireg(i).PO(j,:)));
		end;
		for j=1:size(mireg(i).PO,1),
			spm_get_space(deblank(mireg(i).PO(j,:)), M*MM(:,:,j));
		end;

		spm_print;
	end;
	if p == 2 | p == 3,
		% Write the coregistered images
		P = str2mat(mireg(i).VG.fname,mireg(i).PO);
		spm_reslice(P,struct('mask',0,'mean',0,'hold',1,'which',1));
	end;
end;
spm('FigName','MI Coregister: done',Finter,CmdLine);
spm('Pointer');
return;
