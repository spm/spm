% %W% John Ashburner

function spm_sn3d_ui
set(2,'Name','','Pointer','Arrow'); drawnow;
figure(2);clf

global SWD CWD;
a1 = spm_input('Which option?',1,'m','Determine Parameters Only|Write Normalised Only|Determine Parameters & Write Normalised',[1 2 3]);
nsubjects = spm_input('# Subjects',1);

for i=1:nsubjects
	if (a1 == 1 | a1 == 3)
		P = spm_get(1,'.img',['Image to Normalise - subject ' num2str(i)]);
		eval(['P' num2str(i) '=P;']);
		matname = [spm_str_manip(P,'sd') '_sn3d.mat'];
	else
		matname = spm_get(1,'_sn3d.mat',['select Normalisation Parameter Set']);
	end
	eval(['matname' num2str(i) '=matname;']);

	if (a1 == 2 | a1 == 3)
		P = spm_get(Inf,'.img',['select Images to Write ']);
		eval(['PP' num2str(i) '=P;']);
	end
end

pos = 1;
if (a1 == 1 | a1 == 3)

	% Get template(s)
	ok = 0;
	while (~ok)
		Template = spm_get(Inf,'.img',['select Template(s) '],'', SWD);
		if (size(Template,1)>0)
			dims = zeros(size(Template,1),9);
			for i=1:size(Template,1)
				[dim vox dummy dummy dummy origin dummy] = spm_hread(deblank(Template(i,:)));
				dims(i,:) = [dim vox origin];
			end
			if size(dims,1) == 1 | ~any(diff(dims)) ok = 1; end
		end
	end
	a2 = spm_input('Normalisation Type?',2,'m','Affine Only|Affine & Nonlinear',[1 2]);
	pos = 2;
	if (a2 == 2)
		iterations = spm_input('# Nonlinear Iterations?',3,'m','1  iteration|3  iterations|8  iterations|12 iterations|16 iterations',[1 3 8 12 16]);
		nbasis = [];
		while prod(size(nbasis)) ~= 3 | any(nbasis < 1) | prod(nbasis) > 1000
			nbasis = spm_input('# Basis Functions (x y z)',4);
			nbasis = nbasis(:)';
		end
		smoothness = spm_input('Deformation Smoothness',5,'e',0.001);
		pos = 5;
	else
		nbasis     = [0 0 0];
		iterations = 0;
		smoothness = 0;
	end
end

if (a1 == 2 | a1 == 3)

	Hold = spm_input('Interpolation Method?',pos+1,'m','Nearest Neighbour|Bilinear Interpolation|Rough   Sinc Interpolation|Medium  Sinc Interpolation (slow)|Better Sinc Interpolation (very slow)',[0 1 3 5 8]);
		directions = 'XYZ';
		if (spm_input('Use default Bounding Box?',pos+2,'y/n')=='y')
			bb = [ [-64 64]' [-104 68]' [-28 72]' ];
		else
			bb = zeros(2,1);
			for d=1:3
				bbx = [];
				while size(bbx,1) ~= 2
					bbx = spm_input(['Bounding Box ' directions(d) ], pos+1+d);
					bbx = bbx(:);
				end
				bb(:,d) = bbx;
			end
			pos = pos+2;
		end
		Vox = [];
		while size(Vox,2) ~= 3
			Vox = spm_input('Voxel Sizes ',pos+3);
			Vox = Vox(:)';
		end
else
	bb     = [ [-64 64]' [-104 68]' [-28 72]' ];
	Vox    = [2 2 4];
end

set(2,'Name','Normalising','Pointer','Watch'); drawnow;
if (a1 == 1 | a1 == 3)
	for i=1:nsubjects
		eval(['matname=matname' num2str(i) ';']);
		eval(['P=P' num2str(i) ';']);
		spm_sn3d(P,matname,bb,Vox,[nbasis iterations 7 smoothness],Template);
	end
end
set(2,'Name','Writing     Normalised','Pointer','Watch'); drawnow;

if (a1 == 2 | a1 == 3)
	for i=1:nsubjects
		eval(['matname=matname' num2str(i) ';']);
		eval(['P=PP' num2str(i) ';']);
		spm_write_sn(P,matname,bb,Vox,Hold);
	end
end

set(2,'Name','','Pointer','Arrow'); drawnow;
