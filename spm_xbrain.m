function spm_xbrain
% Brain extraction.
% FORMAT spm_xbrain
% A rough and ready routine for extracting the brain from segmented
% images.  It begins by taking the white matter, and eroding it a
% couple of times to get rid of any odd voxels.  The algorithm continues
% on to do conditional dilations for several iterations, where the
% condition is based upon gray or white matter being present.
%
% The extracted brain is written to "brain.img" in the current
% directory.
%
% The function took less than an afternoon to develop, implement and
% test - so don't expect it to work every time.  It is also very slow
% and could be speeded up considerably if re-coded in C.
%
%_______________________________________________________________________
% %W% John Ashburner %E%

linfun = inline('fprintf([''%-60s%s''],x,[sprintf(''\b'')*ones(1,60)])','x');

P=spm_get(2,'*.img','Select gray and white matter images');

VG=spm_vol(P(1,:));
VW=spm_vol(P(2,:));
br=zeros(VW.dim(1:3));
for i=1:VW.dim(3),
	linfun(['Initialising plane ' num2str(i)]);
	br(:,:,i) = spm_slice_vol(VW,spm_matrix([0 0 i]),VW.dim(1:2),1);
end;

kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;
th = 0.6;
spm_progress_bar('Init',25,'Extracting Brain','Iterations completed');
for j=1:25,
	if j>2, th=0.2; end;
	for i=1:VW.dim(3),
		br(:,:,i) = br(:,:,i)>th;
	end;
	for i=1:VW.dim(3),
		linfun(['Iteration ' num2str(j) ' - reading and multiplying plane ' num2str(i)]);
		w=spm_slice_vol(VW,spm_matrix([0 0 i]),VW.dim(1:2),1);
		g=spm_slice_vol(VG,spm_matrix([0 0 i]),VW.dim(1:2),1);
		br(:,:,i) = br(:,:,i).*(w+g);
	end;
	linfun(['Iteration ' num2str(j) ' - convolving']);
	spm_conv_vol(br,br,kx,ky,kz,-[1 1 1]);
	spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');
linfun('Writing volume');
VO=struct('fname','brain.img',...
	  'dim',[VW.dim(1:3) spm_type('uint8')],...
	  'mat',VW.mat,...
	  'pinfo',[1/255 0 0]',...
	  'descrip','extracted brain');
spm_write_vol(VO,br);
linfun(' ');
return;
