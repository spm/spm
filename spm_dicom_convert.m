function spm_dicom_convert(hdr,opts)
% Convert DICOM images into something that SPM can use
% FORMAT spm_dicom_convert(hdr,opts)
% hdr  - a cell array of DICOM headers from spm_dicom_headers
% opts - options
%        'all'      - all DICOM files (default)
%        'mosaic'   - the mosaic images
%        'standard' - standard DICOM files
%        'raw'      - convert raw FIDs (not implemented)
%
% Converted files are written to the current directory
%_______________________________________________________________________
% @(#)spm_dicom_convert.m	2.22 John Ashburner & Jesper Andersson 05/02/09

if nargin<2, opts = 'all'; end;

[images,guff]     = select_tomographic_images(hdr);
[mosaic,standard] = select_mosaic_images(images);

if (strcmp(opts,'all') || strcmp(opts,'mosaic')) && ~isempty(mosaic),
	convert_mosaic(mosaic);
end;
if (strcmp(opts,'all') || strcmp(opts,'standard')) && ~isempty(standard),
	convert_standard(standard);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function convert_mosaic(hdr)
spm_progress_bar('Init',length(hdr),'Writing Mosaic', 'Files written');

for i=1:length(hdr),

	% Output filename
	%-------------------------------------------------------------------
	if checkfields(hdr{i},'SeriesNumber','AcquisitionNumber') 
		fname = sprintf('f%s-%.4d-%.5d-%.6d.img', strip_unwanted(hdr{i}.PatientID),...
			hdr{i}.SeriesNumber, hdr{i}.AcquisitionNumber, hdr{i}.InstanceNumber);
	else
		fname = sprintf('f%s-%.6d.img',strip_unwanted(hdr{i}.PatientID),hdr{i}.InstanceNumber);
	end;

	fname = fullfile(pwd,fname);

	% Image dimensions and data
	%-------------------------------------------------------------------
	nc = hdr{i}.Columns;
	nr = hdr{i}.Rows;

	dim      = [0 0 0];
	dim(3)   = read_NumberOfImagesInMosaic(hdr{i});
	np       = [nc nr]/ceil(sqrt(dim(3)));
	dim(1:2) = np;
	if ~all(np==floor(np)),
		warning('%s: dimension problem [Num Images=%d, Num Cols=%d, Num Rows=%d].',...
			hdr{i}.Filename,dim(3), nc,nr);
		continue;
	end;

	% Apparently, this is not the right way of doing it.
	%np = read_AcquisitionMatrixText(hdr{i});
	%if rem(nc, np(1)) || rem(nr, np(2)),
	%	warning('%s: %dx%d wont fit into %dx%d.',hdr{i}.Filename,...
	%		np(1), np(2), nc,nr);
	%	return;
	%end;
	%dim    = [np read_NumberOfImagesInMosaic(hdr{i})];

	mosaic = read_image_data(hdr{i});
	volume = zeros(dim);
	for j=1:dim(3),
		img = mosaic((1:np(1))+np(1)*rem(j-1,nc/np(1)), (np(2):-1:1)+np(2)*floor((j-1)/(nc/np(1))));
		if ~any(img(:)),
			volume = volume(:,:,1:(j-1));
			break;
		end;
		volume(:,:,j) = img;
	end;
	dim = [size(volume)];
	dt  = [spm_type('int16') spm_platform('bigend')];

	% Orientation information
	%-------------------------------------------------------------------
	% Axial Analyze voxel co-ordinate system:
	% x increases     right to left
	% y increases posterior to anterior
	% z increases  inferior to superior

	% DICOM patient co-ordinate system:
	% x increases     right to left
	% y increases  anterior to posterior
	% z increases  inferior to superior

	% T&T co-ordinate system:
	% x increases      left to right
	% y increases posterior to anterior
	% z increases  inferior to superior

	analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];

	vox    = [hdr{i}.PixelSpacing hdr{i}.SpacingBetweenSlices];
	pos    = hdr{i}.ImagePositionPatient';
	orient = reshape(hdr{i}.ImageOrientationPatient,[3 2]);
	orient(:,3) = null(orient');
	if det(orient)<0, orient(:,3) = -orient(:,3); end;

	% The image position vector is not correct. In dicom this vector points to
	% the upper left corner of the image. Perhaps it is unlucky that this is
	% calculated in the syngo software from the vector pointing to the center of
	% the slice (keep in mind: upper left slice) with the enlarged FoV.
	dicom_to_patient = [orient*diag(vox) pos ; 0 0 0 1];
	truepos          = dicom_to_patient *[(size(mosaic)-dim(1:2))/2 0 1]';
	dicom_to_patient = [orient*diag(vox) truepos(1:3) ; 0 0 0 1];
	patient_to_tal   = diag([-1 -1 1 1]);
	mat              = patient_to_tal*dicom_to_patient*analyze_to_dicom;



	% Maybe flip the image depending on SliceNormalVector from 0029,1010
	%-------------------------------------------------------------------
	SliceNormalVector = read_SliceNormalVector(hdr{i});
	if det([reshape(hdr{i}.ImageOrientationPatient,[3 2]) SliceNormalVector(:)])<0;
		volume = volume(:,:,end:-1:1);
		mat    = mat*[eye(3) [0 0 -(dim(3)-1)]'; 0 0 0 1];
	end;


	% Possibly useful information
	%-------------------------------------------------------------------
	tim = datevec(hdr{i}.AcquisitionTime/(24*60*60));
	descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g Mosaic',...
		hdr{i}.MagneticFieldStrength, hdr{i}.MRAcquisitionType,...
		deblank(hdr{i}.ScanningSequence),...
		hdr{i}.RepetitionTime,hdr{i}.EchoTime,hdr{i}.FlipAngle,...
		datestr(hdr{i}.AcquisitionDate),tim(4),tim(5),tim(6));

	% descrip = [deblank(descrip) '   ' hdr{i}.PatientsName];

	if ~spm_flip_analyze_images,
		mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
		volume = flipdim(volume,1);
	end;

	if isfield(hdr{i},'RescaleSlope') && hdr{i}.RescaleSlope ~= 1,
		volume = volume*hdr{i}.RescaleSlope;
	end;
	if isfield(hdr{i},'RescaleIntercept') && hdr{i}.RescaleIntercept ~= 0,
		volume = volume + hdr{i}.RescaleIntercept;
	end;

	V = struct('fname',fname, 'dim',dim, 'dt',dt, 'mat',mat, 'descrip',descrip);
	spm_write_vol(V,volume);
	spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function convert_standard(hdr)
hdr = sort_into_volumes(hdr);
for i=1:length(hdr),
	write_volume(hdr{i});
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vol = sort_into_volumes(hdr)

%
% First of all, sort into volumes based on relevant
% fields in the header.
%

vol{1}{1} = hdr{1};
for i=2:length(hdr),
	orient = reshape(hdr{i}.ImageOrientationPatient,[3 2]);
	xy1    = hdr{i}.ImagePositionPatient*orient;
	match  = 0;
	for j=1:length(vol),
		orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
		xy2    = vol{j}{1}.ImagePositionPatient*orient;
		dist2  = sum((xy1-xy2).^2);
		if ((~isfield(hdr{i},'SeriesNumber') && ~isfield(vol{j}{1},'SeriesNumber')) ||...
		         (hdr{i}.SeriesNumber            == vol{j}{1}.SeriesNumber)) &&...
		   ((~isfield(hdr{i},'AcquisitionNumber') && ~isfield(vol{j}{1},'AcquisitionNumber')) ||...
		         (hdr{i}.AcquisitionNumber       == vol{j}{1}.AcquisitionNumber)) &&...
	       ((~isfield(hdr{i},'SequenceName') && ~isfield(vol{j}{1},'SequenceName')) ||...
		   strcmp(hdr{i}.SequenceName,              vol{j}{1}.SequenceName)) &&...
	       ((~isfield(hdr{i},'SeriesInstanceUID') && ~isfield(vol{j}{1},'SeriesInstanceUID')) ||...
		   strcmp(hdr{i}.SeriesInstanceUID,         vol{j}{1}.SeriesInstanceUID)) &&...
	       ((~isfield(hdr{i},'TemporalPositionIdentifier') && ~isfield(vol{j}{1},'TemporalPositionIdentifier')) ||... % For Philips
			  hdr{i}.TemporalPositionIdentifier == vol{j}{1}.TemporalPositionIdentifier) &&...
		          hdr{i}.Rows                    == vol{j}{1}.Rows &&...
		          hdr{i}.Columns                 == vol{j}{1}.Columns &&...
		  all(abs(hdr{i}.ImageOrientationPatient -  vol{j}{1}.ImageOrientationPatient)<=1e-5) &&...
		      all(hdr{i}.PixelSpacing            == vol{j}{1}.PixelSpacing & dist2<1e-4) &&...
	        (~isfield(hdr{i},'EchoNumbers')  || ~isfield(vol{j}{1},'EchoNumbers') || ...
		          hdr{i}.EchoNumbers             == vol{j}{1}.EchoNumbers),
			vol{j}{end+1} = hdr{i};
			match = 1;
			break;
		end;
	end;
	if ~match,
		vol{end+1}{1} = hdr{i};
	end;
end;

%
% Secondly, sort volumes into ascending/descending
% slices depending on .ImageOrientationPatient field.
%

vol2 = {};
for j=1:length(vol),
	orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
	proj   = null(orient');
	if det([orient proj])<0, proj = -proj; end;

	z      = zeros(length(vol{j}),1);
	for i=1:length(vol{j}),
		z(i)  = vol{j}{i}.ImagePositionPatient*proj;
	end;
	[z,index] = sort(z);
	vol{j}    = vol{j}(index);
	if length(vol{j})>1,
		dist      = diff(z);
		if any(diff(z)==0)
			% So, if this is a GE scanner we might have
			% the "DTI-series" problem. Lets assume so.

			% if isfield(vol{j}{1},'Manufacturer') &&... 
			%	strcmpi(vol{j}{1}.Manufacturer,'ge medical systems') &&...
			%	isfield(vol{j}{1},'Private_0021_104f') &&...
			%	isfield(vol{j}{1},'Private_0019_1019') &&...
			%	isfield(vol{j}{1},'Private_0019_101b')
			%	vol = sort_ge_into_volumes(vol);
			%	break;
			% else
				tmp = sort_into_vols_again(vol{j});
				vol{j} = tmp{1};
				vol2 = {vol2{:} tmp{2:end}};
			% end
		end;
	end;
end;
vol = {vol{:} vol2{:}};
for j=1:length(vol),
	if length(vol{j})>1,
		orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
		proj   = null(orient');
		if det([orient proj])<0, proj = -proj; end;
		z      = zeros(length(vol{j}),1);
		for i=1:length(vol{j}),
			z(i)  = vol{j}{i}.ImagePositionPatient*proj;
		end;
		[z,index] = sort(z);
		dist      = diff(z);
		if sum((dist-mean(dist)).^2)/length(dist)>0.00001,
			fprintf('***************************************************\n');
			fprintf('* VARIABLE SLICE SPACING                          *\n');
			fprintf('* This may be due to missing DICOM files.         *\n');
			if checkfields(vol{j}{1},'PatientID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
				fprintf('*    %s / %d / %d / %d \n',...
					deblank(vol{j}{1}.PatientID), vol{j}{1}.SeriesNumber, ...
					vol{j}{1}.AcquisitionNumber, vol{j}{1}.InstanceNumber);
				fprintf('*                                                 *\n');
			end;
			fprintf('*  %20.3g                           *\n', dist);
			fprintf('***************************************************\n');
		end;
	end;
end;
%dcm = vol;
%save('dicom_headers.mat','dcm');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vol2 = sort_into_vols_again(volj)
if ~isfield(volj{1},'InstanceNumber'),
    fprintf('***************************************************\n');
    fprintf('* The slices may be all mixed up and the data     *\n');
    fprintf('* not really usable.  Talk to your physicists     *\n');
    fprintf('* about this.                                     *\n');
    fprintf('***************************************************\n');
    vol2 = {volj};
    return;
end;

fprintf('***************************************************\n');
fprintf('* The AcquisitionNumber counter does not appear   *\n');
fprintf('* to be changing from one volume to another.      *\n');
fprintf('* Another possible explanation is that the same   *\n');
fprintf('* DICOM slices are used multiple times.           *\n');
%fprintf('* Talk to your MR sequence developers or scanner  *\n');
%fprintf('* supplier to have this fixed.                    *\n');
fprintf('* The conversion is having to guess how slices    *\n');
fprintf('* should be arranged into volumes.                *\n');
if checkfields(volj{1},'PatientID','SeriesNumber','AcquisitionNumber'),
    fprintf('*    %s / %d / %d\n',...
        deblank(volj{1}.PatientID), volj{1}.SeriesNumber, ...
        volj{1}.AcquisitionNumber);
end;
fprintf('***************************************************\n');

z      = zeros(length(volj),1);
t      = zeros(length(volj),1);
d      = zeros(length(volj),1);
orient = reshape(volj{1}.ImageOrientationPatient,[3 2]);
proj   = null(orient');
if det([orient proj])<0, proj = -proj; end;

for i=1:length(volj),
    z(i)  = volj{i}.ImagePositionPatient*proj;
    t(i)  = volj{i}.InstanceNumber;
end;
msg = 0;
[t,index] = sort(t);
volj      = volj(index);
z         = z(index);
msk       = find(diff(t)==0);
if any(msk),
    % fprintf('***************************************************\n');
    % fprintf('* These files have the same InstanceNumber:       *\n');
    % for i=1:length(msk),
    %    [tmp,nam1,ext1] = fileparts(volj{msk(i)}.Filename);
    %    [tmp,nam2,ext2] = fileparts(volj{msk(i)+1}.Filename);
    %    fprintf('* %s%s = %s%s (%d)\n', nam1,ext1,nam2,ext2, volj{msk(i)}.InstanceNumber);
    % end;
    % fprintf('***************************************************\n');
    index = [logical(1) ; diff(t)~=0];
    t     = t(index);
    z     = z(index);
    d     = d(index);
    volj  = volj(index);
end;

if any(diff(sort(t))~=1), msg = 1; end;
[z,index] = sort(z);
volj      = volj(index);
t         = t(index);
vol2      = {};
while ~all(d),
    i  = find(~d);
    i  = i(1);
    i  = find(z==z(i));
    [t(i),si] = sort(t(i));
    volj(i)   = volj(i(si));
    for i1=1:length(i),
        if length(vol2)<i1, vol2{i1} = {}; end;
        vol2{i1} = {vol2{i1}{:} volj{i(i1)}};
    end;
    d(i) = 1;
end;

msg = 0;
if any(diff(sort(t))~=1), msg = 1; end;
if ~msg,
    len = length(vol2{1});
    for i=2:length(vol2),
        if length(vol2{i}) ~= len,
            msg = 1;
            break;
        end;
    end;
end;
if msg,
    fprintf('***************************************************\n');
    fprintf('* There are missing DICOM files, so the the       *\n');
    fprintf('* resulting volumes may be messed up.             *\n');
    if checkfields(volj{1},'PatientID','SeriesNumber','AcquisitionNumber'),
        fprintf('*    %s / %d / %d\n',...
            deblank(volj{1}.PatientID), volj{1}.SeriesNumber, ...
            volj{1}.AcquisitionNumber);
    end;
    fprintf('***************************************************\n');
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function vol = sort_ge_into_volumes(pvol)
%
% This is a routine that will take the p(rovisional)
% volumes in pvol and subdivide them into more volumes
% based on circumstancial (at best) evidence in GE
% file headers.
%

%
% First of all check how many slices there "should" be
% in each volume. And that there is agreement.
%

nsl = zeros(1,length(pvol));
for i=1:length(pvol)
	nsl(i) = pvol{i}{1}.Private_0021_104f;
	for j=1:length(pvol{i})
		if pvol{i}{j}.Private_0021_104f ~= nsl(i)
			warning('Looks like there is something wrong with the conversion software.');
			vol = pvol;
			return
		end
	end
end

%
% Next use independent method to assess
% number of slices and check agreement.
%

for i=1:length(pvol)
	insl = abs(pvol{i}{1}.Private_0019_101b - pvol{i}{1}.Private_0019_1019);
	insl = round((insl + pvol{i}{1}.SliceThickness) / pvol{i}{1}.SpacingBetweenSlices);
	if insl ~= nsl(i)
		warning('Looks like there is something wrong with the conversion software.');
		vol = pvol;
		return
	end
end

%
% OK, so we have volumes with multiple slices with the
% same slice-position. They could have been acquired in
% volume order (i.e. complete one volume first and then
% start on the next volume etc.) or they could have been
% acquired in slice order (i.e. do e.g. all diffusion 
% directions for one slice then go to next slice etc.).
% We'll use the .InstanceNumber and .SliceLocation to
% try and figure it out (don't you just love the DICOM
% standard?).
%

in = zeros(1,length(pvol{1}));
for i=1:length(pvol{1})
	in(i) = pvol{1}{i}.InstanceNumber;
end
[in,index] = sort(in);
if pvol{1}{index(1)}.SliceLocation == pvol{1}{index(2)}.SliceLocation % Slice First
	order = 'SF';
elseif pvol{1}{index(1)}.SliceLocation == pvol{1}{index(nsl(1)+1)}.SliceLocation % Volume first
	order = 'VF';
else  
	warning('Looks like there is something wrong with the conversion software.');
	vol = pvol;
	return
end

vol{1}{1} = [];
switch order
case 'SF'
	for i=1:length(pvol)
		z = zeros(1,length(pvol{i}));
		for j=1:length(pvol{i})
			z(j) = pvol{i}{j}.SliceLocation;
		end
		[z,zindex] = sort(z);
		tmp = pvol{i}(zindex);
		nvol = length(pvol{i})/nsl(i);
		for j=1:nsl(i)
			in = zeros(1,nvol);
			for k=1:nvol
				in = tmp{(j-1)*nvol+k}.InstanceNumber;
			end
			[in,index] = sort(in);
			tmp((j-1)*nvol+1:j*nvol) = tmp((j-1)*nvol+index);
		end
		for j=1:nvol
			if isempty(vol{1}{1})
				vol{end} = pvol{i}(j:nvol:length(pvol{i}));
			else
				vol{end+1} = pvol{i}(j:nvol:length(pvol{i}));
			end
		end
	end
case 'VF'
	for i=1:length(pvol) 
		in = zeros(1,length(pvol{i}));
		for j=1:length(pvol{i})
			in(j) = pvol{i}{j}.InstanceNumber;
		end
		[in,index] = sort(in);
		for j=1:(length(pvol{i}) / nsl(i))
			if isempty(vol{1}{1})
				vol{end} = pvol{i}(index((j-1)*nsl(i)+1:j*nsl(i)));
			else
				vol{end+1} = pvol{i}(index((j-1)*nsl(i)+1:j*nsl(i)));
			end
		end
	end
end

%
% Finally, cut and paste from Johns code to sort 
% the volumes acoording to his cunning scheme.
%

for j=1:length(vol),
	orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
	proj   = null(orient');
	if det([orient proj])<0, proj = -proj; end;

	z      = zeros(length(vol{j}),1);
	for i=1:length(vol{j}),
		z(i)  = vol{j}{i}.ImagePositionPatient*proj;
	end;
	[z,index] = sort(z);
	vol{j}    = vol{j}(index);
	if length(vol{j})>1,
		dist      = diff(z);
		if any(diff(z)==0)  % I give up!
			warning('Looks like there is something wrong with the conversion software.');
		end
		if sum((dist-mean(dist)).^2)/length(dist)>0.001,
			warning('Variable slice spacing');
		end
	end
end

return
%_______________________________________________________________________

%_______________________________________________________________________
function write_volume(hdr)

% Output filename
%-------------------------------------------------------------------
if checkfields(hdr{1},'SeriesNumber','AcquisitionNumber')
	fname = sprintf('s%s-%.4d-%.5d-%.6d.img', strip_unwanted(hdr{1}.PatientID),...
		hdr{1}.SeriesNumber, hdr{1}.AcquisitionNumber, hdr{1}.InstanceNumber);
else
	fname = sprintf('s%s-%.6d.img', strip_unwanted(hdr{1}.PatientID),hdr{1}.InstanceNumber);
end;
fname = fullfile(pwd,fname);

% Image dimensions
%-------------------------------------------------------------------
nc = hdr{1}.Columns;
nr = hdr{1}.Rows;

dim    = [nc nr length(hdr)];
dt     = [spm_type('int16') spm_platform('bigend')];

% Orientation information
%-------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];
orient           = reshape(hdr{1}.ImageOrientationPatient,[3 2]);
orient(:,3)      = null(orient');
if det(orient)<0, orient(:,3) = -orient(:,3); end;
if length(hdr)>1,
	z            = zeros(length(hdr),1);
	for i=1:length(hdr),
		z(i)     = hdr{i}.ImagePositionPatient*orient(:,3);
	end;
	z            = mean(diff(z));
else
	if checkfields(hdr{1},'SliceThickness'),
		z = hdr{1}.SliceThickness;
	else
		z = 1;
	end;
end;
vox              = [hdr{1}.PixelSpacing z];
pos              = hdr{1}.ImagePositionPatient';
dicom_to_patient = [orient*diag(vox) pos ; 0 0 0 1];
patient_to_tal   = diag([-1 -1 1 1]);
mat              = patient_to_tal*dicom_to_patient*analyze_to_dicom;

% Possibly useful information
%-------------------------------------------------------------------
if checkfields(hdr{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
	'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
	'AcquisitionDate'),
	tim = datevec(hdr{1}.AcquisitionTime/(24*60*60));
	descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g',...
		hdr{1}.MagneticFieldStrength, hdr{1}.MRAcquisitionType,...
		deblank(hdr{1}.ScanningSequence),...
		hdr{1}.RepetitionTime,hdr{1}.EchoTime,hdr{1}.FlipAngle,...
		datestr(hdr{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
	descrip = hdr{1}.Modality;
end;

if ~spm_flip_analyze_images,
	mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end; 

pinfo = [1 0 0]';
if isfield(hdr{1},'RescaleSlope') || isfield(hdr{1},'RescaleIntercept'),
	pinfo   = repmat(pinfo,1,length(hdr));
	bytepix = spm_type('int16','bits')/8;
	for i=1:length(hdr),
		if isfield(hdr{i},'RescaleSlope'),      pinfo(1,i) = hdr{i}.RescaleSlope;     end;
		if isfield(hdr{i},'RescaleIntercept'),  pinfo(2,i) = hdr{i}.RescaleIntercept; end;
		pinfo(3,i) = dim(1)*dim(2)*(i-1)*bytepix;
	end;
end;

% Write the image volume
%-------------------------------------------------------------------
spm_progress_bar('Init',length(hdr),['Writing ' fname], 'Planes written');
V = struct('fname',fname, 'dim',dim, 'dt',dt, 'pinfo',pinfo, 'mat',mat, 'descrip',descrip);
V = spm_create_vol(V);
for i=1:length(hdr),
	plane = read_image_data(hdr{i});

	if isfield(hdr{i},'RescaleSlope'),      plane = plane*hdr{i}.RescaleSlope;     end;
	if isfield(hdr{i},'RescaleIntercept'),  plane = plane+hdr{i}.RescaleIntercept; end;

	plane = fliplr(plane);
	if ~spm_flip_analyze_images, plane = flipud(plane); end;
	V     = spm_write_plane(V,plane,i);
	spm_progress_bar('Set',i);
end;
spm_close_vol(V);
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [images,guff] = select_tomographic_images(hdr)
images = {};
guff   = {};
for i=1:length(hdr),
	if ~checkfields(hdr{i},'Modality') || ~(strcmp(hdr{i}.Modality,'MR') ||...
		strcmp(hdr{i}.Modality,'PT') || strcmp(hdr{i}.Modality,'CT'))
		disp(['Cant find appropriate modality information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'StartOfPixelData','SamplesperPixel',...
		'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation'),
		disp(['Cant find "Image Pixel" information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient'),
		disp(['Cant find "Image Plane" information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'PatientID','InstanceNumber'),
	%elseif ~checkfields(hdr{i},'PatientID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
		disp(['Cant find suitable filename info for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	else
		images = {images{:},hdr{i}};
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [mosaic,standard] = select_mosaic_images(hdr)
mosaic   = {};
standard = {};
for i=1:length(hdr),
	if ~checkfields(hdr{i},'ImageType','CSAImageHeaderInfo') ||...
		isfield(hdr{i}.CSAImageHeaderInfo,'junk') ||...
		isempty(read_AcquisitionMatrixText(hdr{i})) ||...
		isempty(read_NumberOfImagesInMosaic(hdr{i}))
		standard = {standard{:},hdr{i}};
	else
		mosaic = {mosaic{:},hdr{i}};
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function ok = checkfields(hdr,varargin)
ok = 1;
for i=1:(nargin-1),
	if ~isfield(hdr,varargin{i}),
		ok = 0;
		break;
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function clean = strip_unwanted(dirty)
msk = find((dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
           (dirty>='0'&dirty<='9') | dirty=='_');
clean = dirty(msk);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function img = read_image_data(hdr)
img = [];
if hdr.SamplesperPixel ~= 1,
	warning([hdr.Filename ': SamplesperPixel = ' num2str(hdr.SamplesperPixel) ' - cant be an MRI']);
	return;
end;

prec = ['ubit' num2str(hdr.BitsAllocated) '=>' 'uint32'];

if isfield(hdr,'TransferSyntaxUID') && strcmp(hdr.TransferSyntaxUID,'1.2.840.10008.1.2.2') && strcmp(hdr.VROfPixelData,'OW'),
	fp = fopen(hdr.Filename,'r','ieee-be');
else
	fp = fopen(hdr.Filename,'r','ieee-le');
end;
if fp==-1,
	warning([hdr.Filename ': cant open file']);
	return;
end;

fseek(fp,hdr.StartOfPixelData,'bof');
img = fread(fp,hdr.Rows*hdr.Columns,prec);
fclose(fp);
if length(img)~=hdr.Rows*hdr.Columns,
	error([hdr.Filename ': cant read whole image']);
end;

img = bitshift(img,hdr.BitsStored-hdr.HighBit-1);

if hdr.PixelRepresentation,
	% Signed data - done this way because bitshift only
	% works with signed data.  Negative values are stored
	% as 2s complement.
	neg      = logical(bitshift(bitand(img,uint32(2^hdr.HighBit)),-hdr.HighBit));
	msk      = (2^hdr.HighBit - 1);
	img      = double(bitand(img,msk));
	img(neg) = img(neg)-2^(hdr.HighBit);
else
	% Unsigned data
	msk      = (2^(hdr.HighBit+1) - 1);
	img      = double(bitand(img,msk));
end;

img = reshape(img,hdr.Columns,hdr.Rows);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function nrm = read_SliceNormalVector(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'SliceNormalVector');
for i=1:3,
	nrm(i,1) = sscanf(val(i,:),'%g');
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function n = read_NumberOfImagesInMosaic(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'NumberOfImagesInMosaic');
n   = sscanf(val','%d');
if isempty(n), n=[]; end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function dim = read_AcquisitionMatrixText(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'AcquisitionMatrixText');
dim = sscanf(val','%d*%d')';
if length(dim)==1,
	dim = sscanf(val','%dp*%d')';
end;
if isempty(dim), dim=[]; end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str),
	if strcmp(deblank(str(i).name),name),
		for j=1:str(i).nitems,
			if  str(i).item(j).xx(1),
				val = {val{:} str(i).item(j).val};
			end;
		end;
		break;
	end;
end;
val = strvcat(val{:});
return;
%_______________________________________________________________________

