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
% %W% John Ashburner %E%
if nargin<2, opts = 'all'; end;

[images,guff]     = select_tomographic_images(hdr);
[mosaic,standard] = select_mosaic_images(images);

if (strcmp(opts,'all') | strcmp(opts,'mosaic')) & ~isempty(mosaic),
	convert_mosaic(mosaic);
end;
if (strcmp(opts,'all') | strcmp(opts,'standard')) & ~isempty(standard),
	convert_standard(standard);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function convert_raw(hdr)
disp('*** Ignoring Raw Data DICOM Files ***');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function convert_mosaic(hdr)
spm_progress_bar('Init',length(hdr),['Writing Mosaic'], 'Files written');

for i=1:length(hdr),

	% Output filename
	%-------------------------------------------------------------------
	fname = sprintf('f%s-%.4d-%.5d-%.6d.img', strip_unwanted(hdr{i}.StudyID),...
		hdr{i}.SeriesNumber, hdr{i}.AcquisitionNumber, hdr{i}.InstanceNumber);


	% Image dimensions and data
	%-------------------------------------------------------------------
	nc = hdr{i}.Columns;
	nr = hdr{i}.Rows;
	np = read_AcquisitionMatrixText(hdr{i});

	if rem(nc, np(1)) | rem(nr, np(2)),
		warning(sprintf('%s: %dx%d wont fit into %dx%d.',hdr{i}.Filename,...
			np(1), np(2), nc,nr));
		return;
	end;
	dim    = [np read_NumberOfImagesInMosaic(hdr{i})];
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
	dim = [size(volume) spm_type('int16')];

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

	analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1];

	vox    = [hdr{i}.PixelSpacing hdr{i}.SpacingBetweenSlices];
	pos    = hdr{i}.ImagePositionPatient';
	orient = reshape(hdr{i}.ImageOrientationPatient,[3 2]);
	orient(:,3) = null(orient');
	if det(orient)<0, orient(:,3) = -orient(:,3); end;

	% The image position vector is not correct. In dicom this vector points to
	% the upper left corner of the image. Perhaps it is unlucky that this is
	% calculated in the syngo software from the vector pointing to the center of
	% the slice (keep in mind: upper left slice) with the enlarged FoV.
	dicom_to_patient = [diag(vox)*orient pos ; 0 0 0 1];
	truepos          = dicom_to_patient *[(size(mosaic)-dim(1:2))/2 1 1]';
	dicom_to_patient = [diag(vox)*orient truepos(1:3) ; 0 0 0 1];
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

	V = struct('fname',fname,'dim',dim,'mat',mat,'descrip',descrip);
	V = spm_write_vol(V,volume);
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

vol{1}{1} = hdr{1};
for i=2:length(hdr),
	orient = reshape(hdr{i}.ImageOrientationPatient,[3 2]);
	xy1    = hdr{i}.ImagePositionPatient*orient;
	match  = 0;
	for j=1:length(vol),
		orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
		xy2    = vol{j}{1}.ImagePositionPatient*orient;
		dist2  = sum((xy1-xy2).^2);
		if        hdr{i}.SeriesNumber            == vol{j}{1}.SeriesNumber &...
		          hdr{i}.AcquisitionNumber       == vol{j}{1}.AcquisitionNumber &...
		   strcmp(hdr{i}.SeriesInstanceUID,         vol{j}{1}.SeriesInstanceUID) &...
		          hdr{i}.Rows                    == vol{j}{1}.Rows &...
		          hdr{i}.Columns                 == vol{j}{1}.Columns &...
		      all(hdr{i}.ImageOrientationPatient == vol{j}{1}.ImageOrientationPatient) &...
		      all(hdr{i}.PixelSpacing            == vol{j}{1}.PixelSpacing & dist2<0.00001) &...
	        (~isfield(hdr{i},'EchoNumber')   | ~isfield(vol{j}{1},'EchoNumber') | ...
		          hdr{i}.EchoNumber              == vol{j}{1}.EchoNumber),
			vol{j}{end+1} = hdr{i};
			match = 1;
			break;
		end;
	end;
	if ~match,
		vol{end+1}{1} = hdr{i};
	end;
end;

for j=1:length(vol),
	orient = reshape(hdr{j}.ImageOrientationPatient,[3 2]);
	proj   = null(orient');
	if det([orient proj])<0, proj = -proj; end;

	z         = zeros(length(vol{j}),1);
	for i=1:length(vol{j}),
		z(i)  = vol{j}{i}.ImagePositionPatient*proj;
	end;
	[z,index] = sort(z);
	vol{j}    = vol{j}(index);
	if length(vol{j})>1,
		dist      = diff(z);
		if any(diff(z)==0),
			warning('Looks like there is something wrong with the conversion software.');
		end;
		if sum((dist-mean(dist)).^2)/length(dist)>0.001,
			warning('Variable slice spacing');
		end;
	end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function write_volume(hdr)

% Output filename
%-------------------------------------------------------------------
fname = sprintf('s%s-%.4d-%.5d-%.6d.img', strip_unwanted(hdr{1}.StudyID),...
	hdr{1}.SeriesNumber, hdr{1}.AcquisitionNumber, hdr{1}.InstanceNumber);

% Image dimensions
%-------------------------------------------------------------------
nc = hdr{1}.Columns;
nr = hdr{1}.Rows;

dim    = [nc nr length(hdr) spm_type('int16')];

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

analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1];
orient           = reshape(hdr{1}.ImageOrientationPatient,[3 2]);
orient(:,3)      = null(orient');
if det(orient)<0, orient(:,3) = -orient(:,3); end;
if length(hdr)>1,
	z            = zeros(length(hdr),1);
	for i=1:length(hdr),
		z(i)     = hdr{i}.ImagePositionPatient*orient(:,3);
	end;
	z            = mean(diff(z));
else,
	if checkfields(hdr(1),'SliceThickness'),
		z = hdr{1}.SliceThickness;
	else,
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
else,
	descrip = hdr{1}.Modality;
end;

% Write the image volume
%-------------------------------------------------------------------
spm_progress_bar('Init',length(hdr),['Writing ' fname], 'Planes written');
V = struct('fname',fname, 'dim',dim, 'pinfo',[1 0 0]', 'mat',mat, 'descrip',descrip);
V = spm_create_vol(V);
for i=1:length(hdr),
	plane = read_image_data(hdr{i});
	plane = fliplr(plane);
	V     = spm_write_plane(V,plane,i);
	spm_progress_bar('Set',i);
end;
V = spm_close_vol(V);
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [images,guff] = select_tomographic_images(hdr)
images = {};
guff   = {};
for i=1:length(hdr),
	if ~checkfields(hdr{i},'Modality') | ~(strcmp(hdr{i}.Modality,'MR') |...
		strcmp(hdr{i}.Modality,'PT') | strcmp(hdr{i}.Modality,'CT'))
		disp(['Cant find appropriate modality information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'StartOfPixelData','SamplesperPixel',...
		'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation'),
		disp(['Cant find "Image Pixel" information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient'),
		disp(['Cant find "Image Plane" information for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	elseif ~checkfields(hdr{i},'StudyID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
		disp(['Cant find suitable filename info for "' hdr{i}.Filename '".']);
		guff = {guff{:},hdr{i}};
	else,
		images = {images{:},hdr{i}};
	end;
end;
return;
%_______________________________________________________________________

function [mosaic,standard] = select_mosaic_images(hdr)
mosaic   = {};
standard = {};
for i=1:length(hdr),
	if ~checkfields(hdr{i},'ImageType','CSAImageHeaderInfo') |...
		~strcmp(deblank(hdr{i}.ImageType),'ORIGINAL\PRIMARY\M\MOSAIC'),
		standard = {standard{:},hdr{i}};
	else,
		mosaic = {mosaic{:},hdr{i}};
	end;
end;
return;
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

if strcmp(hdr.TransferSyntaxUID,'1.2.840.10008.1.2.2') & strcmp(hdr.VROfPixelData,'OW'),
	fp = fopen(hdr.Filename,'r','ieee-be');
else,
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
	neg      = logical(bitand(img,uint32(2^hdr.HighBit)));
	msk      = (2^hdr.HighBit - 1);
	img      = double(bitand(img,msk));
	img(neg) = img(neg)-2^(hdr.HighBit);
else,
	% Unsigned data
	msk      = (2^(hdr.HighBit+1) - 1);
	img      = double(bitand(img,msk));
end;

img = reshape(img,hdr.Columns,hdr.Rows);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function nrm = read_SliceNormalVector(hdr)
val = get_numaris4_val(hdr.CSAImageHeaderInfo,'SliceNormalVector');
for i=1:3,
	nrm(i,1) = sscanf(val(i,:),'%g');
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function n = read_NumberOfImagesInMosaic(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(hdr.CSAImageHeaderInfo,'NumberOfImagesInMosaic');
n   = sscanf(val,'%d');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function dim = read_AcquisitionMatrixText(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(hdr.CSAImageHeaderInfo,'AcquisitionMatrixText');
dim = sscanf(val,'%d*%d')';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str),
	if strcmp(deblank(str(i).name),name),
		for j=1:str(i).nitems,
			val{j} = str(i).item(j).val;
		end;
		break;
	end;
end;
val = str2mat(val{:});
return;
%_______________________________________________________________________

%_______________________________________________________________________













