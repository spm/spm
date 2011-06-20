function out = spm_dicom_convert(hdr,opts,root_dir,format)
% Convert DICOM images into something that SPM can use
% FORMAT spm_dicom_convert(hdr,opts,root_dir,format)
% Inputs:
% hdr  - a cell array of DICOM headers from spm_dicom_headers
% opts - options
%        'all'      - all DICOM files [default]
%        'mosaic'   - the mosaic images
%        'standard' - standard DICOM files
%        'spect'    - SIEMENS Spectroscopy DICOMs (some formats only)
%                     This will write out a 5D NIFTI containing real and
%                     imaginary part of the spectroscopy time points at the
%                     position of spectroscopy voxel(s).
%        'raw'      - convert raw FIDs (not implemented)
% root_dir - 'flat' - do not produce file tree [default]
%            With all other options, files will be sorted into
%            directories according to their sequence/protocol names
%            'date_time'  - Place files under ./<StudyDate-StudyTime>
%            'patid'      - Place files under ./<PatID>
%            'patid_date' - Place files under ./<PatID-StudyDate>
%            'patname'    - Place files under ./<PatName>
%            'series'     - Place files in series folders, without
%                           creating patient folders
% format - output format
%          'img' Two file (hdr+img) NIfTI format [default]
%          'nii' Single file NIfTI format
%                All images will contain a single 3D dataset, 4D images
%                will not be created.
% Output:
% out - a struct with a single field .files. out.files contains a
%       cellstring with filenames of created files. If no files are
%       created, a cell with an empty string {''} is returned.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Jesper Andersson
% $Id: spm_dicom_convert.m 4368 2011-06-20 11:59:54Z john $


if nargin<2, opts     = 'all'; end
if nargin<3, root_dir = 'flat';end
if nargin<4, format   = 'img'; end

[images,other]    = select_tomographic_images(hdr);
[spect,guff]      = select_spectroscopy_images(other);
[mosaic,standard] = select_mosaic_images(images);

fmos = {};
fstd = {};
fspe = {};
if (strcmp(opts,'all') || strcmp(opts,'mosaic')) && ~isempty(mosaic),
    fmos = convert_mosaic(mosaic,root_dir,format);
end;
if (strcmp(opts,'all') || strcmp(opts,'standard')) && ~isempty(standard),
    fstd = convert_standard(standard,root_dir,format);
end;
if (strcmp(opts,'all') || strcmp(opts,'spect')) && ~isempty(spect),
    fspe = convert_spectroscopy(spect,root_dir,format);
end;

out.files = [fmos(:); fstd(:); fspe(:)];
if isempty(out.files)
    out.files = {''};
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function fnames = convert_mosaic(hdr,root_dir,format)
spm_progress_bar('Init',length(hdr),'Writing Mosaic', 'Files written');

fnames = cell(length(hdr),1);
for i=1:length(hdr),

    % Output filename
    %-------------------------------------------------------------------
    fnames{i} = getfilelocation(hdr{i},root_dir,'f',format);

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
    %   warning('%s: %dx%d wont fit into %dx%d.',hdr{i}.Filename,...
    %       np(1), np(2), nc,nr);
    %   return;
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
    dim = size(volume);
    dt  = determine_datatype(hdr{1});

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

    analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)-1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];

    vox    = [hdr{i}.PixelSpacing(:); hdr{i}.SpacingBetweenSlices];
    pos    = hdr{i}.ImagePositionPatient(:);
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

    if ~true, % LEFT-HANDED STORAGE
        mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
        volume = flipdim(volume,1);
    end;

    %if isfield(hdr{i},'RescaleSlope') && hdr{i}.RescaleSlope ~= 1,
    %   volume = volume*hdr{i}.RescaleSlope;
    %end;
    %if isfield(hdr{i},'RescaleIntercept') && hdr{i}.RescaleIntercept ~= 0,
    %   volume = volume + hdr{i}.RescaleIntercept;
    %end;
    %V = struct('fname',fname, 'dim',dim, 'dt',dt, 'mat',mat, 'descrip',descrip);
    %spm_write_vol(V,volume);

    % Note that data are no longer scaled by the maximum amount.
    % This may lead to rounding errors in smoothed data, but it
    % will get around other problems.
    RescaleSlope     = 1;
    RescaleIntercept = 0;
    if isfield(hdr{i},'RescaleSlope') && hdr{i}.RescaleSlope ~= 1,
        RescaleSlope     = hdr{i}.RescaleSlope;
    end;
    if isfield(hdr{i},'RescaleIntercept') && hdr{i}.RescaleIntercept ~= 0,
        RescaleIntercept = hdr{i}.RescaleIntercept;
    end;
    N      = nifti;
    N.dat  = file_array(fnames{i},dim,dt,0,RescaleSlope,RescaleIntercept);
    N.mat  = mat;
    N.mat0 = mat;
    N.mat_intent  = 'Scanner';
    N.mat0_intent = 'Scanner';
    N.descrip     = descrip;
    create(N);

    % Write the data unscaled
    dat           = N.dat;
    dat.scl_slope = [];
    dat.scl_inter = [];
    % write out volume at once - see spm_write_plane.m for performance comments
    dat(:,:,:) = volume;
    
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function fnames = convert_standard(hdr,root_dir,format)
hdr = sort_into_volumes(hdr);
fnames = cell(length(hdr),1);
for i=1:length(hdr),
    fnames{i} = write_volume(hdr{i},root_dir,format);
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
   %orient = reshape(hdr{i}.ImageOrientationPatient,[3 2]);
   %xy1    = hdr{i}.ImagePositionPatient(:)*orient;
    match  = 0;
    if isfield(hdr{i},'CSAImageHeaderInfo') && isfield(hdr{i}.CSAImageHeaderInfo,'name')
        ice1 = sscanf( ...
            strrep(get_numaris4_val(hdr{i}.CSAImageHeaderInfo,'ICE_Dims'), ...
            'X', '-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
        dimsel = logical([1 1 1 1 1 1 0 0 1]);
    else
        ice1 = [];
    end;
    for j=1:length(vol),
       %orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
       %xy2    = vol{j}{1}.ImagePositionPatient(:)*orient;
        
        % This line is a fudge because of some problematic data that Bogdan,
        % Cynthia and Stefan were trying to convert.  I hope it won't cause
        % problems for others -JA
        % dist2  = sum((xy1-xy2).^2);
        dist2 = 0;
        
        if strcmp(hdr{i}.Modality,'CT') && ...
                strcmp(vol{j}{1}.Modality,'CT') % Our CT seems to have shears in slice positions
            dist2 = 0;
        end;
        if ~isempty(ice1) && isfield(vol{j}{1},'CSAImageHeaderInfo') && isfield(vol{j}{1}.CSAImageHeaderInfo(1),'name')
            % Replace 'X' in ICE_Dims by '-1'
            ice2 = sscanf( ...
                strrep(get_numaris4_val(vol{j}{1}.CSAImageHeaderInfo,'ICE_Dims'), ...
                'X', '-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
            if ~isempty(ice2)
                identical_ice_dims=all(ice1(dimsel)==ice2(dimsel));
            else
                identical_ice_dims = 0; % have ice1 but not ice2, ->
                % something must be different
            end,
        else
            identical_ice_dims = 1; % No way of knowing if there is no CSAImageHeaderInfo
        end;
        try
            match = hdr{i}.SeriesNumber            == vol{j}{1}.SeriesNumber &&...
                hdr{i}.Rows                        == vol{j}{1}.Rows &&...
                hdr{i}.Columns                     == vol{j}{1}.Columns &&...
                sum((hdr{i}.ImageOrientationPatient - vol{j}{1}.ImageOrientationPatient).^2)<1e-4 &&...
                sum((hdr{i}.PixelSpacing            - vol{j}{1}.PixelSpacing).^2)<1e-4 && ...
                identical_ice_dims && dist2<1e-3;
            %if (hdr{i}.AcquisitionNumber ~= hdr{i}.InstanceNumber) || ...
            %   (vol{j}{1}.AcquisitionNumber ~= vol{j}{1}.InstanceNumber)
            %    match = match && (hdr{i}.AcquisitionNumber == vol{j}{1}.AcquisitionNumber)
            %end;
            % For raw image data, tell apart real/complex or phase/magnitude
            if isfield(hdr{i},'ImageType') && isfield(vol{j}{1}, 'ImageType')
                match = match && strcmp(hdr{i}.ImageType, vol{j}{1}.ImageType);
            end;
            if isfield(hdr{i},'SequenceName') && isfield(vol{j}{1}, 'SequenceName')
                match = match && strcmp(hdr{i}.SequenceName,vol{j}{1}.SequenceName);
            end;
            if isfield(hdr{i},'SeriesInstanceUID') && isfield(vol{j}{1}, 'SeriesInstanceUID')
                match = match && strcmp(hdr{i}.SeriesInstanceUID,vol{j}{1}.SeriesInstanceUID);
            end;
            if isfield(hdr{i},'EchoNumbers')  && isfield(vol{j}{1}, 'EchoNumbers')
                match = match && hdr{i}.EchoNumbers == vol{j}{1}.EchoNumbers;
            end;
        catch
            match = 0;
        end
        if match
            vol{j}{end+1} = hdr{i};
            break;
        end;
    end;
    if ~match,
        vol{end+1}{1} = hdr{i};
    end;
end;

%dcm = vol;
%save('dicom_headers.mat','dcm');

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
        z(i)  = vol{j}{i}.ImagePositionPatient(:)'*proj;
    end;
    [z,index] = sort(z);
    vol{j}    = vol{j}(index);
    if length(vol{j})>1,
        % dist      = diff(z);
        if any(diff(z)==0)
            tmp = sort_into_vols_again(vol{j});
            vol{j} = tmp{1};
            vol2 = {vol2{:} tmp{2:end}};
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
            z(i)  = vol{j}{i}.ImagePositionPatient(:)'*proj;
        end;
        dist = diff(sort(z));
        if sum((dist-mean(dist)).^2)/length(dist)>1e-4,
            fprintf('***************************************************\n');
            fprintf('* VARIABLE SLICE SPACING                          *\n');
            fprintf('* This may be due to missing DICOM files.         *\n');
            if checkfields(vol{j}{1},'PatientID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
                fprintf('*    %s / %d / %d / %d \n',...
                    deblank(vol{j}{1}.PatientID), vol{j}{1}.SeriesNumber, ...
                    vol{j}{1}.AcquisitionNumber, vol{j}{1}.InstanceNumber);
                fprintf('*                                                 *\n');
            end;
            fprintf('*  %20.4g                           *\n', dist);
            fprintf('***************************************************\n');
        end;
    end;
end;
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
    z(i)  = volj{i}.ImagePositionPatient(:)'*proj;
    t(i)  = volj{i}.InstanceNumber;
end;
% msg = 0;
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
    index = [true ; diff(t)~=0];
    t     = t(index);
    z     = z(index);
    d     = d(index);
    volj  = volj(index);
end;

%if any(diff(sort(t))~=1), msg = 1; end;
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
function fname = write_volume(hdr,root_dir,format)

% Output filename
%-------------------------------------------------------------------
fname = getfilelocation(hdr{1}, root_dir,'s',format);

% Image dimensions
%-------------------------------------------------------------------
nc = hdr{1}.Columns;
nr = hdr{1}.Rows;

dim    = [nc nr length(hdr)];
dt     = determine_datatype(hdr{1});

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

analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
patient_to_tal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions

R  = [reshape(hdr{1}.ImageOrientationPatient,3,2)*diag(hdr{1}.PixelSpacing); 0 0];
x1 = [1;1;1;1];
y1 = [hdr{1}.ImagePositionPatient(:); 1];

if length(hdr)>1,
    x2 = [1;1;dim(3); 1];
    y2 = [hdr{end}.ImagePositionPatient(:); 1];
else
    orient           = reshape(hdr{1}.ImageOrientationPatient,[3 2]);
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end;
    if checkfields(hdr{1},'SliceThickness'),
        z = hdr{1}.SliceThickness;
    else
        z = 1;
    end
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
dicom_to_patient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat              = patient_to_tal*dicom_to_patient*analyze_to_dicom;

% Possibly useful information
%-------------------------------------------------------------------
if checkfields(hdr{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate'),
    if isfield(hdr{1},'ScanOptions'),
        ScanOptions = hdr{1}.ScanOptions;
    else
        ScanOptions = 'no';
    end
    tim = datevec(hdr{1}.AcquisitionTime/(24*60*60));
    descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg/SO=%s %s %d:%d:%.5g',...
        hdr{1}.MagneticFieldStrength, hdr{1}.MRAcquisitionType,...
        deblank(hdr{1}.ScanningSequence),...
        hdr{1}.RepetitionTime,hdr{1}.EchoTime,hdr{1}.FlipAngle,...
        ScanOptions,...
        datestr(hdr{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
    descrip = hdr{1}.Modality;
end;

if ~true, % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end;

% Write the image volume
%-------------------------------------------------------------------
spm_progress_bar('Init',length(hdr),['Writing ' fname], 'Planes written');
N      = nifti;
pinfos = [ones(length(hdr),1) zeros(length(hdr),1)];
for i=1:length(hdr)
    if isfield(hdr{i},'RescaleSlope'),     pinfos(i,1) = hdr{i}.RescaleSlope;     end 
    if isfield(hdr{i},'RescaleIntercept'), pinfos(i,2) = hdr{i}.RescaleIntercept; end
end

if any(any(diff(pinfos,1))),
    % Ensure random numbers are reproducible (see later)
    % when intensities are dithered to prevent aliasing effects.
    rand('state',0);
end

volume = zeros(dim);
for i=1:length(hdr),
    plane = read_image_data(hdr{i});

    if any(any(diff(pinfos,1))),
        % This is to prevent aliasing effects in any subsequent histograms
        % of the data (eg for mutual information coregistration).
        % It's a bit inelegant, but probably necessary for when slices are
        % individually rescaled.
        plane = double(plane) + rand(size(plane)) - 0.5;
    end

    if pinfos(i,1)~=1, plane = plane*pinfos(i,1); end;
    if pinfos(i,2)~=0, plane = plane+pinfos(i,2); end;

    plane = fliplr(plane);
    if ~true, plane = flipud(plane); end; % LEFT-HANDED STORAGE
    volume(:,:,i) = plane;
    spm_progress_bar('Set',i);
end;

if ~any(any(diff(pinfos,1))),
    % Same slopes and intercepts for all slices
    pinfo = pinfos(1,:);
else
    % Variable slopes and intercept (maybe PET/SPECT)
    mx = max(volume(:));
    mn = min(volume(:));

    %%  Slope and Intercept
    %%  32767*pinfo(1) + pinfo(2) = mx
    %% -32768*pinfo(1) + pinfo(2) = mn
    % pinfo = ([32767 1; -32768 1]\[mx; mn])';

    % Slope only
    dt    = 'int16-be';
    pinfo = [max(mx/32767,-mn/32768) 0];
end

N.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = descrip;
create(N);
N.dat(:,:,:) = volume;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function fnames = convert_spectroscopy(hdr,root_dir,format)
fnames = cell(length(hdr),1);
for i=1:length(hdr),
    fnames{i} = write_spectroscopy_volume(hdr(i),root_dir,format);
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function fname = write_spectroscopy_volume(hdr,root_dir,format)
% Output filename
%-------------------------------------------------------------------
fname = getfilelocation(hdr{1}, root_dir,'S',format);

% guess private field to use
if isfield(hdr{1}, 'Private_0029_1210')
    privdat = hdr{1}.Private_0029_1210;
elseif isfield(hdr{1}, 'Private_0029_1110')
    privdat = hdr{1}.Private_0029_1110;
else
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end

% Image dimensions
%-------------------------------------------------------------------
nc = get_numaris4_numval(privdat,'Columns');
nr = get_numaris4_numval(privdat,'Rows');
% Guess number of timepoints in file - don't know whether this should be
% 'DataPointRows'-by-'DataPointColumns' or 'SpectroscopyAcquisitionDataColumns'
ntp = get_numaris4_numval(privdat,'DataPointRows')*get_numaris4_numval(privdat,'DataPointColumns');

dim    = [nc nr numel(hdr) 2 ntp];
dt     = spm_type('float32'); % Fixed datatype

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

analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
patient_to_tal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions
shift_vx         = [eye(4,3) [.5; .5; 0; 1]];

orient           = reshape(get_numaris4_numval(privdat,...
                                               'ImageOrientationPatient'),[3 2]);
ps               = get_numaris4_numval(privdat,'PixelSpacing');
if nc*nr == 1
    % Single Voxel Spectroscopy (based on the following information from SIEMENS)
    %---------------------------------------------------------------
    % NOTE: Internally the position vector of the CSI matrix shows to the outer border
    % of the first voxel. Therefore the position vector has to be corrected.
    % (Note: The convention of Siemens spectroscopy raw data is in contrast to the
    %  DICOM standard where the position vector points to the center of the first voxel.)
    %---------------------------------------------------------------
    % SIEMENS decides which definition to use based on the contents of the
    % 'PixelSpacing' internal header field. If it has non-zero values,
    % assume DICOM convention. If any value is zero, assume SIEMENS
    % internal convention for this direction.
    % Note that in SIEMENS code, there is a shift when PixelSpacing is
    % zero. Here, the shift seems to be necessary when PixelSpacing is
    % non-zero. This may indicate more fundamental problems with
    % orientation decoding.
    if ps(1) == 0 % row
        ps(1) = get_numaris4_numval(privdat,...
                                    'VoiPhaseFoV');
        shift_vx(1,4) = 0;
    end
    if ps(2) == 0 % col
        ps(2) = get_numaris4_numval(privdat,...
                                    'VoiReadoutFoV');
        shift_vx(2,4) = 0;
    end
end
pos = get_numaris4_numval(privdat,'ImagePositionPatient');
% for some reason, pixel spacing needs to be swapped
R  = [orient*diag(ps([2 1])); 0 0];
x1 = [1;1;1;1];
y1 = [pos; 1];

if length(hdr)>1,
    error('spm_dicom_convert:spectroscopy',...
        'Don''t know how to handle multislice spectroscopy data.');
else
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end;
    try
        z = get_numaris4_numval(privdat,...
            'VoiThickness');
    catch
        try
            z = get_numaris4_numval(privdat,...
                'SliceThickness');
        catch
            z = 1;
        end
    end;
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
dicom_to_patient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat              = patient_to_tal*dicom_to_patient*shift_vx*analyze_to_dicom;

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

if ~true, % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end;

% Write the image volume
%-------------------------------------------------------------------
N      = nifti;
pinfo  = [1 0];
if isfield(hdr{1},'RescaleSlope'),      pinfo(1) = hdr{1}.RescaleSlope;     end;
if isfield(hdr{1},'RescaleIntercept'),  pinfo(2) = hdr{1}.RescaleIntercept; end;
N.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = descrip;
N.extras      = struct('MagneticFieldStrength',...
                       get_numaris4_numval(privdat,'MagneticFieldStrength'),...
                       'TransmitterReferenceAmplitude',...
                       get_numaris4_numval(privdat,'TransmitterReferenceAmplitude'));
create(N);

% Read data, swap dimensions
data = permute(reshape(read_spect_data(hdr{1},privdat),dim([4 5 1 2 3])), ...
                [3 4 5 1 2]);
% plane = fliplr(plane);

N.dat(:,:,:,:,:) = data;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [images,guff] = select_tomographic_images(hdr)
images = {};
guff   = {};
for i=1:length(hdr),
    if ~checkfields(hdr{i},'Modality') || ~(strcmp(hdr{i}.Modality,'MR') ||...
            strcmp(hdr{i}.Modality,'PT') || strcmp(hdr{i}.Modality,'CT'))
        if checkfields(hdr{i},'Modality'),
            fprintf('File "%s" can not be converted because it is of type "%s", which is not MRI, CT or PET.\n', hdr{i}.Filename, hdr{i}.Modality);
        else
            fprintf('File "%s" can not be converted because it does not encode an image.\n', hdr{i}.Filename);
        end
        guff = [guff(:)',hdr(i)];
    elseif ~checkfields(hdr{i},'StartOfPixelData','SamplesperPixel',...
            'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation'),
        disp(['Cant find "Image Pixel" information for "' hdr{i}.Filename '".']);
        guff = [guff(:)',hdr(i)];
   %elseif isfield(hdr{i},'Private_2001_105f'),
   %    % This field corresponds to: > Stack Sequence 2001,105F SQ VNAP, COPY
   %    % http://www.medical.philips.com/main/company/connectivity/mri/index.html
   %    % No documentation about this private field is yet available.
   %    disp('Cant yet convert Phillips Intera DICOM.');
   %    guff = {guff{:},hdr{i}};
    elseif ~(checkfields(hdr{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient')||isfield(hdr{i},'Private_0029_1110')||isfield(hdr{i},'Private_0029_1210')),
        disp(['Cant find "Image Plane" information for "' hdr{i}.Filename '".']);
        guff = [guff(:)',hdr(i)];
    elseif ~checkfields(hdr{i},'PatientID','SeriesNumber','AcquisitionNumber','InstanceNumber'),
       %disp(['Cant find suitable filename info for "' hdr{i}.Filename '".']);
        if ~isfield(hdr{i},'SeriesNumber')
            disp('Setting SeriesNumber to 1');
            hdr{i}.SeriesNumber=1;
            images = [images(:)',hdr(i)];
        end;
        if ~isfield(hdr{i},'AcquisitionNumber')
            if isfield(hdr{i},'Manufacturer') && ~isempty(strfind(upper(hdr{1}.Manufacturer), 'PHILIPS'))
                % WHY DO PHILIPS DO THINGS LIKE THIS????
                if isfield(hdr{i},'InstanceNumber')
                     hdr{i}.AcquisitionNumber = hdr{i}.InstanceNumber;
                else
                     disp('Setting AcquisitionNumber to 1');
                     hdr{i}.AcquisitionNumber=1;
                end
             else
                disp('Setting AcquisitionNumber to 1');
                hdr{i}.AcquisitionNumber=1;
             end
            images = [images(:)',hdr(i)];
        end;
        if ~isfield(hdr{i},'InstanceNumber')
            disp('Setting InstanceNumber to 1');
            hdr{i}.InstanceNumber=1;
            images = [images(:)',hdr(i)];
        end;
    else
        images = [images(:)',hdr(i)];
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
            isempty(read_NumberOfImagesInMosaic(hdr{i})) ||...
            read_NumberOfImagesInMosaic(hdr{i}) == 0
        % NumberOfImagesInMosaic seems to be set to zero for pseudo images
        % containing e.g. online-fMRI design matrices, don't treat them as
        % mosaics
        standard = {standard{:}, hdr{i}};
    else
        mosaic = {mosaic{:},hdr{i}};
    end;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [spect,images] = select_spectroscopy_images(hdr)
spectsel = zeros(1,numel(hdr));
for i=1:length(hdr),
    if isfield(hdr{i},'SOPClassUID')
        spectsel(i) = strcmp(hdr{i}.SOPClassUID,'1.3.12.2.1107.5.9.1');
    end;
end;
spect  = hdr(logical(spectsel));
images = hdr(~logical(spectsel));
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
msk = (dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
      (dirty>='0'&dirty<='9') | dirty=='_';
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

if isfield(hdr,'TransferSyntaxUID')
    switch(hdr.TransferSyntaxUID)
    case {'1.2.840.10008.1.2.4.50','1.2.840.10008.1.2.4.51','1.2.840.10008.1.2.4.70',...
          '1.2.840.10008.1.2.4.80','1.2.840.10008.1.2.4.90','1.2.840.10008.1.2.4.91'},
        % try to read PixelData as JPEG image - offset is just a guess
        offset = 16;
        fseek(fp,hdr.StartOfPixelData+offset,'bof');
        img = fread(fp,Inf,'*uint8');
        % save PixelData into temp file - imread and its subroutines can only
        % read from file, not from memory
        tfile = tempname;
        tfp = fopen(tfile,'w+');
        fwrite(tfp,img,'uint8');
        fclose(tfp);
        % read decompressed data, transpose to match DICOM row/column order
        img = imread(tfile)';
        delete(tfile);
    otherwise
        fseek(fp,hdr.StartOfPixelData,'bof');
        img = fread(fp,hdr.Rows*hdr.Columns,prec);
    end
else
    fseek(fp,hdr.StartOfPixelData,'bof');
    img = fread(fp,hdr.Rows*hdr.Columns,prec);
end
fclose(fp);
if numel(img)~=hdr.Rows*hdr.Columns,
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
function img = read_spect_data(hdr,privdat)
% Guess number of timepoints in file - don't know whether this should be
% 'DataPointRows'-by-'DataPointColumns' or 'SpectroscopyAcquisitionDataColumns'
ntp = get_numaris4_numval(privdat,'DataPointRows')*get_numaris4_numval(privdat,'DataPointColumns');
% Data is stored as complex float32 values, timepoint by timepoint, voxel
% by voxel. Reshaping is done in write_spectroscopy_volume.
if ntp*2*4 ~= hdr.SizeOfCSAData
    warning([hdr.Filename,': Data size mismatch.']);
end
fp = fopen(hdr.Filename,'r','ieee-le');
fseek(fp,hdr.StartOfCSAData,'bof');
img = fread(fp,2*ntp,'float32');
fclose(fp);
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

%_______________________________________________________________________
function val = get_numaris4_numval(str,name)
val1 = get_numaris4_val(str,name);
val  = zeros(size(val1,1),1);
for k = 1:size(val1,1)
    val(k)=str2num(val1(k,:));
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________

function fname = getfilelocation(hdr,root_dir,prefix,format)

if nargin < 3
    prefix = 'f';
end;

if strncmp(root_dir,'ice',3)
    root_dir = root_dir(4:end);
    imtype = textscan(hdr.ImageType,'%s','delimiter','\\');
    try
        imtype = imtype{1}{3};
    catch
        imtype = '';
    end;
    prefix = [prefix imtype get_numaris4_val(hdr.CSAImageHeaderInfo,'ICE_Dims')];
end;

if strcmp(root_dir, 'flat')
    % Standard SPM file conversion
    %-------------------------------------------------------------------
    if checkfields(hdr,'SeriesNumber','AcquisitionNumber')
        if checkfields(hdr,'EchoNumbers')
            fname = sprintf('%s%s-%.4d-%.5d-%.6d-%.2d.%s', prefix, strip_unwanted(hdr.PatientID),...
                hdr.SeriesNumber, hdr.AcquisitionNumber, hdr.InstanceNumber,...
                hdr.EchoNumbers, format);
        else
            fname = sprintf('%s%s-%.4d-%.5d-%.6d.%s', prefix, strip_unwanted(hdr.PatientID),...
                hdr.SeriesNumber, hdr.AcquisitionNumber, ...
                hdr.InstanceNumber, format);
        end;
    else
        fname = sprintf('%s%s-%.6d.%s',prefix, ...
            strip_unwanted(hdr.PatientID),hdr.InstanceNumber, format);
    end;

    fname = fullfile(pwd,fname);
    return;
end;

% more fancy stuff - sort images into subdirectories
if ~isfield(hdr,'ProtocolName')
    if isfield(hdr,'SequenceName')
        hdr.ProtocolName = hdr.SequenceName;
    else
        hdr.ProtocolName='unknown';
    end;
end;
if ~isfield(hdr,'SeriesDescription')
    hdr.SeriesDescription = 'unknown';
end;
if ~isfield(hdr,'EchoNumbers')
    hdr.EchoNumbers = 0;
end;

m = sprintf('%02d', floor(rem(hdr.StudyTime/60,60)));
h = sprintf('%02d', floor(hdr.StudyTime/3600));
studydate = sprintf('%s_%s-%s', datestr(hdr.StudyDate,'yyyy-mm-dd'), ...
    h,m);
switch root_dir
    case {'date_time','series'}
    id = studydate;
    case {'patid', 'patid_date', 'patname'},
    id = strip_unwanted(hdr.PatientID);
end;
serdes = strrep(strip_unwanted(hdr.SeriesDescription),...
    strip_unwanted(hdr.ProtocolName),'');
protname = sprintf('%s%s_%.4d',strip_unwanted(hdr.ProtocolName), ...
    serdes, hdr.SeriesNumber);
switch root_dir
    case 'date_time',
        dname = fullfile(pwd, id, protname);
    case 'patid',
        dname = fullfile(pwd, id, protname);
    case 'patid_date',
        dname = fullfile(pwd, id, studydate, protname);
    case 'patname',
        dname = fullfile(pwd, strip_unwanted(hdr.PatientsName), ...
            id, protname);
    case 'series',
        dname = fullfile(pwd, protname);
    otherwise
        error('unknown file root specification');
end;
if ~exist(dname,'dir'),
    mkdir_rec(dname);
end;

% some non-product sequences on SIEMENS scanners seem to have problems
% with image numbering in MOSAICs - doublettes, unreliable ordering
% etc. To distinguish, always include Acquisition time in image name
sa = sprintf('%02d', floor(rem(hdr.AcquisitionTime,60)));
ma = sprintf('%02d', floor(rem(hdr.AcquisitionTime/60,60)));
ha = sprintf('%02d', floor(hdr.AcquisitionTime/3600));
fname = sprintf('%s%s-%s%s%s-%.5d-%.5d-%d.%s', prefix, id, ha, ma, sa, ...
        hdr.AcquisitionNumber,hdr.InstanceNumber, ...
        hdr.EchoNumbers,format);

fname = fullfile(dname, fname);

%_______________________________________________________________________

%_______________________________________________________________________

function suc = mkdir_rec(str)
% works on full pathnames only
opwd=pwd;
if str(end) ~= filesep, str = [str filesep];end;
pos = strfind(str,filesep);
suc = zeros(1,length(pos));
for g=2:length(pos)
    if ~exist(str(1:pos(g)-1),'dir'),
        cd(str(1:pos(g-1)-1));
        suc(g) = mkdir(str(pos(g-1)+1:pos(g)-1));
    end;
end;
cd(opwd);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function ret = read_ascconv(hdr)
% In SIEMENS data, there is an ASCII text section with
% additional information items. This section starts with a code
% ### ASCCONV BEGIN ###
% and ends with
% ### ASCCONV END ###
% It is read by spm_dicom_headers into an entry 'MrProtocol' in
% CSASeriesHeaderInfo or into an entry 'MrPhoenixProtocol' in
% Private_0029_1110 or Private_0029_1120.
% The additional items are assignments in C syntax, here they are just
% translated according to
% [] -> ()
% "  -> '
% 0xX -> hex2dec('X')
% and collected in a struct.
ret=struct;

% get ascconv data
if isfield(hdr, 'Private_0029_1110')
    X = get_numaris4_val(hdr.Private_0029_1110,'MrPhoenixProtocol');
elseif isfield(hdr, 'Private_0029_1120')
    X = get_numaris4_val(hdr.Private_0029_1120,'MrPhoenixProtocol');
else
    X=get_numaris4_val(hdr.CSASeriesHeaderInfo,'MrProtocol');
end

ascstart = strfind(X,'### ASCCONV BEGIN ###');
ascend = strfind(X,'### ASCCONV END ###');

if ~isempty(ascstart) && ~isempty(ascend)
    tokens = textscan(char(X((ascstart+22):(ascend-1))),'%s', ...
        'delimiter',char(10));
    tokens{1}=regexprep(tokens{1},{'\[([0-9]*)\]','"(.*)"','0x([0-9a-fA-F]*)'},{'($1+1)','''$1''','hex2dec(''$1'')'});
    % If everything would evaluate correctly, we could use
    % eval(sprintf('ret.%s;\n',tokens{1}{:}));
    for k = 1:numel(tokens{1})
        try
            eval(['ret.' tokens{1}{k} ';']);
        catch
            disp(['AscConv: Error evaluating ''ret.' tokens{1}{k} ''';']);
        end;
    end;
end;
%_______________________________________________________________________

%_______________________________________________________________________
function dt = determine_datatype(hdr)
% Determine what datatype to use for NIfTI images
be = spm_platform('bigend');
if hdr.HighBit>16
    if hdr.PixelRepresentation
        dt  = [spm_type( 'int32') be];
    else
        dt  = [spm_type('uint32') be];
    end
else
    if hdr.PixelRepresentation || hdr.HighBit<=15
        dt  = [spm_type( 'int16') be];
    else
        dt  = [spm_type('uint16') be];
    end
end

