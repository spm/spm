function out = spm_dicom_convert(hdr,opts,root_dir,format,out_dir)
% Convert DICOM images into something that SPM can use (e.g. NIfTI)
% FORMAT spm_dicom_convert(hdr,opts,root_dir,format,out_dir)
% Inputs:
% hdr      - a cell array of DICOM headers from spm_dicom_headers
% opts     - options:
%              'all'      - all DICOM files [default]
%              'mosaic'   - the mosaic images
%              'standard' - standard DICOM files
%              'spect'    - SIEMENS Spectroscopy DICOMs (some formats only)
%                           This will write out a 5D NIFTI containing real
%                           and imaginary part of the spectroscopy time 
%                           points at the position of spectroscopy voxel(s).
%              'raw'      - convert raw FIDs (not implemented)
% root_dir - 'flat'       - do not produce file tree [default]
%              With all other options, files will be sorted into
%              directories according to their sequence/protocol names:
%            'date_time'  - Place files under ./<StudyDate-StudyTime>
%            'patid'      - Place files under ./<PatID>
%            'patid_date' - Place files under ./<PatID-StudyDate>
%            'patname'    - Place files under ./<PatName>
%            'series'     - Place files in series folders, without
%                           creating patient folders
% format   - output format:
%              'nii'      - Single file NIfTI format [default]
%              'img'      - Two file (hdr+img) NIfTI format
%            All images will contain a single 3D dataset, 4D images will
%            not be created.
% out_dir  - output directory name.
%
% Output:
% out      - a struct with a single field .files. out.files contains a
%            cellstring with filenames of created files. If no files are
%            created, a cell with an empty string {''} is returned.
%__________________________________________________________________________
% Copyright (C) 2002-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dicom_convert.m 6899 2016-10-07 08:23:34Z volkmar $


%-Input parameters
%--------------------------------------------------------------------------
if nargin<2, opts     = 'all';  end
if nargin<3, root_dir = 'flat'; end
if nargin<4, format   = spm_get_defaults('images.format'); end
if nargin<5, out_dir  = pwd;    end

%-Select files
%--------------------------------------------------------------------------
[images, other]     = select_tomographic_images(hdr);
[multiframe,images] = select_multiframe(images);
[spect, guff]       = select_spectroscopy_images(other);
[mosaic, standard]  = select_mosaic_images(images);
[standard, guff]    = select_last_guff(standard, guff);

if ~isempty(guff)
    warning('spm:dicom','%d files could not be converted from DICOM.', numel(guff));
end

%-Convert files
%--------------------------------------------------------------------------
fmos = {};
fstd = {};
fspe = {};
if (strcmp(opts,'all') || strcmp(opts,'mosaic')) && ~isempty(mosaic)
    fmos = convert_mosaic(mosaic,root_dir,format,out_dir);
end
if (strcmp(opts,'all') || strcmp(opts,'standard')) && ~isempty(standard)
    fstd = convert_standard(standard,root_dir,format,out_dir);
end
if (strcmp(opts,'all') || strcmp(opts,'spect')) && ~isempty(spect)
    fspe = convert_spectroscopy(spect,root_dir,format,out_dir);
end
if (strcmp(opts,'all') || strcmp(opts,'multiframe')) && ~isempty(multiframe)
    fspe = convert_multiframes(multiframe,root_dir,format,out_dir);
end
out.files = [fmos(:); fstd(:); fspe(:)];
if isempty(out.files)
    out.files = {''};
end


%==========================================================================
% function fnames = convert_mosaic(hdr,root_dir,format,out_dir)
%==========================================================================
function fnames = convert_mosaic(hdr,root_dir,format,out_dir)
spm_progress_bar('Init',length(hdr),'Writing Mosaic', 'Files written');

fnames = cell(length(hdr),1);
for i=1:length(hdr)

    % Output filename
    %----------------------------------------------------------------------
    fnames{i} = getfilelocation(hdr{i},root_dir,'f',format,out_dir);

    % Image dimensions and data
    %----------------------------------------------------------------------
    nc = hdr{i}.Columns;
    nr = hdr{i}.Rows;

    dim      = [0 0 0];
    dim(3)   = read_NumberOfImagesInMosaic(hdr{i});
    np       = [nc nr]/ceil(sqrt(dim(3)));
    dim(1:2) = np;
    if ~all(np==floor(np))
        warning('spm:dicom','%s: dimension problem [Num Images=%d, Num Cols=%d, Num Rows=%d].',...
            hdr{i}.Filename,dim(3), nc,nr);
        continue
    end

    % Apparently, this is not the right way of doing it.
    %np = read_AcquisitionMatrixText(hdr{i});
    %if rem(nc, np(1)) || rem(nr, np(2)),
    %   warning('spm:dicom','%s: %dx%d wont fit into %dx%d.',hdr{i}.Filename,...
    %       np(1), np(2), nc,nr);
    %   return;
    %end;
    %dim    = [np read_NumberOfImagesInMosaic(hdr{i})];

    mosaic = read_image_data(hdr{i});
    volume = zeros(dim);
    snnz   = zeros(dim(3), 1);
    for j=1:dim(3)
        img = mosaic((1:np(1))+np(1)*rem(j-1,nc/np(1)), (np(2):-1:1)+np(2)*floor((j-1)/(nc/np(1))));
        snnz(j) = nnz(img) > 0;
        volume(:,:,j) = img;
    end
    d3 = find(snnz, 1, 'last');
    if ~isempty(d3)
        dim(3) = d3;
        volume = volume(:,:,1:dim(3));
    end
    dt  = determine_datatype(hdr{1});

    % Orientation information
    %----------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    SliceNormalVector = read_SliceNormalVector(hdr{i});
    if det([reshape(hdr{i}.ImageOrientationPatient,[3 2]) SliceNormalVector(:)])<0
        volume = volume(:,:,end:-1:1);
        mat    = mat*[eye(3) [0 0 -(dim(3)-1)]'; 0 0 0 1];
    end


    % Possibly useful information
    %----------------------------------------------------------------------
    if checkfields(hdr{i},'AcquisitionTime','MagneticFieldStrength',...
            'MRAcquisitionType','ScanningSequence','RepetitionTime',...
            'EchoTime','FlipAngle','AcquisitionDate')
        tim = datevec(hdr{i}.AcquisitionTime/(24*60*60));
        descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g Mosaic',...
            hdr{i}.MagneticFieldStrength, hdr{i}.MRAcquisitionType,...
            deblank(hdr{i}.ScanningSequence),...
            hdr{i}.RepetitionTime,hdr{i}.EchoTime,hdr{i}.FlipAngle,...
            datestr(hdr{i}.AcquisitionDate),tim(4),tim(5),tim(6));
    else
        descrip = hdr{1}.Modality;
    end

    % descrip = [deblank(descrip) '   ' hdr{i}.PatientsName];

    if ~true % LEFT-HANDED STORAGE
        mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
        volume = flipud(volume);
    end

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
    if isfield(hdr{i},'RescaleSlope') && hdr{i}.RescaleSlope ~= 1
        RescaleSlope     = hdr{i}.RescaleSlope;
    end
    if isfield(hdr{i},'RescaleIntercept') && hdr{i}.RescaleIntercept ~= 0
        RescaleIntercept = hdr{i}.RescaleIntercept;
    end
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
end
spm_progress_bar('Clear');


%==========================================================================
% function fnames = convert_standard(hdr,root_dir,format,out_dir)
%==========================================================================
function fnames = convert_standard(hdr,root_dir,format,out_dir)
hdr = sort_into_volumes(hdr);
fnames = cell(length(hdr),1);
for i=1:length(hdr)
    fnames{i} = write_volume(hdr{i},root_dir,format,out_dir);
end


%==========================================================================
% function vol = sort_into_volumes(hdr)
%==========================================================================
function vol = sort_into_volumes(hdr)

%
% First of all, sort into volumes based on relevant
% fields in the header.
%

vol{1}{1} = hdr{1};
for i=2:length(hdr)
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
    end
    for j=1:length(vol)
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
        end
        if ~isempty(ice1) && isfield(vol{j}{1},'CSAImageHeaderInfo') && numel(vol{j}{1}.CSAImageHeaderInfo)>=1 && isfield(vol{j}{1}.CSAImageHeaderInfo(1),'name')
            % Replace 'X' in ICE_Dims by '-1'
            ice2 = sscanf( ...
                strrep(get_numaris4_val(vol{j}{1}.CSAImageHeaderInfo,'ICE_Dims'), ...
                'X', '-1'), '%i_%i_%i_%i_%i_%i_%i_%i_%i')';
            if ~isempty(ice2)
                identical_ice_dims=all(ice1(dimsel)==ice2(dimsel));
            else
                identical_ice_dims = 0; % have ice1 but not ice2, ->
                % something must be different
            end
        else
            identical_ice_dims = 1; % No way of knowing if there is no CSAImageHeaderInfo
        end
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
            %end
            % For raw image data, tell apart real/complex or phase/magnitude
            if isfield(hdr{i},'ImageType') && isfield(vol{j}{1}, 'ImageType')
                match = match && strcmp(hdr{i}.ImageType, vol{j}{1}.ImageType);
            end
            if isfield(hdr{i},'SequenceName') && isfield(vol{j}{1}, 'SequenceName')
                match = match && strcmp(hdr{i}.SequenceName,vol{j}{1}.SequenceName);
            end
            if isfield(hdr{i},'SeriesInstanceUID') && isfield(vol{j}{1}, 'SeriesInstanceUID')
                match = match && strcmp(hdr{i}.SeriesInstanceUID,vol{j}{1}.SeriesInstanceUID);
            end
            if isfield(hdr{i},'EchoNumbers')  && isfield(vol{j}{1}, 'EchoNumbers')
                match = match && hdr{i}.EchoNumbers == vol{j}{1}.EchoNumbers;
            end
            if isfield(hdr{i},    'GE_ImageType') && numel(   hdr{i}.GE_ImageType)==1 && ...
               isfield(vol{j}{1}, 'GE_ImageType') && numel(vol{j}{1}.GE_ImageType)==1
                match = match && hdr{i}.GE_ImageType == vol{j}{1}.GE_ImageType;
            end
        catch
            match = 0;
        end
        if match
            vol{j}{end+1} = hdr{i};
            break;
        end
    end
    if ~match
        vol{end+1}{1} = hdr{i};
    end
end

%dcm = vol;
%save('dicom_headers.mat','dcm');

%
% Secondly, sort volumes into ascending/descending
% slices depending on .ImageOrientationPatient field.
%

vol2 = {};
for j=1:length(vol)
    orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
    proj   = null(orient');
    if det([orient proj])<0, proj = -proj; end

    z      = zeros(length(vol{j}),1);
    for i=1:length(vol{j})
        z(i)  = vol{j}{i}.ImagePositionPatient(:)'*proj;
    end
    [z,index] = sort(z);
    vol{j}    = vol{j}(index);
    if length(vol{j})>1
        % dist      = diff(z);
        if any(diff(z)==0)
            tmp = sort_into_vols_again(vol{j});
            vol{j} = tmp{1};
            vol2 = {vol2{:} tmp{2:end}};
        end
    end
end
vol = {vol{:} vol2{:}};
for j=1:length(vol)
    if length(vol{j})>1
        orient = reshape(vol{j}{1}.ImageOrientationPatient,[3 2]);
        proj   = null(orient');
        if det([orient proj])<0, proj = -proj; end
        z      = zeros(length(vol{j}),1);
        for i=1:length(vol{j})
            z(i)  = vol{j}{i}.ImagePositionPatient(:)'*proj;
        end
        dist = diff(sort(z));
        if sum((dist-mean(dist)).^2)/length(dist)>1e-4
            fprintf('***************************************************\n');
            fprintf('* VARIABLE SLICE SPACING                          *\n');
            fprintf('* This may be due to missing DICOM files.         *\n');
            PatientID = 'anon';
            if checkfields(vol{j}{1},'PatientID'), PatientID = deblank(vol{j}{1}.PatientID); end
            if checkfields(vol{j}{1},'SeriesNumber','AcquisitionNumber','InstanceNumber')
                fprintf('*    %s / %d / %d / %d \n',...
                    PatientID, vol{j}{1}.SeriesNumber, ...
                    vol{j}{1}.AcquisitionNumber, vol{j}{1}.InstanceNumber);
                fprintf('*                                                 *\n');
            end
            fprintf('*  %20.4g                           *\n', dist);
            fprintf('***************************************************\n');
        end
    end
end


%==========================================================================
% function vol2 = sort_into_vols_again(volj)
%==========================================================================
function vol2 = sort_into_vols_again(volj)
if ~isfield(volj{1},'InstanceNumber'),
    fprintf('***************************************************\n');
    fprintf('* The slices may be all mixed up and the data     *\n');
    fprintf('* not really usable.  Talk to your physicists     *\n');
    fprintf('* about this.                                     *\n');
    fprintf('***************************************************\n');
    vol2 = {volj};
    return;
end

fprintf('***************************************************\n');
fprintf('* The AcquisitionNumber counter does not appear   *\n');
fprintf('* to be changing from one volume to another.      *\n');
fprintf('* Another possible explanation is that the same   *\n');
fprintf('* DICOM slices are used multiple times.           *\n');
%fprintf('* Talk to your MR sequence developers or scanner  *\n');
%fprintf('* supplier to have this fixed.                    *\n');
fprintf('* The conversion is having to guess how slices    *\n');
fprintf('* should be arranged into volumes.                *\n');
PatientID = 'anon';
if checkfields(volj{1},'PatientID'), PatientID = deblank(volj{1}.PatientID); end
if checkfields(volj{1},'SeriesNumber','AcquisitionNumber')
    fprintf('*    %s / %d / %d\n',...
        PatientID, volj{1}.SeriesNumber, ...
        volj{1}.AcquisitionNumber);
end
fprintf('***************************************************\n');

z      = zeros(length(volj),1);
t      = zeros(length(volj),1);
d      = zeros(length(volj),1);
orient = reshape(volj{1}.ImageOrientationPatient,[3 2]);
proj   = null(orient');
if det([orient proj])<0, proj = -proj; end

for i=1:length(volj),
    z(i)  = volj{i}.ImagePositionPatient(:)'*proj;
    t(i)  = volj{i}.InstanceNumber;
end
% msg = 0;
[t,index] = sort(t);
volj      = volj(index);
z         = z(index);
msk       = find(diff(t)==0);
if any(msk)
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
end

%if any(diff(sort(t))~=1), msg = 1; end;
[z,index] = sort(z);
volj      = volj(index);
t         = t(index);
vol2      = {};
while ~all(d)
    i  = find(~d);
    i  = i(1);
    i  = find(z==z(i));
    [t(i),si] = sort(t(i));
    volj(i)   = volj(i(si));
    for i1=1:length(i)
        if length(vol2)<i1, vol2{i1} = {}; end;
        vol2{i1} = {vol2{i1}{:} volj{i(i1)}};
    end
    d(i) = 1;
end

msg = 0;
len = length(vol2{1});
for i=2:length(vol2)
    if length(vol2{i}) ~= len
        msg = 1;
        break;
    end
end
if msg
    fprintf('***************************************************\n');
    fprintf('* There are missing DICOM files, so the the       *\n');
    fprintf('* resulting volumes may be messed up.             *\n');
    PatientID = 'anon';
    if checkfields(volj{1},'PatientID'), PatientID = deblank(volj{1}.PatientID); end
    if checkfields(volj{1},'SeriesNumber','AcquisitionNumber'),
        fprintf('*    %s / %d / %d\n',...
            PatientID, volj{1}.SeriesNumber, ...
            volj{1}.AcquisitionNumber);
    end;
    fprintf('***************************************************\n');
end


%==========================================================================
% function fname = write_volume(hdr,root_dir,format,out_dir)
%==========================================================================
function fname = write_volume(hdr,root_dir,format,out_dir)

% Output filename
%--------------------------------------------------------------------------
fname = getfilelocation(hdr{1}, root_dir,'s',format,out_dir);

% Image dimensions
%--------------------------------------------------------------------------
nc = hdr{1}.Columns;
nr = hdr{1}.Rows;

if length(hdr) == 1 && isfield(hdr{1},'NumberOfFrames') && hdr{1}.NumberOfFrames > 1
    if isfield(hdr{1},'ImagePositionPatient') &&...
       isfield(hdr{1},'ImageOrientationPatient') &&...
       isfield(hdr{1},'SliceThickness') &&...
       isfield(hdr{1},'StartOfPixelData') &&...
       isfield(hdr{1},'SizeOfPixelData')

       orient           = reshape(hdr{1}.ImageOrientationPatient,[3 2]);
       orient(:,3)      = null(orient');
       if det(orient)<0, orient(:,3) = -orient(:,3); end
       slicevec         = orient(:,3);
   
       hdr_temp = cell(1,hdr{1}.NumberOfFrames); % alternative: NumberofSlices
       hdr_temp{1} = hdr{1};
       hdr_temp{1}.SizeOfPixelData           = hdr{1}.SizeOfPixelData / hdr{1}.NumberOfFrames;
       for sn = 2 : hdr{1}.NumberOfFrames
           hdr_temp{sn}                      = hdr{1};
           hdr_temp{sn}.ImagePositionPatient = hdr{1}.ImagePositionPatient + (sn-1) * hdr{1}.SliceThickness * slicevec;
           hdr_temp{sn}.SizeOfPixelData      = hdr_temp{1}.SizeOfPixelData;
           hdr_temp{sn}.StartOfPixelData     = hdr{1}.StartOfPixelData + (sn-1) * hdr_temp{1}.SizeOfPixelData;
       end
       hdr = hdr_temp;
   else
       error('spm_dicom_convert:write_volume','TAGS missing in DICOM file.');
   end
end

dim    = [nc nr length(hdr)];
dt     = determine_datatype(hdr{1});

% Orientation information
%--------------------------------------------------------------------------
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

if length(hdr)>1
    x2 = [1;1;dim(3); 1];
    y2 = [hdr{end}.ImagePositionPatient(:); 1];
else
    orient           = reshape(hdr{1}.ImageOrientationPatient,[3 2]);
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
    if checkfields(hdr{1},'SliceThickness')
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
%--------------------------------------------------------------------------
if checkfields(hdr{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate')
    if isfield(hdr{1},'ScanOptions')
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
end

if ~true, % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end

% Write the image volume
%--------------------------------------------------------------------------
spm_progress_bar('Init',length(hdr),['Writing ' fname], 'Planes written');
pinfos = [ones(length(hdr),1) zeros(length(hdr),1)];
for i=1:length(hdr)
    if isfield(hdr{i},'RescaleSlope'),     pinfos(i,1) = hdr{i}.RescaleSlope;     end 
    if isfield(hdr{i},'RescaleIntercept'), pinfos(i,2) = hdr{i}.RescaleIntercept; end
end

if any(any(diff(pinfos,1)))
    % Ensure random numbers are reproducible (see later)
    % when intensities are dithered to prevent aliasing effects.
    rand('state',0);
end

volume = zeros(dim);
for i=1:length(hdr)
    plane = read_image_data(hdr{i});

    if any(any(diff(pinfos,1)))
        % This is to prevent aliasing effects in any subsequent histograms
        % of the data (eg for mutual information coregistration).
        % It's a bit inelegant, but probably necessary for when slices are
        % individually rescaled.
        plane = double(plane) + rand(size(plane)) - 0.5;
    end

    if pinfos(i,1)~=1, plane = plane*pinfos(i,1); end
    if pinfos(i,2)~=0, plane = plane+pinfos(i,2); end

    plane = fliplr(plane);
    if ~true, plane = flipud(plane); end % LEFT-HANDED STORAGE
    volume(:,:,i) = plane;
    spm_progress_bar('Set',i);
end

if ~any(any(diff(pinfos,1)))
    % Same slopes and intercepts for all slices
    pinfo = pinfos(1,:);
else
    % Variable slopes and intercept (maybe PET/SPECT)
    mx = max(volume(:));
    mn = min(volume(:));

    %  Slope and Intercept
    %  32767*pinfo(1) + pinfo(2) = mx
    % -32768*pinfo(1) + pinfo(2) = mn
    % pinfo = ([32767 1; -32768 1]\[mx; mn])';

    % Slope only
    dt    = 'int16-be';
    pinfo = [max(mx/32767,-mn/32768) 0];
end

N      = nifti;
N.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = descrip;
create(N);
N.dat(:,:,:) = volume;
spm_progress_bar('Clear');


%==========================================================================
% function fnames = convert_spectroscopy(hdr,root_dir,format,out_dir)
%==========================================================================
function fnames = convert_spectroscopy(hdr,root_dir,format,out_dir)
fnames = cell(length(hdr),1);
for i=1:length(hdr)
    fnames{i} = write_spectroscopy_volume(hdr(i),root_dir,format,out_dir);
end


%==========================================================================
% function fname = write_spectroscopy_volume(hdr,root_dir,format,out_dir)
%==========================================================================
function fname = write_spectroscopy_volume(hdr,root_dir,format,out_dir)
% Output filename
%-------------------------------------------------------------------
fname = getfilelocation(hdr{1}, root_dir,'S',format,out_dir);

% private field to use - depends on SIEMENS software version
if isfield(hdr{1}, 'CSANonImageHeaderInfoVA')
    privdat = hdr{1}.CSANonImageHeaderInfoVA;
elseif isfield(hdr{1}, 'CSANonImageHeaderInfoVB')
    privdat = hdr{1}.CSANonImageHeaderInfoVB;
else
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end

% Image dimensions
%--------------------------------------------------------------------------
nc = get_numaris4_numval(privdat,'Columns');
nr = get_numaris4_numval(privdat,'Rows');
% Guess number of timepoints in file - I don't know for sure whether this should be
% 'DataPointRows'-by-'DataPointColumns', 'SpectroscopyAcquisitionDataColumns'
% or sSpecPara.lVectorSize from SIEMENS ASCII header
% ntp = get_numaris4_numval(privdat,'DataPointRows')*get_numaris4_numval(privdat,'DataPointColumns');
ac = read_ascconv(hdr{1});
try
    ntp = ac.sSpecPara.lVectorSize;
catch
    disp('Don''t know how to handle these spectroscopy data');
    fname = '';
    return;
end
dim    = [nc nr numel(hdr) 2 ntp];
dt     = spm_type('float32'); % Fixed datatype

% Orientation information
%--------------------------------------------------------------------------
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
    %----------------------------------------------------------------------
    % NOTE: Internally the position vector of the CSI matrix shows to the
    % outer border of the first voxel. Therefore the position vector has to
    % be corrected. (Note: The convention of Siemens spectroscopy raw data
    % is in contrast to the DICOM standard where the position vector points
    % to the center of the first voxel.)
    %----------------------------------------------------------------------
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

if length(hdr)>1
    error('spm_dicom_convert:spectroscopy',...
        'Don''t know how to handle multislice spectroscopy data.');
else
    orient(:,3)      = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
    z = get_numaris4_numval(privdat,...
        'VoiThickness');
    if isempty(z)
        z = get_numaris4_numval(privdat,...
            'SliceThickness');
    end
    if isempty(z)
        warning('spm_dicom_convert:spectroscopy',...
            'Can not determine voxel thickness.');
        z = 1;
    end
    x2 = [0;0;1;0];
    y2 = [orient*[0;0;z];0];
end
dicom_to_patient = [y1 y2 R]/[x1 x2 eye(4,2)];
mat              = patient_to_tal*dicom_to_patient*shift_vx*analyze_to_dicom;

% Possibly useful information
%--------------------------------------------------------------------------
if checkfields(hdr{1},'AcquisitionTime','MagneticFieldStrength','MRAcquisitionType',...
        'ScanningSequence','RepetitionTime','EchoTime','FlipAngle',...
        'AcquisitionDate')
    tim = datevec(hdr{1}.AcquisitionTime/(24*60*60));
    descrip = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg %s %d:%d:%.5g',...
        hdr{1}.MagneticFieldStrength, hdr{1}.MRAcquisitionType,...
        deblank(hdr{1}.ScanningSequence),...
        hdr{1}.RepetitionTime,hdr{1}.EchoTime,hdr{1}.FlipAngle,...
        datestr(hdr{1}.AcquisitionDate),tim(4),tim(5),tim(6));
else
    descrip = hdr{1}.Modality;
end

if ~true % LEFT-HANDED STORAGE
    mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
end

% Write the image volume
%--------------------------------------------------------------------------
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
% Store LCMODEL control/raw info in N.extras
N.extras      = struct('MagneticFieldStrength',...
                       get_numaris4_numval(privdat,'MagneticFieldStrength'),...
                       'TransmitterReferenceAmplitude',...
                       get_numaris4_numval(privdat,'TransmitterReferenceAmplitude'),...
                       'ImagingFrequency',...
                       get_numaris4_numval(privdat,'ImagingFrequency'),...
                       'EchoTime',...
                       get_numaris4_numval(privdat,'EchoTime'),...
                       'RealDwellTime',...
                       get_numaris4_numval(privdat,'RealDwellTime'));
create(N);

% Read data, swap dimensions
data = permute(reshape(read_spect_data(hdr{1},ntp),dim([4 5 1 2 3])), ...
                [3 4 5 1 2]);
% plane = fliplr(plane);

N.dat(:,:,:,:,:) = data;


%==========================================================================
% function [images,guff] = select_tomographic_images(hdr)
%==========================================================================
function [images,guff] = select_tomographic_images(hdr)
images = {};
guff   = {};
for i=1:length(hdr)
    if ~checkfields(hdr{i},'Modality') || ...
            ~(strcmp(hdr{i}.Modality,'MR') || ...
              strcmp(hdr{i}.Modality,'PT') || ...
              strcmp(hdr{i}.Modality,'NM') || ...
              strcmp(hdr{i}.Modality,'CT'))
        if checkfields(hdr{i},'Modality')
            fprintf('File "%s" can not be converted because it is of type "%s", which is not MRI, CT, NM or PET.\n', hdr{i}.Filename, hdr{i}.Modality);
        else
            fprintf('File "%s" can not be converted because it does not encode an image.\n', hdr{i}.Filename);
        end
        guff = [guff(:)',hdr(i)];
 
    elseif ~checkfields(hdr{i},'StartOfPixelData','SamplesPerPixel',...
            'Rows','Columns','BitsAllocated','BitsStored','HighBit','PixelRepresentation')
        fprintf('Cant find "Image Pixel" information for "%s".\n',hdr{i}.Filename);
        guff = [guff(:)',hdr(i)];
        
    elseif ~(checkfields(hdr{i},'PixelSpacing','ImagePositionPatient','ImageOrientationPatient') ...
            || isfield(hdr{i},'CSANonImageHeaderInfoVA') || isfield(hdr{i},'CSANonImageHeaderInfoVB'))
        if isfield(hdr{i},'SharedFunctionalGroupsSequence') || isfield(hdr{i},'PerFrameFunctionalGroupsSequence'),
            fprintf(['\n"%s" appears to be multi-frame DICOM.\n'...
                     'Converting these data is still experimental and has only been tested on a very small number of\n'...
                     'multiframe DICOM files. Feedback about problems would be appreciated - particularly if you can\n'...
                     'give us examples of problematic data (providing there are no subject confidentiality issues).\n\n'], hdr{i}.Filename);
            images = [images(:)',hdr(i)];
        else
            fprintf('Can''t find "Image Plane" information for "%s".\n',hdr{i}.Filename);
            guff = [guff(:)',hdr(i)];
        end

    elseif ~checkfields(hdr{i},'SeriesNumber','AcquisitionNumber','InstanceNumber')
       %disp(['Cant find suitable filename info for "' hdr{i}.Filename '".']);
        if ~isfield(hdr{i},'SeriesNumber')
            fprintf('Setting SeriesNumber to 1.\n');
            hdr{i}.SeriesNumber = 1;
            images = [images(:)',hdr(i)];
        end
        if ~isfield(hdr{i},'AcquisitionNumber')
            if isfield(hdr{i},'Manufacturer') && ~isempty(strfind(upper(hdr{1}.Manufacturer), 'PHILIPS'))
                % Philips oddity
                if isfield(hdr{i},'InstanceNumber')
                    hdr{i}.AcquisitionNumber = hdr{i}.InstanceNumber;
                else
                    fprintf('Setting AcquisitionNumber to 1.\n');
                    hdr{i}.AcquisitionNumber = 1;
                end
            else
                fprintf('Setting AcquisitionNumber to 1.\n');
                hdr{i}.AcquisitionNumber = 1;
            end
            images = [images(:)',hdr(i)];
        end
        if ~isfield(hdr{i},'InstanceNumber')
            fprintf('Setting InstanceNumber to 1.\n');
            hdr{i}.InstanceNumber = 1;
            images = [images(:)',hdr(i)];
        end
    %elseif isfield(hdr{i},'Private_2001_105f'),
    %    % This field corresponds to: > Stack Sequence 2001,105F SQ VNAP, COPY
    %    % http://www.medical.philips.com/main/company/connectivity/mri/index.html
    %    % No documentation about this private field is yet available.
    %    disp('Cant yet convert Philips Intera DICOM.');
    %    guff = {guff{:},hdr{i}};
    else
        images = [images(:)',hdr(i)];
    end
end


%==========================================================================
% function [multiframe,other] = select_multiframe(hdr)
%==========================================================================
function [multiframe,other] = select_multiframe(hdr)
multiframe = {};
other      = {};
for i=1:length(hdr)
    if isfield(hdr{i},'SharedFunctionalGroupsSequence') || isfield(hdr{i},'PerFrameFunctionalGroupsSequence'),
        multiframe = [multiframe(:)',hdr(i)];
    else
        other      = [other(:)',hdr(i)];
    end
end


%==========================================================================
% function [mosaic,standard] = select_mosaic_images(hdr)
%==========================================================================
function [mosaic,standard] = select_mosaic_images(hdr)
mosaic   = {};
standard = {};
for i=1:length(hdr)
    if ~checkfields(hdr{i},'ImageType','CSAImageHeaderInfo') ||...
            isfield(hdr{i}.CSAImageHeaderInfo,'junk') ||...
            isempty(read_AcquisitionMatrixText(hdr{i})) ||...
            isempty(read_NumberOfImagesInMosaic(hdr{i})) ||...
            read_NumberOfImagesInMosaic(hdr{i}) == 0
        % NumberOfImagesInMosaic seems to be set to zero for pseudo images
        % containing e.g. online-fMRI design matrices, don't treat them as
        % mosaics
        standard = [standard, hdr(i)];
    else
        mosaic   = [mosaic,   hdr(i)];
    end
end


%==========================================================================
% function [spect,images] = select_spectroscopy_images(hdr)
%==========================================================================
function [spect,images] = select_spectroscopy_images(hdr)
spectsel = false(1,numel(hdr));
for i=1:numel(hdr)
    if isfield(hdr{i},'SOPClassUID')
        spectsel(i) = strcmp(hdr{i}.SOPClassUID,'1.3.12.2.1107.5.9.1');
    end
end
spect  = hdr(spectsel);
images = hdr(~spectsel);


%==========================================================================
% function [standard, guff] = select_last_guff(standard, guff)
%==========================================================================
function [standard, guff] = select_last_guff(standard, guff)
% See https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=spm;5b69d495.1108
i = find(cellfun(@(x) ~isfield(x,'ImageOrientationPatient'),standard));
guff = [guff, standard(i)];
standard(i) = [];


%==========================================================================
% function ok = checkfields(hdr,varargin)
%==========================================================================
function ok = checkfields(hdr,varargin)
ok = 1;
for i=1:(nargin-1)
    if ~isfield(hdr,varargin{i})
        ok = 0;
        break;
    end
end


%==========================================================================
% function clean = strip_unwanted(dirty)
%==========================================================================
function clean = strip_unwanted(dirty)
msk = (dirty>='a'&dirty<='z') | (dirty>='A'&dirty<='Z') |...
      (dirty>='0'&dirty<='9') | dirty=='_';
clean = dirty(msk);


%==========================================================================
% function img = read_image_data(hdr)
%==========================================================================
function img = read_image_data(hdr)
img = [];

if hdr.SamplesPerPixel ~= 1,
    warning('spm:dicom','%s: SamplesPerPixel = %d - cant be an MRI.', hdr.Filename, hdr.SamplesPerPixel);
    return;
end;

prec = ['ubit' num2str(hdr.BitsAllocated) '=>' 'uint32'];

if isfield(hdr,'TransferSyntaxUID') && strcmp(hdr.TransferSyntaxUID,'1.2.840.10008.1.2.2') && strcmp(hdr.VROfPixelData,'OW'),
    fp = fopen(hdr.Filename,'r','ieee-be');
else
    fp = fopen(hdr.Filename,'r','ieee-le');
end;
if fp==-1,
    warning('spm:dicom','%s: Cant open file.', hdr.Filename);
    return;
end;

if isfield(hdr,'NumberOfFrames'),
    NFrames = hdr.NumberOfFrames;
else
    NFrames = 1;
end

if isfield(hdr,'TransferSyntaxUID')
    switch(hdr.TransferSyntaxUID)
    case {'1.2.840.10008.1.2.4.50','1.2.840.10008.1.2.4.51',... % 8 bit JPEG & 12 bit JPEG
          '1.2.840.10008.1.2.4.57','1.2.840.10008.1.2.4.70',... % lossless NH JPEG & lossless NH, 1st order
          '1.2.840.10008.1.2.4.80','1.2.840.10008.1.2.4.81',... % lossless JPEG-LS & near lossless JPEG-LS
          '1.2.840.10008.1.2.4.90','1.2.840.10008.1.2.4.91',... % lossless JPEG 2000 & possibly lossy JPEG 2000, Part 1
          '1.2.840.10008.1.2.4.92','1.2.840.10008.1.2.4.93' ... % lossless JPEG 2000 & possibly lossy JPEG 2000, Part 2
         },
        % try to read PixelData as JPEG image
        fseek(fp,hdr.StartOfPixelData,'bof');
        fread(fp,2,'uint16'); % uint16 encoding 65534/57344 (Item)
        offset = double(fread(fp,1,'uint32')); % followed by 4 0 0 0
        fread(fp,2,'uint16'); % uint16 encoding 65534/57344 (Item)
        fread(fp,offset);
        
        sz  = double(fread(fp,1,'*uint32'));
        img = fread(fp,sz,'*uint8');

        % Next uint16 seem to encode 65534/57565 (SequenceDelimitationItem), followed by 0 0

        % save PixelData into temp file - imread and its subroutines can only
        % read from file, not from memory
        tfile = tempname;
        tfp   = fopen(tfile,'w+');
        fwrite(tfp,img,'uint8');
        fclose(tfp);

        % read decompressed data, transpose to match DICOM row/column order
        img = uint32(imread(tfile)');
        delete(tfile);
    case {'1.2.840.10008.1.2.4.94' ,'1.2.840.10008.1.2.4.95' ,... % JPIP References & JPIP Referenced Deflate Transfer
          '1.2.840.10008.1.2.4.100','1.2.840.10008.1.2.4.101',... % MPEG2 MP@ML & MPEG2 MP@HL
          '1.2.840.10008.1.2.4.102',                          ... % MPEG-4 AVC/H.264 High Profile and BD-compatible
         }
         warning('spm:dicom',[hdr.Filename ': cant deal with JPIP/MPEG data (' hdr.TransferSyntaxUID ')']);
    otherwise
        fseek(fp,hdr.StartOfPixelData,'bof');
        img = fread(fp,hdr.Rows*hdr.Columns*NFrames,prec);
    end
else
    fseek(fp,hdr.StartOfPixelData,'bof');
    img = fread(fp,hdr.Rows*hdr.Columns*NFrames,prec);
end
fclose(fp);
if numel(img)~=hdr.Rows*hdr.Columns*NFrames,
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

img = reshape(img,[hdr.Columns,hdr.Rows,NFrames]);


%==========================================================================
% function img = read_spect_data(hdr,privdat)
%==========================================================================
function img = read_spect_data(hdr,ntp)
% Data is stored as complex float32 values, timepoint by timepoint, voxel
% by voxel. Reshaping is done in write_spectroscopy_volume.
if ntp*2*4 ~= hdr.SizeOfCSAData
    warning('spm:dicom', [hdr.Filename,': Data size mismatch.']);
end
fp = fopen(hdr.Filename,'r','ieee-le');
fseek(fp,hdr.StartOfCSAData,'bof');
img = fread(fp,2*ntp,'float32');
fclose(fp);


%==========================================================================
% function nrm = read_SliceNormalVector(hdr)
%==========================================================================
function nrm = read_SliceNormalVector(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'SliceNormalVector');
for i=1:3,
    nrm(i,1) = sscanf(val(i,:),'%g');
end


%==========================================================================
% function n = read_NumberOfImagesInMosaic(hdr)
%==========================================================================
function n = read_NumberOfImagesInMosaic(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'NumberOfImagesInMosaic');
n   = sscanf(val','%d');
if isempty(n), n=[]; end


%==========================================================================
% function dim = read_AcquisitionMatrixText(hdr)
%==========================================================================
function dim = read_AcquisitionMatrixText(hdr)
str = hdr.CSAImageHeaderInfo;
val = get_numaris4_val(str,'AcquisitionMatrixText');
dim = sscanf(val','%d*%d')';
if length(dim)==1
    dim = sscanf(val','%dp*%d')';
end
if isempty(dim), dim=[]; end


%==========================================================================
% function val = get_numaris4_val(str,name)
%==========================================================================
function val = get_numaris4_val(str,name)
name = deblank(name);
val  = {};
for i=1:length(str)
    if strcmp(deblank(str(i).name),name)
        for j=1:str(i).nitems
            if  str(i).item(j).xx(1)
                val = [val {str(i).item(j).val}];
            end
        end
        break;
    end
end
val = strvcat(val{:});


%==========================================================================
% function val = get_numaris4_numval(str,name)
%==========================================================================
function val = get_numaris4_numval(str,name)
val1 = get_numaris4_val(str,name);
val  = zeros(size(val1,1),1);
for k = 1:size(val1,1)
    val(k)=str2num(val1(k,:));
end


%==========================================================================
% function fname = getfilelocation(hdr,root_dir,prefix,format,out_dir)
%==========================================================================
function fname = getfilelocation(hdr,root_dir,prefix,format,out_dir)

if nargin < 3
    prefix = 'f';
end

if strncmp(root_dir,'ice',3)
    root_dir = root_dir(4:end);
    imtype = textscan(hdr.ImageType,'%s','delimiter','\\');
    try
        imtype = imtype{1}{3};
    catch
        imtype = '';
    end
    prefix = [prefix imtype get_numaris4_val(hdr.CSAImageHeaderInfo,'ICE_Dims')];
end

if isfield(hdr,'PatientID'),         PatientID         = deblank(hdr.PatientID); else PatientID         = 'anon'; end
if isfield(hdr,'EchoNumbers'),       EchoNumbers       = hdr.EchoNumbers;        else EchoNumbers       = 0;      end
if isfield(hdr,'SeriesNumber'),      SeriesNumber      = hdr.SeriesNumber;       else SeriesNumber      = 0;      end
if isfield(hdr,'AcquisitionNumber'), AcquisitionNumber = hdr.AcquisitionNumber;  else AcquisitionNumber = 0;      end
if isfield(hdr,'InstanceNumber'),    InstanceNumber    = hdr.InstanceNumber;     else InstanceNumber    = 0;      end

ImTyp = '';
if isfield(hdr,'GE_ImageType')
    if numel(hdr.GE_ImageType)==1
        switch hdr.GE_ImageType
        case 1
            ImTyp = '-Phase';
        case 2
            ImTyp = '-Real';
        case 3
            ImTyp = '-Imag';
        end
    end
end

if strcmp(root_dir, 'flat')
    % Standard SPM file conversion
    %----------------------------------------------------------------------
    if checkfields(hdr,'SeriesNumber','AcquisitionNumber')
        if checkfields(hdr,'EchoNumbers')
            fname = sprintf('%s%s-%.4d-%.5d-%.6d-%.2d%s.%s', prefix, strip_unwanted(PatientID),...
                SeriesNumber, AcquisitionNumber, InstanceNumber, EchoNumbers, ImTyp, format);
        else
            fname = sprintf('%s%s-%.4d-%.5d-%.6d%s.%s', prefix, strip_unwanted(PatientID),...
                SeriesNumber, AcquisitionNumber, InstanceNumber, ImTyp, format);
        end
    else
        fname = sprintf('%s%s-%.6d%s.%s',prefix, ...
            strip_unwanted(PatientID), InstanceNumber, ImTyp, format);
    end

    fname = fullfile(out_dir,fname);
    return;
end

% more fancy stuff - sort images into subdirectories
if isfield(hdr,'StudyTime')
    m = sprintf('%02d', floor(rem(hdr.StudyTime/60,60)));
    h = sprintf('%02d', floor(hdr.StudyTime/3600));
else
    m = '00';
    h = '00';
end
if isfield(hdr,'AcquisitionTime'),   AcquisitionTime   = hdr.AcquisitionTime;            else AcquisitionTime   = 100;       end;
if isfield(hdr,'StudyDate'),         StudyDate         = hdr.StudyDate;                  else StudyDate         = 100;       end; % Obscure Easter Egg
if isfield(hdr,'PatientsName'),      PatientsName      = deblank(hdr.PatientsName);      else PatientsName      = 'anon';    end
if isfield(hdr,'SeriesDescription'), SeriesDescription = deblank(hdr.SeriesDescription); else SeriesDescription = 'unknown'; end
if isfield(hdr,'ProtocolName'),
    ProtocolName = deblank(hdr.ProtocolName);
else
    if isfield(hdr,'SequenceName')
        ProtocolName = deblank(hdr.SequenceName);
    else
        ProtocolName='unknown';
    end
end

studydate = sprintf('%s_%s-%s', datestr(StudyDate,'yyyy-mm-dd'), h,m);
switch root_dir
    case {'date_time','series'}
        id = studydate;
    case {'patid', 'patid_date', 'patname'},
        id = strip_unwanted(PatientID);
end
serdes   = strrep(strip_unwanted(SeriesDescription), strip_unwanted(ProtocolName),'');
protname = sprintf('%s%s_%.4d',strip_unwanted(ProtocolName), serdes, SeriesNumber);
switch root_dir
    case 'date_time',
        dname = fullfile(out_dir, id, protname);
    case 'patid',
        dname = fullfile(out_dir, id, protname);
    case 'patid_date',
        dname = fullfile(out_dir, id, studydate, protname);
    case 'patname',
        dname = fullfile(out_dir, strip_unwanted(PatientsName), id, protname);
    case 'series',
        dname = fullfile(out_dir, protname);
    otherwise
        error('unknown file root specification');
end
if ~exist(dname,'dir')
    mkdir_rec(dname);
end

% some non-product sequences on SIEMENS scanners seem to have problems
% with image numbering in MOSAICs - doublettes, unreliable ordering
% etc. To distinguish, always include Acquisition time in image name
sa    = sprintf('%02d', floor(rem(AcquisitionTime,60)));
ma    = sprintf('%02d', floor(rem(AcquisitionTime/60,60)));
ha    = sprintf('%02d', floor(AcquisitionTime/3600));
fname = sprintf('%s%s-%s%s%s-%.5d-%.5d-%d%s.%s', prefix, id, ha, ma, sa, ...
        AcquisitionNumber, InstanceNumber, EchoNumbers, ImTyp, format);
fname = fullfile(dname, fname);


%==========================================================================
% function suc = mkdir_rec(str)
%==========================================================================
function suc = mkdir_rec(str)
% works on full pathnames only
if str(end) ~= filesep, str = [str filesep]; end
pos = strfind(str,filesep);
suc = zeros(1,length(pos));
for g=2:length(pos)
    if ~exist(str(1:pos(g)-1),'dir')
        suc(g) = mkdir(str(1:pos(g-1)-1),str(pos(g-1)+1:pos(g)-1));
    end
end


%==========================================================================
% function ret = read_ascconv(hdr)
%==========================================================================
function ret = read_ascconv(hdr)
% In SIEMENS data, there is an ASCII text section with
% additional information items. This section starts with a code
% ### ASCCONV BEGIN <some version string> ###
% and ends with
% ### ASCCONV END ###
% It is read by spm_dicom_headers into an entry 'MrProtocol' in
% CSASeriesHeaderInfo or into an entry 'MrPhoenixProtocol' in
% CSAMiscProtocolHeaderInfoVA or CSAMiscProtocolHeaderVB.
% The additional items are assignments in C syntax, here they are just
% translated according to
% [] -> ()
% "  -> '
% 0xX -> hex2dec('X')
% and collected in a struct.
% In addition, there seems to be "__attribute__" meta information for some
% items. All "field names" starting with "_" are currently silently ignored.
ret=struct;

% get ascconv data
if isfield(hdr, 'CSAMiscProtocolHeaderInfoVA')
    X = get_numaris4_val(hdr.CSAMiscProtocolHeaderInfoVA,'MrProtocol');
elseif isfield(hdr, 'CSAMiscProtocolHeaderInfoVB')
    X = get_numaris4_val(hdr.CSAMiscProtocolHeaderInfoVB,'MrPhoenixProtocol');
elseif isfield(hdr, 'CSASeriesHeaderInfo')
    X = get_numaris4_val(hdr.CSASeriesHeaderInfo,'MrProtocol');
else
    return;
end

X = regexprep(X,'^.*### ASCCONV BEGIN [^#]*###(.*)### ASCCONV END ###.*$','$1');

if ~isempty(X)
    tokens = textscan(char(X),'%s', ...
        'delimiter',char(10));
    tokens{1}=regexprep(tokens{1},{'\[([0-9]*)\]','"(.*)"','^([^"]*)0x([0-9a-fA-F]*)','#.*','^.*\._.*$'},{'($1+1)','''$1''','$1hex2dec(''$2'')','',''});
    % If everything would evaluate correctly, we could use
    % eval(sprintf('ret.%s;\n',tokens{1}{:}));
    for k = 1:numel(tokens{1})
        if ~isempty(tokens{1}{k})
            try
                [tlhrh, un] = regexp(tokens{1}{k}, '(?:=)+', 'split', 'match');
                [tlh, un]   = regexp(tlhrh{1}, '(?:\.)+', 'split', 'match');
                tlh = cellfun(@genvarname, tlh, 'UniformOutput',false);
                tlh = sprintf('.%s', tlh{:});
                eval(sprintf('ret%s = %s;', tlh, tlhrh{2}));
            catch
                disp(['AscConv: Error evaluating ''ret.' tokens{1}{k} ''';']);
            end
        end
    end
end


%==========================================================================
% function dt = determine_datatype(hdr)
%==========================================================================
function dt = determine_datatype(hdr)
% Determine what datatype to use for NIfTI images
be = spm_platform('bigend');
if hdr.BitsStored>16
    if hdr.PixelRepresentation
        dt  = [spm_type( 'int32') be];
    else
        dt  = [spm_type('uint32') be];
    end
else
    if hdr.PixelRepresentation 
        dt  = [spm_type( 'int16') be];
    else
        dt  = [spm_type('uint16') be];
    end
end

%==========================================================================
% function fspe = convert_multiframes(hdr,root_dir,format,out_dir)
%==========================================================================
function fspe = convert_multiframes(hdr,root_dir,format,out_dir)
fspe = {};
dict = load('spm_dicom_dict.mat');
for i=1:numel(hdr)
    out  = convert_multiframe(hdr{i}, dict, root_dir, format,out_dir);
    fspe = [fspe(:); out(:)];
end

%==========================================================================
% function out = convert_multiframe(H, dict, root_dir, format,out_dir)
%==========================================================================
function out = convert_multiframe(H, dict, root_dir, format,out_dir)
out      = {};
diminfo  = read_DimOrg(H,dict);
dat      = read_FGS(H,diminfo);
N        = numel(dat);
fname0   = getfilelocation(H, root_dir,'MF',format,out_dir);
[pth,nam,ext0] = fileparts(fname0);
fname0   = fullfile(pth,nam);

if ~isfield(dat,'ImageOrientationPatient') || ~isfield(dat,'PixelSpacing') || ~isfield(dat,'ImagePositionPatient'),
    fprintf('"%s" does not seem to have positional information.\n', H.Filename);
    return;
end

%if isfield(dat,'DimensionIndexValues')
%    disp(cat(1,dat.DimensionIndexValues))
%end

% Read the atual image data
volume = read_image_data(H);

% Sort into unique files according to image orientations etc
stuff = [cat(2,dat.ImageOrientationPatient); cat(2,dat.PixelSpacing)]';
if isfield(dat,'StackID')
    stuff = [stuff double(cat(1,dat.StackID))];
end
ds    = find(any(diff(stuff,1,1)~=0,2));
ord   = sparse([],[],[],N,numel(ds)+1);
start = 1;
for i=1:numel(ds),
    ord(start:ds(i),i) = 1;
    start = ds(i)+1;
end
ord(start:size(stuff,1),numel(ds)+1) = 1;
for n=1:size(ord,2),
    % Identify all slices that should go into the same output file
    %ind  = find(all(bsxfun(@eq,stuff,u(n,:)),2));
    ind = find(ord(:,n));
    this = dat(ind); 
    if size(ord,2)>1,
        fname = sprintf('%s_%-.3d',fname0,n);
    else
        fname = fname0;
    end


    % Orientation information
    %--------------------------------------------------------------------------
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


    if any(sum(diff(cat(2,this.ImageOrientationPatient),1,2).^2,1)>0.001),
        fprintf('"%s" contains irregularly oriented slices.\n', H.Filename);
        % break % Option to skip writing out the NIfTI image
    end
    ImageOrientationPatient = this(1).ImageOrientationPatient(:);

    if any(sum(diff(cat(2,this.PixelSpacing),1,2).^2,1)>0.001),
        fprintf('"%s" contains slices with irregularly spaced pixels.\n', H.Filename);
        % break % Option to skip writing out the NIfTI image
    end
    PixelSpacing            = this(1).PixelSpacing(:);

    R   = [reshape(ImageOrientationPatient,3,2)*diag(PixelSpacing); 0 0];

    % Determine the order of the slices by sorting their positions according to where they project
    % on a vector orthogonal to the iamge plane
    orient      = reshape(this(1).ImageOrientationPatient,[3 2]);
    orient(:,3) = null(orient');
    if det(orient)<0, orient(:,3) = -orient(:,3); end
        
    Positions  = cat(2,this.ImagePositionPatient)';
    [proj,slice_order,inv_slice_order] = unique(round(Positions*orient(:,3)*100)/100);
    if any(abs(diff(diff(proj),1,1))>0.025),
        problem1 = true;
    else
        problem1 = false;
    end
 
    inv_time_order = ones(numel(this),1);
    for i=1:numel(slice_order),
        ind = find(inv_slice_order == i);
        if numel(ind)>1,
            if isfield(this,'TemporalPositionIndex'),
                sort_on = cat(2,this(ind).TemporalPositionIndex)';
            else
                sort_on = ind;
            end
            [unused,tsort]      = sort(sort_on);
            inv_time_order(ind) = tsort';
        end
    end

    
    % Image dimensions
    %--------------------------------------------------------------------------
    nc   = H.Columns;
    nr   = H.Rows;
    dim  = [nc nr 1 1];
    dim(3) = max(inv_slice_order);
    dim(4) = max(inv_time_order);

    problem2 = false;
    for i=1:max(inv_time_order)
        if sum(inv_time_order==i)~=dim(3)
            problem2 = true;
        end
    end
    if problem1 || problem2
        fname = [fname '-problem'];
    end
    if problem1
        fprintf('"%s" contains irregularly spaced slices with the following spacings:\n', H.Filename);
        spaces = unique(round(diff(proj)*100)/100);
        fprintf(' %g', spaces);
        fprintf('\nSee %s%s\n\n',fname,ext0);
        % break % Option to skip writing out the NIfTI image
    end
    if problem2
        fprintf('"%s" contains missing slices.\nSee %s%s\n\n', H.Filename,fname,ext0);
        % break % Option to skip writing out the NIfTI image
    end
 
    if dim(3)>1
        y1 = [this(slice_order(  1)).ImagePositionPatient(:); 1];
        y2 = [this(slice_order(end)).ImagePositionPatient(:); 1];
        x2 = [1;1;dim(3); 1];
    else
        orient      = reshape(ImageOrientationPatient,[3 2]);
        orient(:,3) = null(orient');
        if det(orient)<0, orient(:,3) = -orient(:,3); end
        if isfield(H,'SliceThickness'), z = H.SliceThickness; else z = 1; end
        y1 = [this(1).ImagePositionPatient(:); 1];
        y2 = [orient*[0;0;z];0];
        x2 = [0;0;1;0];
    end

    x1               = [1;1;1;1];
    dicom_to_patient = [y1 y2 R]/[x1 x2 eye(4,2)];
    analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)+1) 0]'; 0 0 0 1]; % Flip voxels in y
    patient_to_tal   = diag([-1 -1 1 1]); % Flip mm coords in x and y directions
    mat              = patient_to_tal*dicom_to_patient*analyze_to_dicom;
    flip_lr          = det(mat(1:3,1:3))>0;

    if flip_lr,
        mat    = mat*[-1 0 0 (dim(1)+1); 0 1 0 0; 0 0 1 0; 0 0 0 1];
    end

    % Possibly useful information for descrip field
    %--------------------------------------------------------------------------
    if isfield(H,'AcquisitionTime')
        tim = datevec(H.AcquisitionTime/(24*60*60));
    elseif isfield(H,'StudyTime')
        tim = datevec(H.StudyTime/(24*60*60));
    elseif isfield(H,'ContentTime') 
        tim = datevec(H.ContentTime/(24*60*60));      
    else
        tim = '';
    end
    if ~isempty(tim), tim = sprintf(' %d:%d:%.5g', tim(4),tim(5),tim(6)); end

    if isfield(H,'AcquisitionDate') 
        day = datestr(H.AcquisitionDate);
    elseif isfield(H,'StudyDate')
        day = datestr(H.StudyDate);
    elseif isfield(H,'ContentDate')
        day = datestr(H.ContentDate);
    else
        day = '';
    end
    when = [day tim]; 

    if checkfields(H,'MagneticFieldStrength','MRAcquisitionType',...
                     'ScanningSequence','RepetitionTime','EchoTime','FlipAngle')
        if isfield(H,'ScanOptions')
            ScanOptions = H.ScanOptions;
        else
            ScanOptions = 'no';
        end
        modality = sprintf('%gT %s %s TR=%gms/TE=%gms/FA=%gdeg/SO=%s',...
            H.MagneticFieldStrength, H.MRAcquisitionType,...
            deblank(H.ScanningSequence),...
            H.RepetitionTime,H.EchoTime,H.FlipAngle,...
            ScanOptions);
    else
         modality = [H.Modality ' ' H.ImageType];
    end
    descrip = [modality ' ' when];



    % Sort out datatype, as well as any scalefactors or intercepts
    %--------------------------------------------------------------------------
    pinfos = [ones(length(this),1) zeros(length(this),1)];
    for i=1:length(this)
        if isfield(this(i),'RescaleSlope'),     pinfos(i,1) = this(i).RescaleSlope;     end
        if isfield(this(i),'RescaleIntercept'), pinfos(i,2) = this(i).RescaleIntercept; end
    end
    if ~any(any(diff(pinfos,1)))
        % Same slopes and intercepts for all slices
        dt    = determine_datatype(H);
        pinfo = pinfos(1,:);
    else
        % Variable slopes and intercept (maybe PET/SPECT)
        mx = max(volume(:));
        mn = min(volume(:));

        %  Slope and Intercept
        %  32767*pinfo(1) + pinfo(2) = mx
        % -32768*pinfo(1) + pinfo(2) = mn
        % pinfo = ([32767 1; -32768 1]\[mx; mn])';

        % Slope only
        dt    = 'int16-be';
        pinfo = [max(mx/32767,-mn/32768) 0];

        % Ensure random numbers are reproducible (see later)
        % when intensities are dithered to prevent aliasing effects.
        rand('state',0);
    end


    % Create the header
    %--------------------------------------------------------------------------
    Nii      = nifti;
    fname    = sprintf('%s%s',fname,ext0);
    Nii.dat  = file_array(fname,dim,dt,0,pinfo(1),pinfo(2));
    Nii.mat  = mat;
    Nii.mat0 = mat;
    Nii.mat_intent  = 'Scanner';
    Nii.mat0_intent = 'Scanner';
    Nii.descrip     = descrip;
    create(Nii);
    Nii.dat(end,end,end,end,end) = 0;

    % Write the image volume
    %--------------------------------------------------------------------------
    spm_progress_bar('Init',length(this),['Writing ' fname], 'Planes written');

    for i=1:length(this)

        plane = volume(:,:,this(i).number);

        if any(any(diff(pinfos,1)))
            % This is to prevent aliasing effects in any subsequent histograms
            % of the data (eg for mutual information coregistration).
            % It's a bit inelegant, but probably necessary for when slices are
            % individually rescaled.
            plane = double(plane) + rand(size(plane)) - 0.5;
        end

        if pinfos(i,1)~=1, plane = plane*pinfos(i,1); end
        if pinfos(i,2)~=0, plane = plane+pinfos(i,2); end

        plane = fliplr(plane);
        if flip_lr, plane = flipud(plane); end

        z = inv_slice_order(i);
        t = inv_time_order(i);
        Nii.dat(:,:,z,t) = plane;
        spm_progress_bar('Set',i);
    end

    out = [out; {Nii.dat.fname}];
    spm_progress_bar('Clear');
end




%==========================================================================
% function dim = read_DimOrg(H,dict)
%==========================================================================
function dim = read_DimOrg(H,dict)
% This is only sometimes used.
dim = struct('DimensionIndexPointer',{},'FunctionalGroupPointer',{});

if isfield(H,'DimensionIndexSequence'),
    for i=1:numel(H.DimensionIndexSequence),
        DIP = H.DimensionIndexSequence{i}.DimensionIndexPointer;
        ind = find(dict.group==DIP(1) & dict.element == DIP(2));
        if numel(ind)==1,
            dim(i).DimensionIndexPointer = dict.values(ind).name;
        else
            if rem(DIP(1),2)
                dim(i).DimensionIndexPointer = sprintf('Private_%.4x_%.4x',DIP);
            else
                dim(i).DimensionIndexPointer = sprintf('Tag_%.4x_%.4x',DIP);
            end
        end

        FGP = H.DimensionIndexSequence{i}.FunctionalGroupPointer;
        ind = find(dict.group==FGP(1) & dict.element == FGP(2));
        if numel(ind)==1,
            dim(i).FunctionalGroupPointer     = dict.values(ind).name;
        else
            if rem(DIP(1),2)
                dim(i).FunctionalGroupPointer = sprintf('Private_%.4x_%.4x',DIP);
            else
                dim(i).FunctionalGroupPointer = sprintf('Tag_%.4x_%.4x',DIP);
            end
        end
    end
end
return

%==========================================================================
% dat = read_FGS(H,dim)
%==========================================================================
function dat = read_FGS(H,dim)
dat = struct;
if isfield(H,'PerFrameFunctionalGroupsSequence'),
    % For documentation, see C.7.6.16 "Multi-frame Functional Groups Module"
    % of DICOM Standard PS 3.3 - 2003, Page 261.

    % Multiframe DICOM
    N = numel(H.PerFrameFunctionalGroupsSequence);
    if isfield(H,'NumberOfFrames') && H.NumberOfFrames ~= N,
        fprintf('"%s" has incompatible numbers of frames.', dat.Filename);
    end

    macros = {'PixelMeasuresSequence',{'PixelSpacing','SliceThickness'}
              'FrameContentSequence',{'Frame Acquisition Number','FrameReferenceDatetime',...
                                      'FrameAcquisitionDatetime','FrameAcquisitionDuration',...
                                      'CardiacCyclePosition','RespiratoryCyclePosition',...
                                      'DimensionIndexValues','TemporalPositionIndex','Stack ID',...
                                      'InStackPositionNumber','FrameComments'}
              'PlanePositionSequence',{'ImagePositionPatient'}
              'PlaneOrientationSequence',{'ImageOrientationPatient'}
              'ReferencedImageSequence',{'ReferencedSOPClassUID','ReferencedSOPInstanceUID',...
                                         'ReferencedFrameNumber','PurposeOfReferenceCode'}
              'PixelValueTransformationSequence',{'RescaleIntercept','RescaleSlope','RescaleType'}};

    if isfield(H,'SharedFunctionalGroupsSequence')
        for d=1:numel(dim),
            if isfield(H.SharedFunctionalGroupsSequence{1},dim(d).FunctionalGroupPointer) && ...
               isfield(H.SharedFunctionalGroupsSequence{1}.(dim(d).FunctionalGroupPointer){1},dim(d).DimensionIndexPointer),
                for n=N:-1:1
                    dat(n).(dim(d).DimensionIndexPointer) = H.SharedFunctionalGroupsSequence{1}.(dim(d).FunctionalGroupPointer){1}.(dim(d).DimensionIndexPointer);
                end
            end
        end

        for k1=1:size(macros,1)
            if isfield(H.SharedFunctionalGroupsSequence{1},macros{k1,1})
                for k2=1:numel(macros{k1,2})
                    if (numel(H.SharedFunctionalGroupsSequence{1}.(macros{k1,1}))>0) &&...
                            isfield(H.SharedFunctionalGroupsSequence{1}.(macros{k1,1}){1},macros{k1,2}{k2})
                        for n=N:-1:1
                            dat(n).(macros{k1,2}{k2}) = H.SharedFunctionalGroupsSequence{1}.(macros{k1,1}){1}.(macros{k1,2}{k2});
                        end
                    end
                end
            end
        end
    end

    for n=1:N,
        dat(n).number = n;
        F = H.PerFrameFunctionalGroupsSequence{n};
        for d=1:numel(dim),
            if isfield(F,dim(d).FunctionalGroupPointer) && ...
               isfield(F.(dim(d).FunctionalGroupPointer){1},dim(d).DimensionIndexPointer)
                dat(n).(dim(d).DimensionIndexPointer) = F.(dim(d).FunctionalGroupPointer){1}.(dim(d).DimensionIndexPointer);
            end
        end

        for k1=1:size(macros,1)
            if isfield(F,macros{k1,1}) && (numel(F.(macros{k1,1}))>0)
                for k2=1:numel(macros{k1,2})
                    if isfield(F.(macros{k1,1}){1},macros{k1,2}{k2})
                        dat(n).(macros{k1,2}{k2}) = F.(macros{k1,1}){1}.(macros{k1,2}{k2});
                    end
                end
            end
        end
    end
else
    error('"%s" is not multiframe.', H.FileName);
end


