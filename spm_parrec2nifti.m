function N = spm_parrec2nifti(parfile,opts)
% Import PAR/REC images from Philips scanners into NIfTI
% FORMAT N = spm_parrec2nifti(parfile,opts)
% parfile   - name of PAR file
% opts      - options structure
%    .ext     - NIfTI file extension {'img','nii'} [default: spm_file_ext]
%    .outdir  - output directory [default: pwd]
%
% N         - NIfTI object
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_parrec2nifti.m 6703 2016-01-28 17:36:10Z guillaume $


%-Display warning
%--------------------------------------------------------------------------
fprintf('*************************************************************\n');
fprintf('**         PRELIMINARY SUPPORT OF PHILIPS PAR/REC          **\n');
fprintf('**                  USE AT YOUR OWN RISK                   **\n');
fprintf('*************************************************************\n');

%-Get options
%--------------------------------------------------------------------------
if nargin < 2, opts = struct; end

if ~isfield(opts,'ext'), opts.ext = spm_file_ext; end
if opts.ext(1) ~= '.', opts.ext = ['.' opts.ext]; end

if ~isfield(opts,'outdir') || isempty(opts.outdir)
    opts.outdir = pwd;
end

%-Read PAR header file
%--------------------------------------------------------------------------
hdr = spm_par_hdr(parfile);

%-Find REC data file
%--------------------------------------------------------------------------
if strcmp(spm_file(parfile,'ext'),'PAR')
    recfile = spm_file(parfile,'ext','REC');
else
    recfile = spm_file(parfile,'ext','rec');
end
d = dir(recfile);
if isempty(d), error('Cannot find REC file "%s".',recfile); end
bytes = d.bytes;

%-Write NIfTI header
%--------------------------------------------------------------------------
dim   = [hdr.ImageInfo(1).recon_resolution hdr.MaxNumberOfSlices hdr.MaxNumberOfDynamics];
dim   = double(dim);
dtype = double(hdr.ImageInfo(1).image_pixel_size);
if prod(dim)*dtype/8 ~= bytes
    % can happen for corrupted files or non dynamic data
    dim(4) = floor(bytes*8/dtype/prod(dim(1:3)));
    fprintf('Size of REC file does not match PAR file. ');
    fprintf('Setting dim(4) to %d.\n',dim(4));
end
% Two scaling options [using DV here]:
% DV = displayed value on console
% FP = floating point value 
scale = hdr.ImageInfo(1).rescale_slope;
inter = hdr.ImageInfo(1).rescale_intercept;
mat   = spm_par_orient(hdr);
switch dtype
    case 8
        dtype = spm_type('int8'); % uint8?
    case 16
        dtype = spm_type('int16'); % uint16?
    case 32
        dtype = spm_type('float32');
    case 64
        dtype = spm_type('float64');
    otherwise
        error('Unknown data type.');
end

ofile    = spm_file(parfile,'path',opts.outdir,'ext',opts.ext);
dato     = file_array(ofile,dim,[dtype 0],0,scale,inter);

N        = nifti;
N.dat    = dato;
N.mat    = mat;
N.mat0   = mat;
N.mat_intent  = 'Scanner';
N.mat0_intent = 'Scanner';
N.descrip     = hdr.ProtocolName;
%N.timing = struct('toffset',[],'tspace',[]); % store TR
create(N);

%-Write NifTI data
%--------------------------------------------------------------------------
dati     = file_array(recfile,dim,[dtype 0],0,scale,inter);
dato     = N.dat;

% should handle interleaved data
% data should not be scaled/unscaled
for i=1:dim(4)
    slice_order = [hdr.ImageInfo([hdr.ImageInfo.dynamic_scan_number]==i).slice_number];
    dato(:,:,:,i) = dati(:,:,slice_order,i);
end


%==========================================================================
% FUNCTION hdr = spm_par_hdr(fname)
%==========================================================================
function hdr = spm_par_hdr(fname)

% could yse fileread instead
fid = fopen(fname,'rt');
if fid == -1
    error('Cannot open "%s".',fname);
end

iid = [];
j = 1;

while 1
    l = fgetl(fid);
    if ~ischar(l), break, end
    l = strtrim(l);
    if isempty(l), continue, end
    switch l(1)
        case '#'
            if strncmp(l,'# CLINICAL TRYOUT',17)
                hdr.Version = fliplr(strtok(fliplr(l)));
            end
            if strncmp(l,'# Dataset name:',15)
                hdr.DatasetName = strtrim(l(16:end));
            end
            if strncmp(l,'# === IMAGE INFORMATION DEFINITION',34)
                l = fgetl(fid); % #  The rest of this file contains...
                l = fgetl(fid); % # 
                iid = spm_par_get_iid(fid);
            end
        case '.'
            [key,val] = strtok(l(2:end),':');
            e = spm_par_findindict(strtrim(key));
            if ~isempty(e)
                hdr.(e.name) = feval(e.fcn,sscanf(strtrim(val(2:end)),e.fmt));
            else
                fprintf('Unknown key "%s" - ignored.\n',strtrim(key));
            end
        otherwise
            if isempty(iid)
                fclose(fid);
                error('Cannot read image information without definition.');
            end
            C = textscan(l,[iid.typ]);
            o = 1;
            for i=1:numel(iid)
                for k=1:iid(i).n
                    hdr.ImageInfo(j).(iid(i).name)(k) = C{o};
                    o = o + 1;
                end
            end
            j = j + 1;
    end
end

fclose(fid);


%==========================================================================
% FUNCTION iid = spm_par_get_iid(fid)
%==========================================================================
function iid = spm_par_get_iid(fid)
% get Image Information Definition

n = 1;
while 1
    l = fgetl(fid);
    if ~ischar(l), break, end
    l = strtrim(l(2:end));
    if isempty(l), break, end
    i = strfind(l,'(');
    if isempty(i), continue, else i = i(end); end
    
    name = strrep(strrep(strtrim(strtok(l(1:i-1),'(<')),' ','_'),'-','_');
    typ  = strtrim(l(i:end));
    [m,typ] = strtok(typ(2:end-1),'*');
    if isempty(typ)
        typ = m; m = 1;
    else
        m = str2double(m);
        typ = typ(2:end);
    end
    
    iid(n).name = name;
    switch typ
        case 'integer'
            iid(n).typ = strtrim(repmat('%d',1,m));
        case 'float'
            iid(n).typ = strtrim(repmat('%f',1,m));
        case 'string'
            iid(n).typ = strtrim(repmat('%s',1,m));
        otherwise
            error('Unknown data type "%s".',typ);
    end
    iid(n).n = m;
    n = n + 1;
end


%==========================================================================
% FUNCTION e = spm_par_findindict(key)
%==========================================================================
function e = spm_par_findindict(key)

dict = spm_par_dict;
for i=1:numel(dict)
    if strcmp(dict(i).key,key)
        e = dict(i);
        return
    end
end
e = [];


%==========================================================================
% FUNCTION dict = spm_par_dict 
%==========================================================================
function d = spm_par_dict

persistent dict;
if ~isempty(dict)
    d = dict;
    return;
end

n = @(x) x;
dict = {...
    'Patient name',                        'PatientName', '%s', n
    'Examination name',                    'ExaminationName', '%s', n
    'Protocol name',                       'ProtocolName', '%s', n
    'Examination date/time',               'ExaminationDate', '%s', @(x) datevec(x,'yyyy.mm.dd / HH:MM:SS')
    'Series Type',                         'SeriesType', '%s', n
    'Series_data_type',                    'SeriesDataType', '%s', n
    'Acquisition nr',                      'AcquisitionNr', '%d', n
    'Reconstruction nr',                   'ReconstructionNr', '%d', n
    'Scan Duration [sec]',                 'ScanDuration', '%d', n
    'Max. number of cardiac phases',       'MaxNumberCardiacPhases', '%d', n
    'Max. number of echoes',               'MaxNumberOfEchoes', '%d', n
    'Max. number of slices/locations',     'MaxNumberOfSlices', '%d', n
    'Max. number of dynamics',             'MaxNumberOfDynamics', '%d', n
    'Max. number of mixes',                'MaxNumberOfMixes', '%d', n
    'Patient Position',                    'PatientPosition', '%s', n
    'Patient position',                    'PatientPosition', '%s', n
    'Preparation direction',               'PreparationDirection', '%s', n
    'Technique',                           'Technique', '%s', n
    'Scan resolution  (x, y)',             'ScanResolution', '%d  %d', n
    'Scan mode',                           'ScanMode', '%s', n
    'Repetition time [ms]',                'RepetitionTime', '%f', n
    'Repetition time [msec]',              'RepetitionTime', '%f', n
    'FOV (ap,fh,rl) [mm]',                 'FOV', '%f  %f  %f', n
    'Water Fat shift [pixels]',            'WaterFatShift', '%f', n
    'Angulation midslice(ap,fh,rl)[degr]', 'AngulationMidslice', '%f  %f  %f', n
    'Off Centre midslice(ap,fh,rl) [mm]',  'OffCentreMidslice', '%f  %f  %f', n
    'Flow compensation <0=no 1=yes> ?',    'FlowCompensation', '%d', @logical
    'Presaturation     <0=no 1=yes> ?',    'Presaturation', '%d', @logical
    'Phase encoding velocity [cm/sec]',    'PhaseEncodingVelocity', '%f  %f  %f', n
    'MTC               <0=no 1=yes> ?',    'MTC', '%d', @logical
    'SPIR              <0=no 1=yes> ?',    'SPIR', '%d', @logical
    'EPI factor        <0,1=no EPI>',      'EPIFactor', '%d', n
    'Dynamic scan      <0=no 1=yes> ?',    'DynamicScan', '%d', @logical
    'Diffusion         <0=no 1=yes> ?',    'Diffusion', '%d', @logical
    'Diffusion echo time [ms]',            'DiffusionEchoTime', '%f', n
    'Diffusion echo time [msec]',          'DiffusionEchoTime', '%f', n
    'Max. number of diffusion values',     'MaxNumberOfDiffusionValues', '%d', n
    'Max. number of gradient orients',     'MaxNumberOfGradientsOrients', '%d', n
    'Number of label types   <0=no ASL>',  'NumberOfLabelTypes', '%d', n
    };

dict = struct(...
    'key' , dict(:,1),...
    'name', dict(:,2),...
    'fmt' , dict(:,3),...
    'fcn' , dict(:,4));

d = dict;


%==========================================================================
% FUNCTION M = spm_par_orient(hdr)
%==========================================================================
function M = spm_par_orient(hdr)

dim = double([hdr.ImageInfo(1).recon_resolution hdr.MaxNumberOfSlices]);
to_center = [eye(3) -dim(:)/2];
to_center(4,4) = 1;

vox = diag([hdr.ImageInfo(1).pixel_spacing hdr.ImageInfo(1).slice_thickness+hdr.ImageInfo(1).slice_gap 1]);

switch hdr.ImageInfo(1).slice_orientation
    case 1   % transversal
        to_pls = [0 1 0 0; 0 0 1 0; 1 0 0 0; 0 0 0 1];
    case 2   % sagittal
        to_pls = diag([1 -1 -1 1]);
    case 3   % coronal
        to_pls = [0 0 1 0; 0 -1 0 0; 1 0 0 0; 0 0 0 1];
    otherwise
        error('Unknown slice orientation.');
end

ang = hdr.ImageInfo(1).image_angulation * pi / 180; % {radians} [AP FH RL]
rz  = [cos(ang(3)) -sin(ang(3)) 0; sin(ang(3)) cos(ang(3)) 0; 0 0 1];
ry  = [cos(ang(2)) 0 sin(ang(2)); 0 1 0; -sin(ang(2)) 0 cos(ang(2))];
rx  = [1 0 0; 0 cos(ang(1)) -sin(ang(1)); 0 sin(ang(1)) cos(ang(1))];
rot = rx * ry *rz;
rot(4,4) = 1;

pls = rot * to_pls * vox * to_center;
pls(1:3,4) = pls(1:3,4) + hdr.ImageInfo(1).image_offcentre(:);

pls_to_ras = [0 0 -1 0; -1 0 0 0; 0 1 0 0; 0 0 0 1];
          
M = pls_to_ras * pls;
