function out = FieldMap_applyvdm(job)
% Apply VDM and reslice images 
% FORMAT FieldMap_applyvdm(job)
% job.data(sessnum).scans   - images for session/run sessnum
% job.data(sessnum).vdmfile - VDM file for session/run sessnum      
% job.roptions.rinterp - interpolation method
% job.roptions.wrap    - perform warp around in specified dimensions
% job.roptions.mask    - perform masking
% job.roptions.which(1) - reslice images in time series only
% job.roptions.which(2) - reslice images in time series and mean
% job.roptions.prefix   - prefix for vdm applied files
% job.roptions.pedir    - phase encode direction (i.e. aplly vdm file along
% this dimension
%__________________________________________________________________________
%
% A VDM (voxel displacement map) created using the FieldMap toolbox 
% can be used to resample and reslice realigned images to the original 
% subdirectory with the same (prefixed) filename. 
%
% Voxels in the images will be shifted according to the values in the VDM 
% file along the direction specified by job.roptions.pedir (i.e. this is
% usually the phase encode direction) and resliced to the space of the 
% first image in the time series.
%
% Inputs:
% A job structure containing fields for the input data and the processing
% options. The input data contains the series of images conforming to 
% SPM data format (see 'Data Format'), the relative displacement of the images 
% is stored in their header and a VDM which has (probably) been created 
% using the FieldMap toolbox and matched to the first image in the time 
% series (this can also be done via the FieldMap toolbox).
%
% Outputs:
% The resampled and resliced images resliced to the same subdirectory with a prefix.
%__________________________________________________________________________
% Copyright (C) 2011-2014 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton
% $Id: FieldMap_applyvdm.m 6258 2014-11-07 18:15:40Z guillaume $

tiny = 5e-2;

% assemble roptions
%--------------------------------------------------------------------------
flags.interp    = job.roptions.rinterp;
flags.wrap      = job.roptions.wrap;
flags.mask      = job.roptions.mask;
flags.which     = job.roptions.which(1);
flags.mean      = job.roptions.which(2);
flags.prefix    = job.roptions.prefix;
flags.pedir     = job.roptions.pedir;
hold = [repmat(flags.interp,1,3) flags.wrap];

% Determine dimension along which to apply vdm
applydim=flags.pedir;
if applydim~=1 &&  applydim~=2 && applydim~=3
    applydim=2;
end

% Gather up data into ds structure which holds images and vdm file
%--------------------------------------------------------------------------
P = {};
for i = 1:numel(job.data)
    P{i} = strvcat(job.data(i).scans{:});
    ds(i).P=spm_vol(P{i});
    if ~isempty(job.data(i).vdmfile)
        sfP{i} = job.data(i).vdmfile{1};
        ds(i).sfP=spm_vol(sfP{i}); 
    else
        sfP{i} = [];       
    end  
    ds(i).hold = [1 1 1 0 1 0];
end

ntot = 0;
for i=1:length(ds)
   ntot = ntot + length(ds(i).P);
end

% Set up  x y z for resampling
%--------------------------------------------------------------------------
[x,y,z] = ndgrid(1:ds(1).P(1).dim(1),1:ds(1).P(1).dim(2),1:ds(1).P(1).dim(3));
xyz = [x(:) y(:) z(:) ones(prod(ds(1).P(1).dim(1:3)),1)]; clear x y z;

% Create mask if required (usually default and required to create mean)
%--------------------------------------------------------------------------
if flags.mask || flags.mean,
    spm_progress_bar('Init',ntot,'Computing available voxels',...
        'volumes completed');
    
    if flags.mean
        Count    = zeros(prod(ds(1).P(1).dim(1:3)),1);
        Integral = zeros(prod(ds(1).P(1).dim(1:3)),1);
    end
    
    % if flags.mask
    msk = zeros(prod(ds(1).P(1).dim(1:3)),1);
    % end
    
    % To create mean, read each session specific vdmfile in
    % to the space of the first image of first session
    tv = 1;
    for s=1:length(ds)
        T = ds(s).sfP.mat\ds(1).P(1).mat;
        txyz = xyz * T';
        c = spm_bsplinc(ds(s).sfP,ds(1).hold);
        ds(s).sfield = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds(1).hold);
        ds(s).sfield = ds(s).sfield(:);
        clear c txyz;
        sess_msk = zeros(prod(ds(s).P(1).dim(1:3)),1);
        % Read in each images in space of first image of first session
        for i = 1:numel(ds(s).P)
            T = inv(ds(s).P(i).mat) * ds(1).P(1).mat;
            txyz = xyz * T';
            txyz(:,applydim) = txyz(:,applydim) + ds(s).sfield;
            tmp       = false(size(txyz,1),1);
            if ~flags.wrap(1), tmp = tmp | txyz(:,1) < (1-tiny) | txyz(:,1) > (ds(s).P(i).dim(1)+tiny); end
            if ~flags.wrap(2), tmp = tmp | txyz(:,2) < (1-tiny) | txyz(:,2) > (ds(s).P(i).dim(2)+tiny); end
            if ~flags.wrap(3), tmp = tmp | txyz(:,3) < (1-tiny) | txyz(:,3) > (ds(s).P(i).dim(3)+tiny); end
            sess_msk = sess_msk + real(tmp);
            spm_progress_bar('Set',tv);
            tv = tv+1;
        end
        msk = msk + sess_msk;
        if flags.mean, Count = Count + repmat(length(ds(s).P),prod(ds(s).P(1).dim(1:3)),1) - sess_msk; end
        
        %
        % Include static field in estmation of mask.
        %
        if isfield(ds(s),'sfP') && ~isempty(ds(s).sfP)
            T = inv(ds(1).sfP.mat) * ds(1).P(1).mat;
            txyz = xyz * T';
            tmp  = false(size(txyz,1),1);
            if ~flags.wrap(1), tmp = tmp | txyz(:,1) < (1-tiny) | txyz(:,1) > (ds(1).sfP.dim(1)+tiny); end
            if ~flags.wrap(2), tmp = tmp | txyz(:,2) < (1-tiny) | txyz(:,2) > (ds(1).sfP.dim(2)+tiny); end
            if ~flags.wrap(3), tmp = tmp | txyz(:,3) < (1-tiny) | txyz(:,3) > (ds(1).sfP.dim(3)+tiny); end
            msk = msk + real(tmp);
        end

        if isfield(ds(1),'sfield') && ~isempty(ds(1).sfield)
            ds(1).sfield = [];
        end
    end
    if flags.mask, msk = find(msk ~= 0); end
end

% Apply fieldmap to all files, looping through sessions
%--------------------------------------------------------------------------
spm_progress_bar('Init',ntot,'Reslicing','volumes completed');
tv = 1;
for s = 1:numel(ds)
    
    % Get transformation between distortion field and first image
    T = ds(s).sfP.mat\ds(1).P(1).mat;
    txyz = xyz * T';
    c = spm_bsplinc(ds(s).sfP,ds(1).hold);
    ds(s).sfield = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),ds(1).hold);
    ds(s).sfield = ds(s).sfield(:);
    
    % Read in each images in space of first image of first session
    for i = 1:numel(ds(s).P)
        T = inv(ds(s).P(i).mat) * ds(1).P(1).mat;
        txyz = xyz * T';
        txyz(:,applydim) = txyz(:,applydim) + ds(s).sfield;       
        c = spm_bsplinc(ds(s).P(i),hold);
        ima = spm_bsplins(c,txyz(:,1),txyz(:,2),txyz(:,3),hold);
        %
        % Write out resliced images
        %
        if flags.which
            PO         = ds(s).P(i);
            [pth,nm,xt,vr] = spm_fileparts(deblank(PO.fname));
            PO.fname       = fullfile(pth,[flags.prefix nm xt vr]);
            PO.mat     = ds(1).P(1).mat;
            PO.descrip = sprintf('spm - applied vdm');
            ivol       = ima;
            if flags.mask
                ivol(msk) = NaN;
            end
            ivol = reshape(ivol,PO.dim(1:3));
            PO   = spm_create_vol(PO);
            spm_write_vol(PO,ivol);
            if nargout > 0
                  out.sess(s).rfiles{i} = PO.fname;   
            end
        end
        %
        % Build up mean image if so required.
        %
        if flags.mean
            Integral = Integral + nan2zero(ima);
        end
        spm_progress_bar('Set',tv);
        tv = tv+1;
    end
    if isfield(ds(s),'sfield') && ~isempty(ds(s).sfield)
        ds(s).sfield = [];
    end
end

% Write mean image
%--------------------------------------------------------------------------
if flags.mean
   % Write integral image (16 bit signed)
   %-----------------------------------------------------------------------
   sw = warning('off','MATLAB:divideByZero'); 
   Integral   = Integral./Count;
   warning(sw);
   PO         = ds(1).P(1);
   [pth,nm,xt,vr] = spm_fileparts(deblank(ds(1).P(1).fname));
   PO.fname       = fullfile(pth,['mean' flags.prefix nm xt vr]);
   PO.pinfo   = [max(max(max(Integral)))/32767 0 0]';
   PO.descrip = 'spm - mean applied vdm image';
   PO.dt      = [spm_type('int16') spm_platform('bigend')];
   ivol = reshape(Integral,PO.dim);
   spm_write_vol(PO,ivol);
end

if nargout > 0
   out.rmean{1} = PO.fname;
end

spm_figure('Clear','Interactive');


%==========================================================================
function vo = nan2zero(vi)
vo = vi;
vo(~isfinite(vo)) = 0;
return;
