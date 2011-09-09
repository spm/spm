function R = spm_searchlight(SPM,searchopt,fun,varargin)
% Local mass-multivariate (c.f., searchlight) facility
% FORMAT R = spm_searchlight(SPM,searchopt,fun,varargin)
% SPM       - structure with fields:
%    .xY.VY - filenames char array or spm_vol struct array of images
%    .VM    - filename or spm_vol structure to a mask (binary) image
%             Mask image can have any orientation, voxel size or data type.
%             It is interpolated using nearest neighbour interpolation to
%             the voxel locations of the data.
%             If empty, all voxels are used.
% searchopt - searchlight options using VOI structure (xY) from spm_ROI
%    .def   - searchlight definition {['sphere'] 'box'}
%    .spec  - searchlight parameters [sphere radius {mm}]
% fun       - function handle to a function that takes three input arguments:
%               a [n x v] matrix (nb images x nb voxels within searchlight)
%               a [3 x v] matrix of voxels location within searchlight {vox}
%               a list of parameters provided in varargin
%             and returns a vector value [1 x N]
% varargin  - list of parameters sent to fun
%
% R         - a [N x 1] cell array with each output (fun nargout) reshaped
%             to a volume or directly a volume if N == 1
%             Values outside the mask are attributed NaN.
%__________________________________________________________________________
%
% References:
%
% [1] Adaptive Analysis of fMRI Data. Friman O, Borga M, Lundberg P and 
% Knutsson H. (2003) NeuroImage 19(3):837-845.
%
% [2] Information-based functional brain mapping. Kriegeskorte N, Goebel R,
% Bandettini P. (2006) PNAS 103: 3863-3868.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_searchlight.m 4475 2011-09-09 17:53:14Z guillaume $

spm('FnBanner',mfilename);                                              %-#

%-Get input images
%--------------------------------------------------------------------------
try
    VY = SPM.xY.VY;
catch
    [VY, sts] = spm_select([1 Inf],'image','Select images');
    if ~sts, R = {}; return; end
end
if iscellstr(VY), VY = char(VY); end
if ~isstruct(VY)
    VY = spm_vol(VY);
end

%-Check dimensions and orientations of all images
%--------------------------------------------------------------------------
spm_check_orientations(VY);

%-Get mask image
%--------------------------------------------------------------------------
try
    VM = SPM.VM;
catch
    [VM, sts] = spm_select([0 1],'image','Select mask');
    if ~sts, VM = []; end
end
if ~isstruct(VM) && ~isempty(VM)
    VM = spm_vol(VM);
end

%-Get space details
%--------------------------------------------------------------------------
N            = numel(VY);                          %-number of images
M            = VY(1).mat;                          %-voxels to mm matrix
iM           = inv(M);                             %-mm to voxels matrix
DIM          = VY(1).dim;                          %-image dimensions
NDIM         = prod([DIM N]);                      %-overall dimension
[x,y,z]      = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ          = [x(:)';y(:)';z(:)']; clear x y z    %-voxel coordinates {vx}
XYZmm        = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];%-voxel coordinates {mm}
XYZmm_cpy    = XYZmm;                              %-copy without masking

%-Strategy for reading data: memory [1] / block [2] / disk [3]
%--------------------------------------------------------------------------
fprintf('%-40s: ','Read data');                                         %-#
try
    ds         = NDIM*8;                           %-data size {bytes}
    MAXMEM     = spm('memory');                    %-max amnt usable memory
    if ds > MAXMEM
        YY     = struct('y',[], 'i',[Inf Inf]);
        blk    = floor(MAXMEM/(prod([DIM(1:2) N])*8));%-blk size {# slices}
        if blk == 0, error('revert to disk.'); end
        Ystrtg = 2;
        fprintf('%30s\n','...per block');                               %-#
    else
        YY     = spm_read_vols(VY);
        Ystrtg = 1;
        fprintf('%30s\n','...in memory');                               %-#
    end
catch
    Ystrtg     = 3;
    fprintf('%30s\n','...from disk');                                   %-#
end

%-Search volume (from mask)
%--------------------------------------------------------------------------
fprintf('%-40s: ','Read mask');                                         %-#
if ~isempty(VM)
    if any(DIM-VM.dim) || any(any(abs(VM.mat-M)>1e-4))
        MM   = spm_get_data(VM,VM.mat\[XYZmm;ones(1,size(XYZmm,2))],false);
    else
        MM   = spm_read_vols(VM);
    end
    MM       = logical(MM);
    XYZmm    = XYZmm(:,MM(:));
    XYZ      = XYZ(:,MM(:));
    fprintf('%30s\n','...done');                                        %-#
else
    MM       = true(DIM);
    fprintf('%30s\n', '...none');                                       %-#
end

%-Searchlight options (clique definition)
%--------------------------------------------------------------------------
try, xY      = searchopt; end
xY.xyz       = [NaN NaN NaN];
xY.rej       = {'cluster','mask'};
xY           = spm_ROI(xY);

%-Evaluated function
%--------------------------------------------------------------------------
if ischar(fun)
    fun      = str2func(fun);
end
if ~isa(fun, 'function_handle')
    error('''fun'' must be a function handle with two input parameters.');
end

%-Get local clique and perform searchlight over voxels
%==========================================================================

%-Build local clique
%--------------------------------------------------------------------------
fprintf('%-40s: ','Construct clique');                                  %-#
c            = round(DIM(:)/2);
xY.xyz       = M(1:3,:) * [c;1];
[xY, clique] = spm_ROI(xY,XYZmm_cpy);
clique       = round(iM(1:3,:) * [clique;ones(1,size(clique,2))]);
clique       = bsxfun(@minus, clique, c);
dc           = (max(clique,[],2) - min(clique,[],2) + 1)';
fprintf('%30s\n',sprintf('%d voxels - [%dx%dx%d]',size(clique,2),dc));  %-#
if Ystrtg == 2 && blk < dc(3)
    fprintf('%-40s: %30s\n','Read data','revert to disk');              %-#
    Ystrtg   = 3;
end

%-Initialise progress bar
%--------------------------------------------------------------------------
spm_figure('GetWin','Interactive');
spm_progress_bar('Init',size(XYZ,2));
Ibar = floor(linspace(1,size(XYZ,2), 100));
fprintf('%-40s: %30s','Searchlight','...computing');                    %-#

SLR  = [];

%-Searchlight
%--------------------------------------------------------------------------
for i=1:size(XYZ,2)

    %-Local clique (handle image boundaries and mask)
    %----------------------------------------------------------------------
    xyz          = bsxfun(@plus,XYZ(:,i),clique);
    xyz(:,any(bsxfun(@lt,xyz,[1 1 1]') | bsxfun(@gt,xyz,DIM'))) = [];
    idx          = sub2ind(DIM,xyz(1,:),xyz(2,:),xyz(3,:));
    j            = MM(idx);
    idx          = idx(j);
    xyz          = xyz(:,j);
    
    %-Read data
    %----------------------------------------------------------------------
    if Ystrtg == 3
        Y        = spm_get_data(VY,xyz,false);
    elseif Ystrtg == 2
        if min(xyz(3,:)) < min(YY.i) || max(xyz(3,:)) > max(YY.i)
            YY.i = min(xyz(3,:)):min(xyz(3,:))+blk;
            YY.i(YY.i>DIM(3)) = [];
            YY.y = zeros(DIM(1),DIM(2),numel(YY.i),N);
            for v=1:N
                for p=1:numel(YY.i)
                    YY.y(:,:,p,v) = ...
                        spm_slice_vol(VY(v), ...
                        spm_matrix([0 0 YY.i(p)]), ...
                        VY(v).dim(1:2),0);
                end
            end
        end
        idx      = idx - (min(YY.i)-1) * prod(DIM(1:2));
        k        = prod(DIM(1:2))*numel(YY.i);
        idx      = bsxfun(@plus, idx(:), 0:k:k*N-1);
        Y        = YY.y(idx)';
    elseif Ystrtg == 1
        idx      = bsxfun(@plus, idx(:), 0:prod(DIM):NDIM-1);
        Y        = YY(idx)';
    end
    
    %-Call user-specified function
    %----------------------------------------------------------------------
    if isempty(SLR)
        t        = fun(Y,xyz,varargin{:});
        SLR      = zeros(size(XYZ,2),length(t));
        SLR(1,:) = t;
    else
        SLR(i,:) = fun(Y,xyz,varargin{:});
    end
    
    %-Update progress bar
    %----------------------------------------------------------------------
    if any(Ibar == i), spm_progress_bar('Set', i); end
    
end

%-Clear progress bar
%--------------------------------------------------------------------------
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#
spm_progress_bar('Clear');

%-Return computations - reshaped as volumes in a cell array
%--------------------------------------------------------------------------
R        = cell(size(SLR,2),1);
for i=1:size(SLR,2)
    MV   = NaN(DIM);
    MV(sub2ind(DIM,XYZ(1,:),XYZ(2,:),XYZ(3,:))) = SLR(:,i);
    R{i} = MV;
end

%-Write images if required
%--------------------------------------------------------------------------
fprintf('%-40s: %30s','Output images','...writing');                    %-#
VO(1:size(SLR,2)) = deal(struct(...
    'fname',   [],...
    'dim',     DIM,...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', ['spm_searchlight: ' func2str(fun)]));

for i=1:size(SLR,2)
    VO(i).fname   = sprintf('%s_%04d%s', ...
        'searchlight',i,spm_file_ext);
    VO(i).descrip = sprintf('%s (%04d)',VO(i).descrip,i);
end
for i=1:size(SLR,2)
    spm_write_vol(VO(i),R{i});
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#

%-And exit
%--------------------------------------------------------------------------
if size(SLR,2) == 1, R = R{1}; end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
