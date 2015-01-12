function out = spm_run_voi(job)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_run_voi.m 6301 2015-01-12 17:23:08Z guillaume $


%-Load SPM.mat
%--------------------------------------------------------------------------
swd     = spm_file(job.spmmat{1},'fpath');
load(fullfile(swd,'SPM.mat'));
SPM.swd = swd;

%-Initialise VOI voxels coordinates
%--------------------------------------------------------------------------
[x,y,z] = ndgrid(1:SPM.xVol.DIM(1),1:SPM.xVol.DIM(2),1:SPM.xVol.DIM(3));
XYZ     = [x(:),y(:),z(:)]'; clear x y z
XYZmm   = SPM.xVol.M(1:3,:) * [XYZ;ones(1,size(XYZ,2))];

%-Estimate VOIs
%--------------------------------------------------------------------------
voi     = cell(1,numel(job.roi));
for i=1:numel(job.roi)
    voi = roi_estim(XYZmm,i,job,SPM,voi);
end

%-Evaluate resulting VOI
%--------------------------------------------------------------------------
voi     = roi_eval(voi,job.expression);

%-Save VOI as image
%--------------------------------------------------------------------------
Vm = struct('fname', fullfile(swd, ['VOI_' job.name '_mask' spm_file_ext]), ...
     'dim',     SPM.xVol.DIM', ...
     'dt',      [spm_type('uint8') spm_platform('bigend')], ...
     'mat',     SPM.xVol.M, ...
     'pinfo',   [1 0 0]', ...
     'descrip', 'VOI');
Vm = spm_write_vol(Vm,voi);

%-Extract VOI time-series
%--------------------------------------------------------------------------
xY.name    = job.name;
xY.Ic      = job.adjust;
xY.Sess    = job.session;
xY.xyz     = []'; % irrelevant here
xY.def     = 'mask';
xY.spec    = Vm;

xSPM.XYZmm = XYZmm;
xSPM.XYZ   = XYZ;
xSPM.M     = SPM.xVol.M; % irrelevant here

if ~isempty(xY.Ic), cwd = pwd; cd(SPM.swd); end % to find beta images
[Y,xY]     = spm_regions(xSPM,SPM,[],xY);
if  ~isempty(xY.Ic), cd(cwd); end

%-Save first eigenimage
%--------------------------------------------------------------------------
Ve = struct('fname', fullfile(swd, ['VOI_' job.name '_eigen' spm_file_ext]), ...
     'dim',     SPM.xVol.DIM', ...
     'dt',      [spm_type('float32') spm_platform('bigend')], ...
     'mat',     SPM.xVol.M, ...
     'pinfo',   [1 0 0]', ...
     'descrip', 'VOI: first eigenimage');
Ve = spm_create_vol(Ve);
eigimg = double(voi);
eigimg(voi) = xY.v;
Ve = spm_write_vol(Ve,eigimg);

%-Export results
%--------------------------------------------------------------------------
assignin('base','Y',Y);
assignin('base','xY',xY);

if isfield(SPM,'Sess'), s = sprintf('_%i',xY.Sess); else s = ''; end
out.voimat = cellstr(fullfile(swd,['VOI_' job.name s '.mat']));
out.voiimg = cellstr(Vm.fname);
out.voieig = cellstr(Ve.fname);

%==========================================================================
function voi = roi_estim(xyz,n,job,SPM,voi)

if ~isempty(voi{n}), return; end

voi{n} = false(SPM.xVol.DIM');

Q      = ones(1,size(xyz,2));

switch char(fieldnames(job.roi{n}))
    
    case 'sphere'
    %----------------------------------------------------------------------
        c              = get_centre(xyz,n,job,SPM,voi);
        r              = job.roi{n}.sphere.radius;
        voi{n}(sum((xyz - c*Q).^2) <= r^2) = true;
    
    case 'box'
    %----------------------------------------------------------------------
        c              = get_centre(xyz,n,job,SPM,voi);
        d              = job.roi{n}.box.dim(:);
        voi{n}(all(abs(xyz - c*Q) <= d*Q/2)) = true;
    
    case 'spm'
    %----------------------------------------------------------------------
        if isempty(job.roi{n}.spm.spmmat{1})
            job.roi{n}.spm.spmmat = job.spmmat;
        end
        [SPM1,xSPM]    = getSPM(job.roi{n}.spm);
        voi1           = zeros(SPM1.xVol.DIM');
        voi1(sub2ind(SPM1.xVol.DIM',xSPM.XYZ(1,:),xSPM.XYZ(2,:),xSPM.XYZ(3,:))) = 1;
        XYZ            = SPM1.xVol.iM(1:3,:)*[xyz; Q];
        voi{n}(spm_sample_vol(voi1, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0) = true;
    
    case 'mask'
    %----------------------------------------------------------------------
        v              = spm_vol(job.roi{n}.mask.image{1});
        t              = job.roi{n}.mask.threshold;
        iM             = inv(v.mat);
        XYZ            = iM(1:3,:)*[xyz; Q];
        voi{n}(spm_sample_vol(v, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > t) = true;
    
    case 'label'
    %----------------------------------------------------------------------
        v              = spm_vol(job.roi{n}.label.image{1});
        l              = job.roi{n}.label.list;
        iM             = inv(v.mat);
        XYZ            = iM(1:3,:)*[xyz; Q];
        voi{n}(ismember(spm_sample_vol(v, XYZ(1,:), XYZ(2,:), XYZ(3,:),0),l)) = true;
    
end

%==========================================================================
function voi = roi_eval(voi,expr)
for i=1:numel(voi)
    eval(sprintf('i%d=voi{%d};',i,i));
end
try
    eval(['voi=' expr ';']);
catch
    error('The expression cannot be evaluated.');
end

%==========================================================================
function idx = roi_expr(expr)
e   = regexp(expr,'i\d+','match');
idx = zeros(1,numel(e));
for i=1:numel(e)
    idx(i) = str2num(e{i}(2:end));
end

%==========================================================================
function [SPM, xSPM] = getSPM(s)
xSPM.swd       = spm_file(s.spmmat{1},'fpath');
xSPM.Ic        = s.contrast;
xSPM.n         = s.conjunction;
xSPM.u         = s.thresh;
xSPM.thresDesc = s.threshdesc;
xSPM.k         = s.extent;
xSPM.title     = '';
xSPM.Im        = [];
if ~isempty(s.mask)
    xSPM.Im    = s.mask.contrast;
    xSPM.pm    = s.mask.thresh;
    xSPM.Ex    = s.mask.mtype;
end
[SPM,xSPM]     = spm_getSPM(xSPM);

%==========================================================================
function c = get_centre(xyz,n,job,SPM,voi)
t             = char(fieldnames(job.roi{n}));
c             = job.roi{n}.(t).centre(:);
mv            = char(fieldnames(job.roi{n}.(t).move));
if strcmp(mv,'fixed'), return; end

m             = job.roi{n}.(t).move.(mv).spm;
e             = job.roi{n}.(t).move.(mv).mask;
k             = union(roi_expr(e), m);
for i=1:numel(k)
    voi       = roi_estim(xyz,k(i),job,SPM,voi);
end
try
    if isempty(job.roi{m}.spm.spmmat{1})
        job.roi{m}.spm.spmmat = job.spmmat;
    end
catch
    error('The SPM index does not correspond to a Thresholded SPM ROI.');
end
[mySPM, xSPM] = getSPM(job.roi{m}.spm);
XYZmm         = xSPM.XYZmm;
XYZ           = SPM.xVol.iM(1:3,:)*[XYZmm;ones(1,size(XYZmm,2))];
Z             = xSPM.Z;
if ~isempty(e)
    R         = spm_sample_vol(uint8(roi_eval(voi,e)), ...
                  XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0;
    XYZ       = XYZ(:,R);
    XYZmm     = xSPM.XYZmm(:,R);
    Z         = xSPM.Z(R);
end
[N, Z, M]     = spm_max(Z,XYZ);
if isempty(Z)
    warning('No voxel survived. Default to user-specified centre.');
    return
end

str           = '[%3.0f %3.0f %3.0f]';
switch mv
    case 'global'
        [i,j] = max(Z);
        nc    = SPM.xVol.M(1:3,:)*[M(:,j);1];
        str   = sprintf(['centre moved to global maximum ' str],nc);
    case 'local'
        XYZmm = SPM.xVol.M(1:3,:)*[M;ones(1,size(M,2))];
        nc    = spm_XYZreg('NearestXYZ',c,XYZmm);
        str   = sprintf(['centre moved from ' str ' to ' str],c,nc);
    case 'supra'
        nc    = spm_XYZreg('NearestXYZ',c,XYZmm);
        str   = sprintf(['centre moved from ' str ' to ' str],c,nc);
    otherwise
        error('Unknown option: ''%s''.',mv);
end
c             = nc;
fprintf(['   ' upper(t(1)) t(2:end) ' ' str '\n']);                     %-#
