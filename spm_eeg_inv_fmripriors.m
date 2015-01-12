function D = spm_eeg_inv_fmripriors(S)
% Generate fMRI priors for the M/EEG source reconstruction
% FORMAT D = spm_eeg_inv_fmripriors(S)
%
% S        - optional input struct
% (optional) fields of S:
%  .D      - MEEG object or filename of M/EEG mat-file
%  .fmri   - filename of prior (SPM) image to be used
%  [.gm    - filename of grey matter (GM) image] {unused}
%  .space  - native (0) or MNI (1) space (must be same for SPM and GM images)
%  .hthr   - height threshold of prior image [defaults: 0.5]
%  .ethr   - extent threshold of clusters in prior image [default: 1]
%  .ncomp  - maximal number of priors component to be extracted [default: Inf]
%  .smooth - variance of the smoothing kernel onto the surface [default: 0.2] {unused}
%  .disp   - whether to display priors on mesh [default: 0]
%
% D.inv{D.val}.inverse.fmri.priors   - MAT filename containing a variable 'pQ' that
%            is a [ncomp] cell array of [nb vertices] vectors describing spatial priors
% D.inv{D.val}.inverse.fmri.texture  - GIfTI texture filename containing all
%            spatial priors
% D.inv{D.val}.inverse.fmri.clusters - image filename containing clusters as labels
%__________________________________________________________________________
%
% Reference:
%
% A Parametric Empirical Bayesian framework for fMRI-constrained MEG/EEG 
% source reconstruction. Henson R, Flandin G, Friston K & Mattout J.
% Human Brain Mapping (in press).
%__________________________________________________________________________
% Copyright (C) 2008-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin and Rik Henson
% $Id: spm_eeg_inv_fmripriors.m 6301 2015-01-12 17:23:08Z guillaume $

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts]   = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

[D, val] = spm_eeg_inv_check(D);

%-Input parameters
%--------------------------------------------------------------------------
try
    S.fmri;
catch
    [S.fmri, sts] = spm_select(1,'image','Select prior (eg thresholded SPM) image...');
    if ~sts, return; end
end

% try
%     S.gm;
% catch
%     [S.gm, sts] = spm_select([0 1],'image','Select grey matter (GM) image...');
% end

try
    space = S.space;   % 0 = native, 1 = MNI
catch
    space = spm_input('Image space','+1','b',{'Native|MNI'},[0 1],1);
end

if space 
    try
        m = D.inv{val}.mesh.tess_mni;
        %fprintf('Note: assuming SPM and GM images are in MNI space...\n');
    catch
        error('This D structure has no MNI cortical mesh stored.');
    end
else
    try
        m = export(gifti(D.inv{val}.mesh.tess_ctx),'spm');
        %fprintf('Note: assuming SPM and GM images are in subject''s native MRI space...\n');
    catch
        error('This D structure has no cortical mesh stored.');
    end  
end

try
    S.hthr;
catch
    S.hthr = 0.5; % assume binary prior
end

try
    S.ethr;
catch
    S.ethr = 1; % no threshold on extent
end

try
    S.ncomp;
catch
    S.ncomp = Inf; % all components
end

try
    S.bincomp;
catch
    S.bincomp = 1; % default to binary priors
end

try
    S.varcomp;
catch
    S.varcomp = 1; % default to variance priors (vectors)
end

% try
%     S.smooth;
% catch
%     S.smooth = 0.2; 
% end

try
    S.disp;
catch
    S.disp = 0; 
end

%-Extracting clusters from functional image
%==========================================================================
V     = spm_vol(S.fmri);
prior = spm_read_vols(V);

%-Height threshold
%--------------------------------------------------------------------------
prior = prior > S.hthr;

%-Connected Component labelling
%--------------------------------------------------------------------------
[l2, num] = spm_bwlabel(double(prior),26);
if ~num
    fprintf('No suprathreshold clusters available.\n');
    return
end

%-Extent threshold, and sort clusters according to their extent
%--------------------------------------------------------------------------
[n, ni] = sort(histc(l2(:),0:num), 1, 'descend');
l  = zeros(size(l2));
n  = n(2:end);      ni = ni(2:end)-1;
ni = ni(n>=S.ethr); n  = n(n>=S.ethr);
S.ncomp = min(S.ncomp, length(n));
for i=1:S.ncomp
    l(l2==ni(i)) = i;
end
clear l2 ni
fprintf('Selected %d clusters (out of %d) in prior image.\n',S.ncomp,num);

%-Projecting volumetric clusters on surface mesh
%==========================================================================
q = zeros(S.ncomp, size(m.vert,1));
for i=1:S.ncomp
    q(i,:) = spm_mesh_project(m.vert,struct('dat',double(l==i),'mat',V.mat),'nn');
end
q(~any(q,2),:) = [];
fprintf('After projection, %d clusters remaining.\n',size(q,1));
if isempty(q), return; end

%-Smooth, binarize and save in output variable
%--------------------------------------------------------------------------
pQ = cell(1,size(q,1));
for i = 1:size(q,1)
    qq    = q(i,:)';
    %qq   = spm_mesh_smooth(struct('faces',double(m.face),'vertices',m.vert),qq, S.smooth);
    qq    = qq .* (qq > exp(-8));
    if S.bincomp
        pQ{i} =  double(qq > 0)';  % binarise
    end
end

%-Display and export clusters
%==========================================================================
if S.disp
    for i=1:numel(pQ)
        spm_eeg_render(struct('faces',double(m.face),'vertices',m.vert),...
            struct('texture',pQ{i}));
    end
end

%-Save clusters as an image of labels
%--------------------------------------------------------------------------
[pth,name] = fileparts(S.fmri);
D.inv{val}.inverse.fmri.clusters = fullfile(pth,['cluster_' name spm_file_ext]);
V = struct('fname',   D.inv{val}.inverse.fmri.clusters, ...
           'dim',     V.dim, ...
           'dt',      [spm_type('uint16') spm_platform('bigend')], ...
           'mat',     V.mat, ...
           'pinfo',   [1 0 0]', ...
           'descrip', 'clusters');
V = spm_write_vol(V,l);

%-Save spatial priors vectors as GIfTI file
%--------------------------------------------------------------------------
[pth,name] = fileparts(D.fname);
D.inv{val}.inverse.fmri.texture = fullfile(pth,['priors_' name '_' num2str(val) '.func.gii']);
G          = gifti;
G.cdata    = cat(1, pQ{:})';
save(G,D.inv{val}.inverse.fmri.texture);

%-Save spatial priors vectors as MAT-file
%--------------------------------------------------------------------------
[pth,name] = fileparts(D.fname);
D.inv{val}.inverse.fmri.priors = fullfile(pth,['priors_' name '_' num2str(val) '.mat']);
if ~S.varcomp
    pQ = pQ*pQ';
end

save(D.inv{val}.inverse.fmri.priors,'pQ', spm_get_defaults('mat.format'));

%-Save D structure
%--------------------------------------------------------------------------
%D.save;
