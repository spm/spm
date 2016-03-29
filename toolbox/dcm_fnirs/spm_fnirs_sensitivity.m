function A = spm_fnirs_sensitivity(DCM)
% Calculate sensitivity matrix which corresponds to the effective
% pathlength of detected photons for the channel measurements in the
% hemodynamic source. 
% FORMAT [A] = spm_fnirs_sensitivity(DCM)
% 
% DCM - DCM structure or its filename
%
% A - sensitivity matrix 
% 
% Green's function (see \dcm_fnirs\mmclab\estimate_greens_mmclab.m)
%--------------------------------------------------------------------------
% G.s - estimated Green's function from sensor (light emitter) positions
% into source positions [# sensor x # voxels x # wavelengths] 
% G.d - estimated Green's function from sensor (light detector) positions
% into source positions [# sensor x # voxels x # wavelengths] 
% G.xyz - MNI locations [3 x # voxels] 
% G.elem - tissue types of voxels [3 x # voxels] 
% 1-scalp, 2-CSF, 3-gray matter, 4-white matter 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $id$

%--------------------------------------------------------------------------
if ~nargin || isempty(DCM) 
    [DCM, sts] = spm_select(1, '^DCM*.*\.mat$', 'Select DCM.mat');
    if ~sts, S = []; return; end 
elseif ischar(DCM) 
    load(DCM); 
end

%--------------------------------------------------------------------------
load(DCM.Y.P.fname.pos); % load channel positions 
load(DCM.Y.P.fname.g);   % load Green's function 

xyzh = [DCM.xY.xyz];     % hemodynamic source positions 
xyzd = R.d.xyzH;         % detector positions 
ns = R.s.ns;             % number of sources (emitter)
nd = R.d.nd;             % number of detectors 
nh = size(xyzh, 2);      % number of hemodynamic sources 

nnode = size(G.xyz, 2);  % number of nodes 
r = DCM.options.rs .* 3; % maximum distance between point and distributed sources 
rois = DCM.Y.P.rois; 

%--------------------------------------------------------------------------
%-Find indices of Green's function of detector and hemodynamic source 
node_d = zeros(nd, 1);
for i = 1:nd
    diff = G.xyz - xyzd(:, i) * ones(1, nnode);
    diff = sum(diff.^2);
    node_d(i, 1) = find(diff == min(diff));
end

if isfield(G, 'elem') 
    [indx] = find(G.elem(5,:) == 3); % grey matter 
    elem = G.elem(1:4, indx); 
    node_g = sort(unique(elem(:))); % indx_node_grey
else
    node_g = 1:size(G.xyz, 2); 
end 

ng = length(node_g); 
node_h = []; 
A = []; 

for i = 1:nh
    diff = G.xyz(:, node_g) - xyzh(:,i) * ones(1, ng); 
    diff = sqrt(sum(diff.^2)); 
    
    if r == 0,
        A.node_h(i, 1) = node_g(find(diff == min(diff)));
    else
        indx_r = find(diff < r); 
        [dist_h, indx_h] = sort(diff(indx_r));
        node_h = node_g(indx_r(indx_h));
        
        A(i).dist_h = dist_h;
        A(i).node_h = node_h;
    end
end

% Calculate sensitivity matrix using Green's function
%--------------------------------------------------------------------------
nwav = size(G.s, 3); % number of wavelengths 
nroi = length(rois); 

for i = 1:nwav
    niter = size(A, 2);
    for j = 1:niter 
        S = []; 
        for k = 1:nroi,
            indx_s = find(R.s.label == R.ch_sd(rois(k), 2));
            indx_d = find(R.d.label == R.ch_sd(rois(k), 3));
            S(k,:) = (G.s(indx_s, A(j).node_h,i) .* G.d(indx_d, A(j).node_h,i))./(G.s(indx_s, node_d(indx_d),i)); 
        end
        A(j).(sprintf('S%d',i)) = S; 
    end
end
