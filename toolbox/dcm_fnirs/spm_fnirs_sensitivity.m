function A = spm_fnirs_sensitivity(G, pos_o, pos_h, r, roi_ch) 
% Calculate sensitivity matrix which corresponds to the effective
% pathlength of detected photons for the channel measurements in the
% hemodynamic source. 
% FORMAT [A] = spm_fnirs_sensitivity(G, pos_o, pos_h, r, roi_ch) 
% 
% G        - estimated Green's function 
% pos_o    - MNI coordinates of source and detector positions 
% pos_h    - MNI coordinates of hemodynamic source positions 
% r        - radius of spatially distributed sources 
% roi_ch   - channels of interets 
%
% A        - sensitivity matrix 
% 
% Green's function
%--------------------------------------------------------------------------
% G.g      - estimated Green's function [# channels x # voxels x # wavelengths]
% G.xyz    - MNI locations [3 x # voxels] 
% G.tissue - tissue types of voxels [3 x # voxels] (optional)
%            1-scalp, 2-CSF, 3-gray matter, 4-white matter 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $id$

if ~exist('roi_ch'), roi_ch = 1:pos.ch.nch; end

if ~exist('r') r = 0; end 

ns = pos_o.s.ns; % # of source
nd = pos_o.d.nd; % # of detector
xyz_d = pos_o.d.xyzH; 

node_mm = G.xyz; 
n_node = size(G.xyz, 2); 

%-Find indices of Green's function of detector and hemodynamic source 
%--------------------------------------------------------------------------
node_d = zeros(nd, 1);
for i = 1:nd
    diff = node_mm - xyz_d(:, i) * ones(1, n_node);
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

nh = size(pos_h, 2); 
ng = length(node_g); 
node_h = []; 
A = []; 

for i = 1:nh
    diff = node_mm(:, node_g) - pos_h(:,i) * ones(1, ng); 
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
nwav = 2; 
nroi = length(roi_ch); 

for i = 1:nwav
    niter = size(A, 2);
    for j = 1:niter 
        S = []; 
        for k = 1:nroi,
            indx_s = find(pos_o.s.label == pos_o.ch_sd(roi_ch(k), 1));
            indx_d = find(pos_o.d.label == pos_o.ch_sd(roi_ch(k), 2));
            S(k,:) = (G.s(indx_s, A(j).node_h,i) .* G.d(indx_d, A(j).node_h,i))./(G.s(indx_s, node_d(indx_d),i)); 
        end
        A(j).(sprintf('S%d',i)) = S; 
    end
end
