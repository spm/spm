function C = spm_mesh_isoline(M, T, t)
% Compute isolines on a triangular mesh
% FORMAT C = spm_mesh_isoline(M, T, t)
% M   - a GIfTI object or patch structure
% T   - [vx1] data vector
% t   - isovalue [Default: 0]
%
% C   - struct array of isolines with fields 'xdata', 'ydata', 'zdata' and
%       'isopen'
%__________________________________________________________________________
%
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
% M = export(M,'patch');
% M = spm_mesh_inflate(M);
% T = randn(size(M.vertices,1),1);
% T = spm_mesh_smooth(M,T,100);
% H = spm_mesh_render('Disp',M);
% H = spm_mesh_render('Overlay',H,T);
% hold on
% t = linspace(min(T),max(T),20);
% for i=1:numel(t)
%   C = spm_mesh_isoline(M,T,t(i));
%   for j=1:numel(C)
%     plot3(C(j).xdata,C(j).ydata,C(j).zdata,'k-');
%   end
% end
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_isoline.m 7618 2019-06-17 12:29:46Z guillaume $


if nargin < 3, t = 0; end

C = spm_mesh_contour(M,struct('T',T,'t',t));
