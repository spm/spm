function normal = spm_eeg_inv_normals(vert,face)
% Compute the normals of a mesh
% FORMAT normal = spm_eeg_inv_normals(vert,face)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_normals.m 2696 2009-02-05 20:29:48Z guillaume $


m = struct('Vertices',vert,'Faces',face);
h = figure('Visible','off');
n = get(patch(m),'VertexNormals');
close(h);

f = sqrt(sum(n.^2,2));
I = find(f == 0);
for i = 1:length(I)
    n(I(i)) = n(I(i) - 1);
end
f           = sqrt(sum(n.^2,2));
normal(:,1) = n(:,1)./f;
normal(:,2) = n(:,2)./f;
normal(:,3) = n(:,3)./f;

return
%==========================================================================
