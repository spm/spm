function T = spm_mesh_inside(M,XYZ)
% Test whether a point is inside or outside a watertight triangle mesh
% FORMAT T = spm_mesh_inside(M,XYZ)
% M        - a patch structure or GIfTI object 
% XYZ      - a 1 x 3 vector of point coordinates {mm}
%
% T        - logical scalar indicating inside/outside mesh test
%__________________________________________________________________________
%
% Compute generalised winding number from equations 5 and 6 from:
% "Robust Inside-Outside Segmentation using Generalized Winding Numbers"
% Alec Jacobson, Ladislav Kavan, Olga Sorkine-Hornung, ACM SIGGRAPH 2013.
% https://igl.ethz.ch/projects/winding-number/
%__________________________________________________________________________
%
% M = gifti('mesh.gii');
% M = export(M,'patch');
%
% m = max(M.vertices,[],1);
% n = min(M.vertices,[],1);
% P = (m-n).*rand(4096,3) + n;
%
% for i=1:size(P,1)
%     T(i) = spm_mesh_inside(M,P(i,:));
% end
%
% figure, plot3(P(T,1), P(T,2), P(T,3), '.')
% H = spm_mesh_render(M);
% hold(H.axis,'on');
% plot3(P(T,1), P(T,2), P(T,3), '.','Parent',H.axis)
% plot3(P(~T,1), P(~T,2), P(~T,3), '.r','Parent',H.axis)
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


V = M.vertices - XYZ;
F = M.faces;

%-Vertices for all faces
V = V(F',:);

%-Lengths [3 x nf]
L = squeeze(sqrt(sum(V.^2,2)));
L = reshape(L,3,[]);

%-Vertices [3 x 3 x nf]
V = permute(reshape(V',3,3,[]),[2 1 3]);

%-Winding number (solid angle)
w = 0;
for i=1:size(F,1)
    A = V(1,:,i)';
    B = V(2,:,i)';
    C = V(3,:,i)';
    a = L(1,i);
    b = L(2,i);
    c = L(3,i);

    w = w + atan2(det([A B C]), a*b*c + A'*B*c + B'*C*a + C'*A*b);
end

T = w >= 2*pi;
