function S = spm_mesh_contour(M,mat)
% Compute contour lines of a triangular mesh
% FORMAT S = spm_mesh_contour(M,z)
% M   - a GIfTI object or patch structure
% z   - height of z-plane
%
% FORMAT S = spm_mesh_contour(M,mat)
% mat - 4 x 4 transformation matrix
%       (use z-plane at z = 0 after linear transformation according to mat)
%
% S   - struct array of contour lines with fields 'xdata', 'ydata',
%       'zdata' and 'isopen'
%__________________________________________________________________________
%
% figure, hold on, axis equal
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
% z = linspace(min(M.vertices(:,3)),max(M.vertices(:,3)),20);
% for i=1:numel(z)
%   S = spm_mesh_contour(M,z(i));
%   for j=1:numel(S)
%     plot3(S(j).xdata,S(j).ydata,S(j).zdata);
%   end
% end
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_contour.m 7392 2018-08-13 11:50:28Z guillaume $


%-Input and output arguments
%==========================================================================

%-Check input M
%--------------------------------------------------------------------------
if isa(M,'gifti')
    M = export(M,'patch');
end
if ~all(isfield(M,{'vertices','faces'}))
    error('First input has to be a patch structure.');
end
if isinteger(M.faces), M.faces = double(M.faces); end

%-Check input mat
%--------------------------------------------------------------------------
if numel(mat) == 1
    mat = spm_matrix([0 0 -mat]);
end

%-Initialise output
%--------------------------------------------------------------------------
S = struct('xdata',{},'ydata',{},'zdata',{},'isopen',{});


%-Only consider triangles intersecting the z-plane at z = 0
%==========================================================================
M.vertices = (mat(1:3,:) * [M.vertices';ones(1,size(M.vertices,1))])';
X = M.vertices(:,1);
Y = M.vertices(:,2);
Z = M.vertices(:,3);
I = Z(M.faces) > 0;
J = sum(I,2);
J = J > 0 & J < 3;
M.faces = M.faces(J,:);


%-Marching squares (https://en.wikipedia.org/wiki/Marching_squares)
%==========================================================================
I = I(J,:) * [4 2 1]'; % binary index in base 10
F = true(size(I,1),1); % available triangles
E = [2 1 2 3 3 1]; % forward lookup table of the 6 possibilities

while ~isempty(F)
    i      = 1; % index of current triangle
    F(i)   = false;
    isopen = false;
    
    %-Initialise contour
    %----------------------------------------------------------------------
    j  = 0; % number of points in contour
    C  = zeros(2*size(M.faces,1),1); % contour indices
    
    %-Store all edges
    %----------------------------------------------------------------------
    ed = [M.faces(:,[2 3]);M.faces(:,[1 3]);M.faces(:,[1 2])];
    
    %-Follow contour
    %----------------------------------------------------------------------
    while true
        Ei = E(I(i));
        
        %-Store edge index
        %------------------------------------------------------------------
        j = j + 1; C(j) = i + (Ei-1)*size(M.faces,1);
        
        %-Find next triangle (see also triangulation.neighbors)
        %------------------------------------------------------------------
        try
            ii = spm_mesh_utils('neighbouringfaces',M.faces,i);
            i = ii(Ei);
        catch
            % non-MEX implementation
            f = [1 2 3]; f(Ei) = [];
            ii = find(sum(M.faces == M.faces(i,f(1)) | ...
                          M.faces == M.faces(i,f(2)),2)==2);
            ii(ii==i) = []; i = ii;
            if isempty(i), i = NaN; end
        end
        if ~isnan(i) && i ~= 1 && ~F(i), i = NaN; end

        %-Detect dead end or loop
        %------------------------------------------------------------------
        if isnan(i)
            if isopen, break; end
            % try to go from start backwards
            isopen = true;
            i = 1;
            E = fliplr(E);
            C(1:j) = C(j:-1:1);
        elseif i == 1
             % loop the loop
            j = j + 1; C(j) = C(1);
            break;
        end
        F(i) = false;
    end
    
    %-Discard used triangles
    %----------------------------------------------------------------------
    M.faces = M.faces(F,:);
    I = I(F);
    F = F(F);
    
    %-Linear interpolation of coordinates
    %----------------------------------------------------------------------
    ed = ed(C(C>0),:);
    xe = X(ed); ye = Y(ed); ze = Z(ed);
    a  = ze(:,1) ./ diff(ze,1,2);
    xc = xe(:,1) - a .* diff(xe,1,2);
    yc = ye(:,1) - a .* diff(ye,1,2);
    XYZ = mat\[xc,yc,zeros(size(xc)),ones(size(xc))]';
    S(end+1) = struct('xdata',XYZ(1,:),'ydata',XYZ(2,:),'zdata',XYZ(3,:),'isopen',isopen);
end
