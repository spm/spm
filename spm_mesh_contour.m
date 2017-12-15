function S = spm_mesh_contour(M,z)
% Compute contour lines of a triangular mesh
% FORMAT S = spm_mesh_contour(M,z)
% M - a GIfTI object or patch structure
% z - height of z-plane
%
% S - structure of contour levels
%__________________________________________________________________________
%
% figure, hold on, axis equal
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
% z = linspace(min(M.vertices(:,3)),max(M.vertices(:,3)),20);
% for i=1:numel(z)
%   S = spm_mesh_contour(M,z(i));
%   for i=1:numel(S)
%     plot3(S(i).xdata,S(i).ydata,repmat(S(i).level,1,numel(S(i).xdata)));
%   end
% end
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_contour.m 7239 2017-12-15 17:14:33Z guillaume $


%-Check inputs and initialise output
%--------------------------------------------------------------------------
if isa(M,'gifti')
    M = export(M,'patch');
end
if ~all(isfield(M,{'vertices','faces'}))
    error('First input has to be a patch structure.');
end
if isinteger(M.faces), M.faces = double(M.faces); end

S = struct('xdata',{},'ydata',{},'level',{},'numel',{},'isopen',{});

%-Only consider triangles intersecting the z-plane
%--------------------------------------------------------------------------
X = M.vertices(:,1);
Y = M.vertices(:,2);
Z = M.vertices(:,3);
I = Z(M.faces) > z;
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
            if Ei == 1,     f = [2 3];
            elseif Ei == 2, f = [1 3];
            else            f = [1 2];
            end
            ii = find(sum(M.faces == M.faces(i,f(1)) | ...
                          M.faces == M.faces(i,f(2)),2)==2);
            ii(ii==i) = []; i = ii;
            if isempty(i), i = NaN; end
        end

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
    a  = (z-ze(:,1))./(ze(:,2)-ze(:,1));
    xc = a.*(xe(:,2)-xe(:,1))+xe(:,1);
    yc = a.*(ye(:,2)-ye(:,1))+ye(:,1);
    S(end+1) = struct('xdata',xc','ydata',yc','level',z,'numel',numel(xc),'isopen',isopen);
end
