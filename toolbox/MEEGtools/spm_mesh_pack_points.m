function [Pout,ms2s ,ims2s,n] = spm_mesh_pack_points(S)
% Place approximately equally spaced points over a convex (ideally) mesh
% FORMAT [Pout,ms2s ,ims2s,n] = spm_mesh_pack_points(S)
%   S          - input structure
% Fields of S:
%   S.g        - gifti mesh                  - Default: mni scalp template
%   S.niter    - number of iterations        - Default: 2000
%   S.p        - initial points (nx3 matrix) - Default: guesses...
%   S.space    - desired spacing (mm)        - Default: 10
%   S.division - number of mesh subdivisions - Default: 3
%   S.nDens    - number of density checks    - Default: 40
% Output:
%   Pnew        - N x 3 matrix containing new points
%   ms2s        - nearest neighbour distances
%   ims2s       - initial nearest neighbour distances
%   n           - number of sensors at each iteration
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_mesh_pack_points.m 7740 2019-12-02 14:06:41Z guillaume $


%-Set default values and check arguments
%--------------------------------------------------------------------------
spm('FnBanner', mfilename);

if ~nargin, S = struct(); end

% mesh check
if ~isfield(S, 'g')
    mniScalp = fullfile(spm('Dir'),'canonical','scalp_2562.surf.gii');
    S.g   = gifti(mniScalp);
else
    if ~isa(S.g,'gifti')
        error('S.g should be of class gifti.');
    end
end

% space check
if ~isfield(S, 'space')
    S.space  = 10;
else
    if ~isnumeric(S.space)
        error('S.space should be numeric.');
    end
end

% niter check
if ~isfield(S, 'niter')
    S.niter  = 2000;
else
    if ~isnumeric(S.niter)
        error('S.niter should be numeric.');
    end
end

% point check
if ~isfield(S, 'p')
    disp('Setting up initial points.')
    S.p = poiGridSurface(S.g,S.space,S.space*.75);
else
    numer = isnumeric(S.p);
    if ~numer
        error('S.p should be numeric with 3 columns.');
    end
end

% divisions check
if ~isfield(S, 'division')
    S.division  = 3;
else
    if ~isnumeric(S.division)
        error('S.division should be numeric.');
    end
end

% divisions check
if ~isfield(S, 'nDens')
    S.nDens  = 40;
else
    if ~isnumeric(S.nDens)
        error('S.division should be numeric.');
    end
end

%- Subdivide surface and compute neighbourhood
%--------------------------------------------------------------------------
v = S.g.vertices;
f = S.g.faces;

for i=1:S.division
    [v,f]      = linearSubdivision(v,f);
    g          = gifti(struct('vertices',v,'faces',f));
    [dummy,Di] = spm_mesh_neighbours(g,1);
    medVertexSpacing = median(Di(:));
end

msg = 'Median point separation of solution space  is %3.3f mm.\n';
fprintf(msg,medVertexSpacing);

%- Assign sensors to nearest vertex
%--------------------------------------------------------------------------
% 4th column will hold index of sensors on surface
Pnew = zeros(length(S.p),4);
nSensors= length(Pnew);

for i = 1:nSensors
    tmp = S.p(i,:);
    di = sum(bsxfun(@minus,v,tmp).^2,2);
    [dummy,miInd] = min(di);
    Pnew(i,1:3)= v(miInd,:);
    Pnew(i,4)= miInd;
end

%- Compute initial sensor to sensor distance matrix
%--------------------------------------------------------------------------
disp('Computing initial point to point distance matrix.');

s2s = eye(size(Pnew,1));
for i = 1:size(Pnew,1)
    temp = Pnew(i,1:3);
    s2s(:,i) = sqrt(sum(bsxfun(@minus,Pnew(:,1:3),temp).^2,2));
    s2s(i,i) = Inf;
end
ims2s = min(s2s);

%- optimiisation
%--------------------------------------------------------------------------
disp('Optimising distance matrix.');
Pnew = optim_sens_dist(g,S.space,Pnew,10000,1:nSensors);

%- optimise sensor density
%--------------------------------------------------------------------------
disp('Optimising point density.');
n = zeros(S.nDens+1,1);
n(1) = size(Pnew,1);
for i = 1:S.nDens
    Pnew = increase_sens_dens(g,Pnew,S.space,S.space*.85);
    Pnew = optim_sens_dist(g,Inf,Pnew,S.niter,1:size(Pnew,1));
    Pnew = optim_sens_dist(g,S.space,Pnew,S.niter,1:size(Pnew,1));
    n(i+1)=size(Pnew,1);
    
    s2s = eye(size(Pnew,1));
    for j = 1:size(Pnew,1)
        temp=Pnew(j,1:3);
        s2s(:,j) = sqrt(sum(bsxfun(@minus,Pnew(:,1:3),temp).^2,2));
        s2s(j,j) = Inf;
    end
    msg = '%4d points with %3.3f mm spacing.\n';
    fprintf(msg,n(i+1),median(min(s2s)));
end

%- Compute penultimate sensor to sensor distance matrix
%--------------------------------------------------------------------------
msg = ['optimising points with nearest neighbour >'...
    num2str(S.space+.3) 'mm away.'];
disp(msg)

s2s = eye(size(Pnew,1));
for i = 1:size(Pnew,1)
    temp = Pnew(i,1:3);
    s2s(:,i) = sqrt(sum(bsxfun(@minus,Pnew(:,1:3),temp).^2,2));
    s2s(i,i) = Inf;
end
ms2s = min(s2s);

tooFar = find(ms2s>(S.space+.3),1);
if ~isempty(tooFar)
    Pnew = optim_sens_dist(g,S.space,Pnew,10000,1:size(Pnew,1));
end

%- Compute final sensor to sensor distance matrix
%--------------------------------------------------------------------------
s2s = eye(size(Pnew,1));
for i = 1:size(Pnew,1)
    temp = Pnew(i,1:3);
    s2s(:,i) = sqrt(sum(bsxfun(@minus,Pnew(:,1:3),temp).^2,2));
    s2s(i,i) = Inf;
end
ms2s = min(s2s);

msg = 'Final median sensor nearest neighbour is %3.3f mm.\n';
fprintf(msg,median(ms2s));
msg = 'Final sd of nearest neighbour is %3.3f mm.\n';
fprintf(msg,std(ms2s));

Pout = Pnew(:,1:3);
fprintf('%-40s: %30s\n','Completed',spm('time'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Linear Mesh Subdivision                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v1, f1] =  linearSubdivision(xyz, faces)

numverts = size(xyz,1);
numfaces = size(faces,1);
fk1 = faces(:,1);
fk2 = faces(:,2);
fk3 = faces(:,3);
% average the vertices
m1x = (xyz( fk1,1) + xyz( fk2,1) )/2;
m1y = (xyz( fk1,2) + xyz( fk2,2) )/2;
m1z = (xyz( fk1,3) + xyz( fk2,3) )/2;

m2x = (xyz( fk2,1) + xyz( fk3,1) )/2;
m2y = (xyz( fk2,2) + xyz( fk3,2) )/2;
m2z = (xyz( fk2,3) + xyz( fk3,3) )/2;

m3x = (xyz( fk3,1) + xyz( fk1,1) )/2;
m3y = (xyz( fk3,2) + xyz( fk1,2) )/2;
m3z = (xyz( fk3,3) + xyz( fk1,3) )/2;

vnew = [ [m1x m1y m1z]; [m2x m2y m2z]; [m3x m3y m3z] ];
[uvnew, dummy, jj] = unique(vnew, 'rows' );

m1 = jj(1:numfaces)+numverts;
m2 = jj(numfaces+1:2*numfaces)+numverts;
m3 = jj(2*numfaces+1:3*numfaces)+numverts;
tri1 = [fk1 m1 m3];
tri2 = [fk2 m2 m1];
tri3 = [m1 m2 m3];
tri4 = [m2 fk3 m3];

v1 = [xyz; uvnew]; % the new vertices
f1 = [tri1; tri2; tri3; tri4]; % the new faces


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     sensor distance optimisation                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sens = optim_sens_dist(surface,space,sens,niter,vec)
v=surface.vertices;
N=spm_mesh_neighbours(surface,1);
minmax=length(vec);
for i = 1:niter
    
    % randomly select 1 sensor
    worst = vec(randi(minmax));
    
    % sensor index on surface
    worstSurf = sens(worst,4);
    
    % Neighbour indices on surface
    worstNeighbours = N(worstSurf,:);
    
    % include current vertex as potential solution
    worstNeighbours=[nonzeros(worstNeighbours); worstSurf];
    
    % Get vertex coordinates
    neighP = v(worstNeighbours,:);
    
    % compute distance between current sensors and potential new sensors
    hmm = zeros(length(sens),length(worstNeighbours));
    for j =1:length(worstNeighbours)
        hmm(:,j)= sqrt(sum(bsxfun(@minus,sens(:,1:3),neighP(j,:)).^2,2));
    end
    
    %remove self distance
    n2s =hmm;
    n2s(worst,:)=[];
    
    % Pick neighbour closest to desired spacing or maximise
    if isinf(space)
        [dummy, miLoc] =  max(min(n2s));
        newPind = worstNeighbours(miLoc);
    else
        [dummy, miLoc] =  min(abs(min(n2s-space)));
        newPind = worstNeighbours(miLoc);
    end
    
    % check if invading a vertex ocupied by  a sensor already
    alreadyP = any(sens(:,4)==newPind);
    
    % update sensors
    if ~alreadyP
        sens(worst,1:3) = neighP(miLoc,:);
        sens(worst,4) = newPind;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Increase Density                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Pnew = increase_sens_dens(surface,Pnew,space,thresh)
%- Compute vertex to sensor matrix
%--------------------------------------------------------------------------
Pout = Pnew;
v=surface.vertices;
thresh=(thresh).^2;
P = Pout(:,1:3);
for i =1:size(v,1)
    v2s= sum((bsxfun(@minus,P,v(i,:))).^2,2);
    
    mi= min(v2s);
    if mi>thresh
        tmp=v(i,:);
        
        Pout(end+1,:)=[tmp,i];
        P(end+1,:)=tmp;
    end
end
Pnew = Pout;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Grid Surface Intersection                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outs = poiGridSurface(surface,space,lowerThresh)
% Find where a grid intersects a surface (ideally convex)

% Create Grid
%--------------------------------------------------------------------------
bbMin = min(surface.vertices);
bbMax = max(surface.vertices);
x = (bbMin(1)-space):space:(bbMax(1)+space);
y = (bbMin(2)-space):space:(bbMax(2)+space);
z = (bbMin(3)-space):space:(bbMax(3)+space);
[X,Y,Z] = ndgrid(x,y,z);
grid   = [X(:) Y(:) Z(:)];

% Keep only outer layer of grid
%--------------------------------------------------------------------------
mi = min(grid);
ma = max(grid);

miX = grid(find(abs(grid(:,1) - mi(1))<=.01),:);
miY = grid(find(abs(grid(:,2) - mi(2))<=.01),:);
miZ = grid(find(abs(grid(:,3) - mi(3))<=.01),:);

maX = grid(find(abs(grid(:,1) - ma(1))<=.01),:);
maY = grid(find(abs(grid(:,2) - ma(2))<=.01),:);
maZ = grid(find(abs(grid(:,3) - ma(3))<=.01),:);

% ny
%--------------------------------------------------------------------------
Pa = maY;
Pb  = bsxfun(@plus,maY,200.*[0,-1,0]);
outsny = ray_mesh_intersect2(Pa,Pb,surface);

% py
%--------------------------------------------------------------------------
Pa = miY;
Pb  = bsxfun(@plus,miY,200.*[0,1,0]);
outspy = ray_mesh_intersect2(Pa,Pb,surface);

% nx
%--------------------------------------------------------------------------
Pa = maX;
Pb  = bsxfun(@plus,maX,200.*[-1,0,0]);
outsnx = ray_mesh_intersect2(Pa,Pb,surface);
% px
%--------------------------------------------------------------------------
Pa = miX;
Pb  = bsxfun(@plus,miX,200.*[1,0,0]);
outspx = ray_mesh_intersect2(Pa,Pb,surface);

% nz
%--------------------------------------------------------------------------
Pa = maZ;
Pb  = bsxfun(@plus,maZ,200.*[0,0,-1]);
outsnz = ray_mesh_intersect2(Pa,Pb,surface);
% pz
%--------------------------------------------------------------------------
Pa = miZ;
Pb  = bsxfun(@plus,miZ,200.*[0,0,1]);
outspz = ray_mesh_intersect2(Pa,Pb,surface);


% Constrained adding of points
%--------------------------------------------------------------------------
uspz = unique(round(outsnz,3),'rows');

other = [unique(round(outsnx,3),'rows');...
    unique(round(outspy,3),'rows');...
    unique(round(outsny,3),'rows');...
    unique(round(outspz,3),'rows');...
    unique(round(outspx,3),'rows');...
    ];
us = uspz;

for i = 1:length(other)
    ot = other(i,:);
    di = sqrt(sum((bsxfun(@minus,us,ot)).^2,2));
    
    if(all(di>lowerThresh))
        us = [us;ot];
    end
    
end

% return
%--------------------------------------------------------------------------
outs = us;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     ray mesh intersection                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [poi] = ray_mesh_intersect2(Pa,Pb,surface)

% out put vector containing points
outs = zeros(length(Pb),3);

% gifti is slow to access multiple times in loops
ve = surface.vertices;
f  = surface.faces;

% preallocate some things
verts = double(zeros(3,3));
Y = zeros(3,1);
m = zeros(3,3);

% compute face centre coordinate
posi = zeros(size(f,1),3);
for i = 1:length(posi)
    whichVerts = f(i,:);
    cos = ve(whichVerts,:);
    posi(i,:) = mean(cos);
end


ninter=ones(size(Pa,1),1);
for j = 1:length(Pa)
    
    % initial point, final point and direction vector
    Ia   = Pa(j,:);
    Ib   = Pb(j,:);
    IaIb = [Ia;Ib];
    dir  = Ia-Ib;
    
    % define a box containing all possible solution faces
    up   = max(IaIb)+10;
    down = min(IaIb)-10;
    
    % check what faces are within box
    fup        = all(bsxfun(@lt,posi,up),2);
    fdown      = all(bsxfun(@gt,posi,down),2);
    smallFaces = f(fup&fdown,:);
    nfaces     = size(smallFaces,1);
    t          = ones(5,1)*Inf;
    
    % iterate over possible solutions
    for i = 1:nfaces
        % and  finds possible solusions
        % see https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/moller-trumbore-ray-triangle-intersection
        verts = ve(smallFaces(i,:),:);
        v0    = verts(1,:);
        v0v1  = verts(2,:)-v0;
        v0v2  = verts(3,:)-v0;
        pvec  = cross_product(dir,v0v2);
        deter = pvec*v0v1';
        
        if deter>0
            continue;
        end
        
        ideter = 1/deter;
        tvec   = (Ia-v0)';
        u      = pvec*tvec*ideter;
        
        if u<0 || u>1
            continue;
        end
        
        qvec   = cross_product(tvec,v0v1);
        v      = qvec * dir' * ideter;
        
        if v<0 || (u+v)>1
            continue;
        end
        t(ninter(j)) = v0v2 * qvec' * ideter;
        ninter(j) = ninter(j) + 1;
    end
    
    % algorithm always gives result but checks ensure result is valid
    [mul] = min(abs(t));
    
    % update output vector
    if(isfinite(mul))
        outs(j,:) = Ia-dir*mul;
    end
    
end

% subset for points that actually intersect the mesh
ind = sum(abs(outs),2)>0;
poi = outs(ind,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Cross product                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cp = cross_product(V1,V2)

cp = [(V1(2)*V2(3) - V2(2)*V1(3)),...
    (V1(3)*V2(1) - V2(3)*V1(1)),...
    (V1(1)*V2(2) - V2(1)*V1(2))];
