function vol = triangle4pt(vol)

% FORMAT vol = triangle4pt(vol)
%
% Takes the volume model and estimates the 4th point of each triangle of
% each mesh.
% 
% In each vol.bnd sub-structure, a field '.pnt4' is added. The '.pnt4'
% field is a Ntri*3 matrix, with the coordinates of a point for each
% triangle in the meshed surface.
%
% Explanations:
% The point is that for some BEM, specifically 'solid angle', calculation
% it is necessary to estimate the local curvature of the true surface which
% is approximated by the flat triangle. One way to proceed is to use 
% "close by" vertices to estimate the overall area's curvature. 
% A more elegant(?) way uses a 4th point for each triangle: the "centroid" 
% of the triangle is simply pusehd away from the triangle surface to fix 
% the local surface curvature (assuming the surface is smooth enough).
% This 4th point is thus hovering above/under the triangle and can be used 
% to fit a sphere on the triangle in a realistic way. 
% 
% Method:
% - The 4th point can/could be defined at the tessalation stage, based on
%   the anatomical images directly. 
% - With any model, the curvature can be estimated/approximated by looking
%   at the vertices around the triangle considered and fit a sphere on
%   those few vertices, assuming the surface is smooth enough
% The latter option is the one followed here.
% The extra-vertices considered here are those 3 which are linked to the
% triangle by 2 edges.
%__________________________________________________________________________
%
% written by Christophe Phillips, 2009/01/19
% Cyclotron Research Centre, University of liège, belgium

Ns = length(vol.bnd);
for ii=1:Ns % treat each mesh one at a time
    tri = vol.bnd(ii).tri;
    pnt = vol.bnd(ii).pnt;
    Nt = size(tri,1);
    pnt4 = zeros(Nt,3);
    for jj=1:Nt % treat each triangle on a t a time
        lt1 = find(tri(:,1)==tri(jj,1) | tri(:,2)==tri(jj,1) | ...
                                                    tri(:,3)==tri(jj,1));
        lt2 = find(tri(:,1)==tri(jj,2) | tri(:,2)==tri(jj,2) | ...
                                                    tri(:,3)==tri(jj,2));
        lt3 = find(tri(:,1)==tri(jj,3) | tri(:,2)==tri(jj,3) | ...
                                                    tri(:,3)==tri(jj,3));
        lt = [intersect(lt1,lt2); intersect(lt2,lt3); intersect(lt3,lt1)];
        lt(lt==jj) = [];
        % list of 3 directly surrounding triangles
        lv = tri(lt,:);
        lv = setxor(lv(:),tri(jj,:));
        % list of 3 voxels connected by 2 edges to the jj_th triangle.
        sph_pnt = pnt([tri(jj,:) lv],:);
        [center,radius] = sphfit(sph_pnt);
        % best fitting sphere radius & centre, for the 6 points chosen
        pnt_c = sum(pnt(tri(jj,:),:))/3;
        % centroid of the triangle treated
        tmp = pnt_c-center;
        vd = tmp/norm(tmp);
        % unit vector from sphere center in direction of middle triangle
        pnt4(jj,:) = center+vd*radius ;
        % projection of the centroid, at 'radius' from the center
    end
    vol.bnd(ii).pnt4 = pnt4;
end

return

% figure, plot3(pnt(:,1),pnt(:,2),pnt(:,3),'.')
% axis equal, hold on
% plot3(pnt(tri(jj,:),1),pnt(tri(jj,:),2),pnt(tri(jj,:),3),'*')
% plot3(pnt(lv,1),pnt(lv,2),pnt(lv,3),'o')
% plot3(pnt4(jj,1),pnt4(jj,2),pnt4(jj,3),'rs')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [center,radius]=sphfit(vc,Ni,stopth)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fits a sphere to a set of surface points
% usage: [center,radius]=sphfit(vc,Ni,stopth)
%
% input:
% vc        nx3 matrix, where each row represents the location
%           of one surface point. vc can have more than 3 columns
%           (e.g. orientations) - then only the first 3 columns are used
% Ni        number of iteration requested, 20 by default
% stopth    stopping threshold used for the relative difference between 2
%           succeessive estimates of the radius. Fixed 10^-6 by default
%           If stopth<0, then no stopping till the end of the Ni iterations
%
% center  1x3 vector denoting the center
% radius  scalar denoting the radius
%
% written by Guido Nolte
% updated by Christophe Phillips, 2009/1/19
%   - add the number of iterations as input, and use 20 as default
%       instead of 5
%   - add a stopping criteria based on the relative difference between 2
%       successive estimates of the radius. 
%       If rel_diff<1e-6 (default of use fixed), then break out.

if nargin<3, stopth = 1e-6; end
if nargin<2, Ni = 20; end
if isempty(Ni), Ni = 20; end

vc=vc(:,1:3);
[nvc,ndum]=size(vc);

center_0=mean(vc);
vcx=vc-repmat(center_0,nvc,1);
radius_0=mean(sqrt(vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2));

alpha=1;
rd = 1;
err_0=costfun(vc,center_0,radius_0);

for k=1:Ni;

    [center_new,radius_new]=lm1step(vc,center_0,radius_0,alpha);

    err_new=costfun(vc,center_new,radius_new);
%     disp([k,err_0,err_new,center_new,radius_new]);

    if err_new<err_0;
        rd = abs(radius_0-radius_new)/radius_0;
%         disp(rd)
        center_0=center_new;
        radius_0=radius_new;
        err_0=err_new;
        alpha=alpha/5;
    else
        alpha=alpha*5;
    end

    radius=radius_0;
    center=center_0;
    if rd<stopth, break, end % stop if

end

return;

%%
function err=costfun(vc,center,radius)
[nvc,ndum]=size(vc);

vcx=vc-repmat(center,nvc,1);

err=sqrt(sqrt(mean( (vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2).^2)));

return

%%
function  [center_new,radius_new]=lm1step(vc,center,radius,alpha)

[nvc,ndum]=size(vc);
vcx=vc-repmat(center,nvc,1);
f=vcx(:,1).^2+vcx(:,2).^2+vcx(:,3).^2-radius^2;

L=2*[vcx,repmat(radius,nvc,1)];

par_new=inv(L'*L+alpha*eye(4))*L'*f;

center_new=center+par_new(1:3)';
radius_new=radius+par_new(4);

return;



