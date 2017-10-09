function S = spm_mesh_contour(M,z)
% Compute contour lines of a triangular mesh
% FORMAT S = spm_mesh_contour(M,z)
% M - a GIfTI object or patch structure
% z - height of z-plane
%
% S - structure of contour levels
%__________________________________________________________________________
%
% figure
% hold on
% M = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
% S = spm_mesh_contour(M,linspace(min(M.vertices(:,3)),max(M.vertices(:,3)),20));
% for i=1:numel(S)
%     plot3(S(i).xdata,S(i).ydata,repmat(S(i).level,1,numel(S(i).xdata)));
% end
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_contour.m 7183 2017-10-09 15:26:47Z guillaume $


if numel(z) == 1, z = [z z]; end
[C,H] = tricontour(M.faces, ...
    M.vertices(:,1), M.vertices(:,2), M.vertices(:,3), z);
if ~isempty(C)
    S = contourdata(C);
else
    S = struct('xdata',{},'ydata',{},'level',{},'numel',{},'isopen',{});
end


%==========================================================================
% Copyright (c) 2006, Duane Hanselman
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==========================================================================

function [c,h]=tricontour(tri,x,y,z,nv)
%TRICONTOUR Triangular Contour Plot.
% TRICONTOUR(TRI,X,Y,Z,N) draws scalar N contour lines treating the values
% in Z as heights above a plane. TRI,X,Y,and Z define a triangulation where
% the triangles are defined by the M-by-3 face matrix TRI, such as that
% returned by DELAUNAY. Each row of TRI contains indices into the X,Y, and
% Z vertex vectors to define a single triangular face. Contours are
% computed directly from the triangulation rather than interpolating back
% to a cartesian grid using GRIDDATA.
% TRICONTOUR(TRI,X,Y,Z,V) draws length(V) contour lines at the values
% specified in vector V.
% TRICONTOUR(TRI,X,Y,Z,[v v]) draws a single contour line at the level v.
%
% [C,H] = TRICONTOUR(...) returns contour matrix C as described in CONTOURC
% and a vector of handles H to the created patch objects.
% H can be used to set patch properties.
% CLABEL(C) or CLABEL(C,H) labels the contour levels.
%
% Example:
%           x=linspace(-3,3,39);
%           y=linspace(-2.5,2.5,49);
%           [xx,yy]=meshgrid(x,y);
%           zz=peaks(xx,yy);
%           v=-3:1:5; % contour levels
%           subplot(1,2,1)
%           [C,h]=contour(xx,yy,zz,v);   % standard contour for comparison
%           clabel(C)
%           title Contour
%
%           idx=randperm(numel(zz));     % grab some scattered indices
%           n=idx(1:ceil(numel(zz)/2))'; % one half of them
%           x=xx(n);                     % get scattered data
%           y=yy(n);
%           z=zz(n);
%           tri=delaunay(x,y);           % triangulate scattered data
%           subplot(1,2,2)
%           [C,h]=tricontour(tri,x,y,z,v);
%           clabel(C,h)
%           title TriContour
%
% view(3) displays the contour in 3-D.
%
% See also DELAUNAY, CONTOUR, TRIMESH, TRISURF, TRIPLOT, PATCH.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-05-07, 2006-05-16, 2006-07-25

% https://www.mathworks.com/matlabcentral/fileexchange/38858

if nargin<5
    error('Not Enough Input Arguments.')
end
x=x(:);	% convert input data into column vectors
y=y(:);
z=z(:);
xlen=length(x);
if ~isequal(xlen,length(y),length(z))
    error('X, Y, and Z Must Have the Same Number of Elements.')
end
if size(tri,2)~=3 || any(tri(:)<0) || any(tri(:)>xlen)
    error('TRI Must Be a Valid Triangulation of the Data in X, Y, Z.')
end

zs=z(tri);
zmax=max(max(zs));              % find max and min in z data that is in tri
zmin=min(min(zs));

if length(nv)==1                                 % nv is number of contours
    zlev=linspace(zmax,zmin,nv+2);
elseif length(nv)==2 && nv(1)==nv(2)              % nv is one contour level
    zlev=nv(1);
else                                       % nv is vector of contour levels
    zlev=sort(nv,'descend');
end
zlev(zlev>=zmax | zlev<=zmin)=[];  % eliminate contours outside data limits
nlev=length(zlev);

if nlev==0
    %warning('No Contours to Plot. Chosen Contours Outside Limits of Data.');
    c = []; h = []; return;
end

% precondition the input data
[zs,zidx]=sort(zs,2);         % sort vertices by z value ascending
%for k=1:size(zs,1)            % shuffle triangles to match
%    tri(k,:)=tri(k,zidx(k,:));
%end
tri2 = tri'; tri = tri2(zidx+repmat(3*(0:size(zidx,1)-1)',1,3)); clear tri2;

if nargout<2 % Added by Torben to suppress graphics
    hax=newplot;                  % create new axis if needed
end
h=[];                         % patch handle storage

C=zeros(2,0);                 % Clabel data storage
cs=[2 1];                     % column swap vector cs(1)=2, cs(2)=1;

% Main Loop ---------------------------------------------------------------
for v=1:nlev                  % one contour level at a time
    zc=zlev(v);                % chosen level
    above=zs>=zc;              % true for vertices above given contour
    numabove=sum(above,2);     % number of triangle vertices above contour
    tri1=tri(numabove==1,:);   % triangles with one vertex above contour
    tri2=tri(numabove==2,:);   % triangles with two vertices above contour
    n1=size(tri1,1);           % number with one vertex above
    n2=size(tri2,1);           % number with two vertices above
    
    edge=[tri1(:,[1 3])        % first column is indices below contour level
        tri1(:,[2 3])        % second column is indices above contour level
        tri2(:,[1 2])
        tri2(:,[1 3])];
    if n1==0                   % assign edges to triangle number
        n=[1:n2 1:n2]';
    elseif n2==0
        n=[1:n1 1:n1]';
    else
        n=[1:n1 1:n1 n1+(1:n2) n1+(1:n2)]';
    end
    
    [edge,idx]=sortrows(edge);    % put shared edges next to each other
    n=n(idx);                     % shuffle triangle numbers to match
    
    idx=all(diff(edge)==0,2);     % find shared edges
    idx=[idx;false]|[false;idx];  % True for all shared edges
    
    % eliminate redundant edges, two triangles per interior edge
    edgeh=edge(~idx,:);           % hull edges
    nh=n(~idx);                   % hull triangle numbers
    if ~isempty(nh)
        nh(end,2)=0;               % zero second column for hull edges
    end
    edges=edge(idx,:);            % shared edges
    edges=edges(1:2:end-1,:);     % take only unique edges
    ns=n(idx);                    % interior triangle numbers
    ns=[ns(1:2:end) ns(2:2:end)]; % second column is second triangle
    edge=[edgeh;edges];           % unique edges
    nn=[nh;ns];                   % two columns of triangle numbers
    ne=size(edge,1);              % number of edges
    
    flag=true(ne,2);              % true for each unused edge per triangle
    tmp=zeros(ne+1,1);            % contour data temporary storage
    
    xe=x(edge);                   % x values at vertices of edges
    ye=y(edge);                   % y values at  vertices of edges
    ze=z(edge);                   % z data at  vertices of edges
    
    alpha=(zc-ze(:,1))./(ze(:,2)-ze(:,1)); % interpolate all edges
    xc=alpha.*(xe(:,2)-xe(:,1)) + xe(:,1); % x values on this contour
    yc=alpha.*(ye(:,2)-ye(:,1)) + ye(:,1); % y values on this contour
    
    while any(flag)	% while there are still unused edges -----------------
        
        xtmp=tmp;
        ytmp=tmp;
        [ir,ic]=find(flag,1);            % find next unused edge
        flag(ir,ic)=false;               % mark this edge used
        
        k=1;                             % first data point in subcontour
        xtmp(k)=xc(ir);                  % store data from this edge
        ytmp(k)=yc(ir);
        
        while true     % complete this subcontour ---------------------------
            
            [ir,ic]=find(flag&nn(ir,ic)==nn,1);% find other edge of triangle
            flag(ir,ic)=false;            % mark this edge used
            k=k+1;
            xtmp(k)=xc(ir);               % store data from this edge
            ytmp(k)=yc(ir);
            
            ic=cs(ic);                    % other triangle that shares edge
            
            if nn(ir,ic)==0               % reached hull, subcontour complete
                k=k+1;
                xtmp(k)=nan;               % don't let subcontour close
                ytmp(k)=nan;
                break
            elseif ~flag(ir,ic)           % complete closed subcontour
                break
            else                          % more points remain on subcontour
                flag(ir,ic)=false;         % mark this edge used
            end
        end % while true ----------------------------------------------------
        xtmp(k+1:end)=[];                % throw away unused storage
        ytmp(k+1:end)=[];                % xtmp,ytmp contain subcontour
        
        if nargout<2                     % plot the subcontour
            patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
                'Parent',hax,'FaceColor','none','EdgeColor','flat',...
                'UserData',zc)
            C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
        else %MODIFIED by TEL to suppress plot % plot subcontour and create output
            %          h=[h;patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
            %          'Parent',hax,'FaceColor','none','EdgeColor','flat',...
            %          'UserData',zc)]; %#ok
            C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
        end
    end % while any(flag) --------------------------------------------------
end % for v=1:nlev
if nargout
    c=C;
end

function s = contourdata(c)
%CONTOURDATA Extract Contour Data from Contour Matrix C.
% CONTOUR, CONTOURF, CONTOUR3, and CONTOURC all produce a contour matrix
% C that is traditionally used by CLABEL for creating contour labels.
%
% S = CONTOURDATA(C) extracts the (x,y) data pairs describing each contour
% line and other data from the contour matrix C. The vector array structure
% S returned has the following fields:
%
% S(k).level contains the contour level height of the k-th line.
% S(k).numel contains the number of points describing the k-th line.
% S(k).isopen is True if the k-th contour is open and False if it is closed.
% S(k).xdata contains the x-axis data for the k-th line as a column vector.
% S(k).ydata contains the y-axis data for the k-th line as a column vector.
%
% For example: PLOT(S(k).xdata,S(k).ydata)) plots just the k-th contour.
%
% See also CONTOUR, CONTOURF, CONTOUR3, CONTOURC.

% From the help text of CONTOURC:
%   The contour matrix C is a two row matrix of contour lines. Each
%   contiguous drawing segment contains the value of the contour, 
%   the number of (x,y) drawing pairs, and the pairs themselves.  
%   The segments are appended end-to-end as
% 
%       C = [level1 x1 x2 x3 ... level2 x2 x2 x3 ...;
%            pairs1 y1 y2 y3 ... pairs2 y2 y2 y3 ...]

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2007-05-22

% https://www.mathworks.com/matlabcentral/fileexchange/38863

if nargin<1 || ~isfloat(c) || size(c,1)~=2 || size(c,2)<4
   error('CONTOURDATA:rhs',...
         'Input Must be the 2-by-N Contour Matrix C.')
end

tol=1e-12;
k=1;     % contour line number
col=1;   % index of column containing contour level and number of points

while col<size(c,2); % while less than total columns in c
   s(k).level = c(1,col); %#ok
   s(k).numel = c(2,col); %#ok
   idx=col+1:col+c(2,col);
   s(k).xdata = c(1,idx); %#ok
   s(k).ydata = c(2,idx); %#ok
   s(k).isopen = abs(diff(c(1,idx([1 end]))))>tol || ...
                 abs(diff(c(2,idx([1 end]))))>tol; %#ok
   k=k+1;
   col=col+c(2,col)+1;
end
