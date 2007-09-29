function    [Center,Radius] = spm_eeg_inv_BestFitSph(points,varargin);

%=======================================================================
% FORMAT [Center,Radius] = spm_eeg_inv_BestFitSph(points,face)
% 
% This function computes the Center location and Radius (in mm) of the Sphere
% that fits best the 3D surface defined by the Nx3 coordinate matrix
% points
%
% optional: plot the surface and the fitted sphere (Gflag == 1)
%
% Uses a Least-Square approach:
%   - first compute an initial estimate of the radius and center
%   - then iterate a Gauss-Newton scheme to optimize the two values
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_BestFitSph.m 932 2007-09-29 18:01:16Z karl $

if nargin > 2
    disp('Wrong input format');
    return
elseif nargin == 2
    Gflag = varargin{1};
else
    Gflag = 1;
end

TOL = 1e8*eps;
Sv  = 40; % #vertices in sphere plot

N = size(points,1);
X = points(:,1);
Y = points(:,2);
Z = points(:,3);


% Initialisation (cost function: F1 = Ri^2 - R^2)
A  = [-2*X -2*Y -2*Z ones(N,1)];
B  = X.^2 + Y.^2 + Z.^2;
Op = inv(A'*A)*A';
P  = Op*B;

C0 = P(1:3);                    % initial center
R0 = sqrt( abs( sum(C0.^2) - P(4) ) ); % initial radius


% Gauss-Newton solution (cost function: F2 = Ri - R)
C  = C0;
R  = R0;
g  = TOL + 1;
k  = 0;

if Gflag == 1
    Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
    plot3(X,Y,Z,'ok');
    set(gcf,'color','white'); axis off; view(-90,0);
    drawnow
    hold on;
    tsph = spm_eeg_inv_TesSph(R,Sv);
    Nv = length(tsph.vert);
    for i = 1:3
        tsph.vert(i,:) = tsph.vert(i,:) + C(i)*ones(1,Nv);
    end
    ph = patch('Vertices',tsph.vert','Faces',tsph.tri,'EdgeColor','r','FaceColor','none');
    axis image
    drawnow
end

while (g > TOL)    
    Di = [ (X - C(1)*ones(N,1)) (Y - C(2)*ones(N,1)) (Z - C(3)*ones(N,1))];
    Ri = sqrt( sum( (Di.^2)' ) )';
    J  = - [ Di(:,1)./Ri Di(:,2)./Ri Di(:,3)./Ri ones(N,1) ];
    D  = Ri - R*ones(N,1);

    Op = inv(J'*J)*J';
    P  = - Op*D;

    C  = C + P(1:3);
    R  = R + P(4);
    
    g  = norm(J'*D);
    
    k  = k + 1;
    
    if Gflag == 1
        tsph = spm_eeg_inv_TesSph(R,Sv);
        Nv = length(tsph.vert);
        for i = 1:3
            tsph.vert(i,:) = tsph.vert(i,:) + C(i)*ones(1,Nv);
        end
        set(ph,'Vertices',tsph.vert','Faces',tsph.tri);
        drawnow
    end
end

Center = C;
Radius = R;

disp(['Scalp best fitting sphere computed (in ' num2str(k) ' iterations)']);
disp(sprintf('Centre = [%5.4f %5.4f %5.4f] (Radius = %5.4f)',C(1),C(2),C(3),R));
