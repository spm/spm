function spm_mountaincar_movie(DEM)
% makes a move for mountain car problem
% FORMAT spm_mountaincar_movie(DEM)
%
% see:
% Gaussian Processes in Reinforcement Learning
% Carl Edward Rasmussen and Malte Kuss
% Max Planck Institute for Biological Cybernetics
% Spemannstraﬂe 38, 72076 T®ubingen, Germany
% {carl,malte.kuss}@tuebingen.mpg.de
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_mountaincar_movie.m 2033 2008-09-02 18:32:14Z karl $
 
 
% get postion
%==========================================================================
p     =  DEM.qU.v{1}(1,:);
 
% get hieght: H(x)
%--------------------------------------------------------------------------
dx    = 1/1024;
x     = linspace(-2,2,1/dx);
xx    = x.^2;
dHdx  = (x < 0).*(2*x + 1);
dHdx  = (x > 0).*(1./(1 + 5*xx).^(1/2) - 5*xx./(1 + 5*xx).^(3/2) + (x/2).^4) + dHdx;
H     = cumsum(dHdx)*dx;
H     = H - min(H);

% get image
%--------------------------------------------------------------------------
car   = imread('cable.jpg');
[h j] = min(abs(x - 1));
X0    = .4;
Y0    = .5;
H0    = H(j);
for i = 1:length(p)
    
    % postiotn and scling
    %----------------------------------------------------------------------
    X     = p(i);
    S     = exp(-X/4);
    SX    = S*200;
    SY    = S*1000;
    
    [h j] = min(abs(x - X));
    Y     = H(j);
    
    % draw image and action
    %----------------------------------------------------------------------
    image(car); hold on
    plot((1 + X0 - X)*SX,(Y - H0)*SY,'g.','MarkerSize',24)
    plot((x + X0 - X)*SX,(Y - H)*SY,'k')
    if DEM.qU.a{2}(i) > 0
        plot((-1.8 + X0 - X)*SX,Y*SY,'r.','MarkerSize',32)
    end
    axis(([-2 2 0 .7] + [(X0*[1 1] - X) (Y - Y0*[1 1])]).*[SX SX SY SY])
    hold off
    axis square off
    drawnow
    
    % save
    %----------------------------------------------------------------------
    M(i) = getframe(gca,[1 1 239 240]);
    
end

% set ButtonDownFcn
%--------------------------------------------------------------------------
h = findobj(gca,'type','image');
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

