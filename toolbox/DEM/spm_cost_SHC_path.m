function spm_cost_SHC_path(qU,A)
% plots path for cost_SHC demo's
% FORMAT spm_cost_SHC_path(qU,A)
%
% qU  - DEM condotioal esimates of states
% A.x - locations of attrcuor
% A.d - radius
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cost_SHC_path.m 3757 2010-03-08 11:41:53Z guillaume $

% plot
%==========================================================================
[X Y Z] = sphere;

% plot path
%--------------------------------------------------------------------------
x     = qU.v{1}(1:2,:);
plot(A.x(1,:),A.x(2,:),'c.','MarkerSize',64), hold on
plot(x(1,:),x(2,:));    
plot(x(1,:),x(2,:),'.','MarkerSize',8);     

% plot locations
%--------------------------------------------------------------------------
% for i = 1:size(A.x,2)
%     surf(X*A.d + A.x(1,i),Y*A.d + A.x(2,i),Z*A.d - 8)
% end
% shading interp


axis([-4 4 -4 4])
axis square
title('trajectory','FontSize',16),hold off

% occupancy
%--------------------------------------------------------------------------
for i = 1:size(A.x,2)
    for j = 1:2
        d(j,:) = x(j,:) - A.x(j,i);
    end
    b(i,:) = sqrt(sum(d.^2)) < 2*A.d;
end

% plot precent occupancy
%--------------------------------------------------------------------------
p     = 100*sum(b,2)/length(b);
for i = 1:length(p)
    text(A.x(1,i)*2,A.x(2,i)*1.4,sprintf('%2.0f%%',p(i)),'Fontsize',12)
end
hold off


