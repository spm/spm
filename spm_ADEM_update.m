function [DEM] = spm_ADEM_update(DEM)
% Updates ADEM structure using conditional expectations
% FORMAT [DEM] = spm_ADEM_update(DEM)
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ADEM_update.m 6039 2014-06-04 18:50:28Z karl $
 
% update states and action
%--------------------------------------------------------------------------
n     = length(DEM.M);
for i = 1:(n - 1)
    DEM.M(i).x  = spm_unvec(DEM.qU.x{i}(:,end),DEM.M(i).x);
    DEM.M(i).pE = DEM.qP.P{i};
end
for i = 1:n
    if ~isempty(DEM.M(i).v)
        DEM.M(i).v  = spm_unvec(DEM.qU.v{i}(:,end),DEM.M(i).v);
    end
end

n     = length(DEM.G);
for i = 1:(n - 1)
    DEM.G(i).x  = spm_unvec(DEM.pU.x{i}(:,end),DEM.G(i).x);
end
for i = 1:n
    if ~isempty(DEM.G(i).v)
        DEM.G(i).v  = spm_unvec(DEM.pU.v{i}(:,end),DEM.G(i).v);
    end
end
DEM.G(n).a      = spm_unvec(DEM.qU.a{n}(:,end),DEM.G(n).a);
