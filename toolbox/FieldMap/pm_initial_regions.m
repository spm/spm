function [irima,cn] = pm_initial_regions(pm,mask,nstep)
%
% Divides 2 or 3D phasemap (pm) into nstep equally wide
% angle ranges and returns the connected components
% of those.
% FORMAT [irima,cn] = pm_initial_regions(pm,mask,nstep)
%
% Input
% pm      : Non-unwrapped phase-map.
% mask    : Tells us what regions of pm to consider.
% nstep   : Defines the number of equi-wide angle ranges
%           between -pi and pi that we should use.
%           If linear phase-ramps have been removed from
%           the data we may have values outside the 
%           -pi->pi range. We will then simply divide
%           the observed range into nstep steps.
%
% Output:
% irima   : Image with connected regions of phase-values
%           within each range.
% cn      : Total number of conncted regions.
%
% This routine is used to make the initial division into
% a set of regions, which within each it is very unlikely that
% a phase-wrap has occurred, that is the preamble for Mark
% J's method. A higher value for nstep makes it less likely
% that a wrap is included within a region, but will also
% lead to more regions->longer execution time.
%
% N.B. The interval > phi <= is based on the observation that
% angle(-1) returns pi (rather than -pi).
%
% Jesper Andersson 1/10-03
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_initial_regions.m 1317 2008-04-08 16:16:38Z chloe $
if max(abs(pm(:))) > pi
   nstep = round(nstep * (max(pm(:))-min(pm(:))) / (2*pi));
   bins = linspace(min(pm(:))-eps,max(pm(:)),nstep+1);
else
   bins = linspace(-pi,pi,nstep+1);
end

cn = 0;
irima = zeros(size(pm));
for i=1:nstep
   tmp = double((pm > bins(i)) & (pm <= bins(i+1))).*mask;
   [lltmp,num] = spm_bwlabel(tmp,6);
   irima = irima+((lltmp+cn).*(lltmp~=0));
   cn = max(irima(:));
end

return
