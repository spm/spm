function [varargout] = pm_merge_regions_m(opm,N,P,rima,gd)
%
% Merges regions as defined in connectogram to minimise 
% total costfunction (sum of phase-differences across 
% region borders).
% FORMAT: pm = pm_merge_regions_m(pm,N,P,rima)
% or
% FORMAT: [pm,rima] = pm_merge_regions_m(pm,N,P,rima);
%
% Input:
% pm       : Phase-map
% N        : Sparse matrix where N(i,j) for i<j signify the number
%            of voxels along the common border between regions labelled
%            i and j in rima. Use pm_create_connectogram to get N.
% P        : Sparse matrix where P(i,j) for i<j signify the sum of phase
%            differences along the common border between regions labelled
%            i and j in rima. Use pm_create_connectogram to get N.
% rima     : Label map consisting of connected regions indentified
%            by unique labels. Use pm_initial_regions to get rima.
% gd       : g(raphical)d(isplay) if exist and equals 1 will produce
%            a graphical display of the merging process. It might be useful
%            for getting an understanding of what happens, but the 
%            excitment wears pretty thin pretty soon.
%
% Output:
% pm       : Phase-map after merging of all regions in rima that
%            are connected.
% rima     : Label map after merging of all possible regions (regions 
%        with a common border). Note that if there are dissconnected
%            regions in the original rima (e.g. in 2D where the temporal
%            lobes may be disconnected from the rest of the brain) there
%            will still be more than one label in rima.
%
% This routine is based on the MRM paper by Mark J. Very briefly it will
% use the summary statistic in the matrices N and P, where each entry in
% N signifies the number of voxels along the common border of the regions
% whose labels correspond to row and column of the matrix. E.g. N(i,j) (for i<j) 
% signifies the number of voxels along the border between regions labelled
% i and j. The matrix P is organised in the same manner, with the difference
% that the numbers correspond to the sum of differences of phase values
% across voxel-faces along that border. The direction of the difference has
% been (arbitrarily) chosen such that we take phi(i)-phi(j) where i<j.
%
% Now we want to merge all these regions, such that after merging all
% phase-wraps will have been resolved. An assumption here is that any
% phase-wraps will always be along borders of the initial regions,
% something that is (almost) guaranteed by the way in which we create them.
% 
% There are two aspects to the merging
% 1. We want to detect and correct for any phase-wrap between regions
%    i and j when merging them.
% 2. We want to merge the regions in such an order that more "important"
%    regions are merged first. This is functionally similar to the 
%    progression of wrapping from low->high varinace areas in region-growing
%    approches.
%
% The first goal is easily reached by noting that (P(i,j)/N(i,j))/2pi
% is a good guess for the number of wraps that differ between regions
% i and j.
%
% The second goal is reached by merging the pairs of regions that have
% the largest border (i.e. the largest N(i,j)) first (it is a little
% more elaborate, but basically like that). 
%
% The rest is really just about being really careful when updating the
% stats regarding all the connections between a newly merged regions
% and all the regions that bordered to one or both of the regions
% constituting the new region.
%
% Jenkinson M. 2003. Fast, automated, N-dimensional phase-unwrapping 
% algorithm. MRM 49:193-197.
% 
% This is a .m version of pm_merge_regions.c. It is a fare bit slower
% and produces identical results. Due to its relative simplicity and
% its graphical output capabilities it might however be useful for
% understanding the process per se and for understanding what happens
% if/when unwrapping fails in a certain data set.
%
% If one wants to use the .m versions one should change in pm_unwrap.m
% so that
%
% [ii,jj,nn,pp] = pm_create_connectogram(rima,pm);
% rs = histc(rima(:),[0:max(rima(:))]+0.5);
% rs = rs(1:end-1);
% upm = pm_merge_regions(pm,rima,ii,jj,nn,pp,rs);
%
% changes to
%
% [N,P] = pm_create_connectogram(rima,pm);
% upm = pm_merge_regions_m(pm,N,P,rima);
%
%_________________________________________________________________________
% Jesper Andersson 2/10-03

if nargin < 5
   gd = 0;
end

%
% Get stats on sizes of initial regions.
%
rs = histc(rima(:),[0:max(rima(:))]+0.5);
rs = rs(1:end-1);

cn = length(N);
[i,j,n] = find(N);
[i,j,p] = find(P);
k = -p./(2*pi*n);
l = round(k);
c = 8*pi^2*n.*(.5-abs(k-l));

if gd
   fn = figure('Position',[100 100 900 300]);
end

upm = opm;
for cc=1:(cn-1)

   if gd
      figure(fn);
      subplot(1,3,1); spy(N);
      if length(size(rima))==3
         subplot(1,3,2); imagesc(rot90(squeeze(rima(round(size(rima,1)/2),:,:)),-1)); axis('image');
         subplot(1,3,3); imagesc(rot90(squeeze(upm(round(size(rima,1)/2),:,:)),-1)); axis('image');
         title(sprintf('cc = %d',cc));
      else
         subplot(1,3,2); imagesc(rima); axis('image');
         subplot(1,3,3); imagesc(upm); axis('image');
         title(sprintf('cc = %d',cc));
      end
      drawnow;
   end

   %
   % Find pair of regions with largest
   % "cost" if mis-wrapped.
   %
   [mv,mi] = max(c);

   %
   % "Merge" the regions, using the optimal
   % integer 2pi offset (l). Use the "label" 
   % of the largest region of those constituting 
   % the pair as the label for the new (merged) region.
   %
   if rs(i(mi)) > rs(j(mi))
      mlbl = i(mi); olbl = j(mi); ooff = -l(mi); 
      rs(i(mi)) = rs(i(mi)) + rs(j(mi)); rs(j(mi)) = 0;
   else     
      mlbl = j(mi); olbl = i(mi); ooff = l(mi); 
      rs(j(mi)) = rs(j(mi)) + rs(i(mi)); rs(i(mi)) = 0;
   end
   n(mi) = 0; p(mi) = 0;
   
   %
   % Update stats of interfaces to old label
   %
   olbli = find(i==olbl);
   olblj = find(j==olbl);
   i(olbli) = mlbl;
   j(olblj) = mlbl;
   tmp = sort([i([olbli; olblj]) j([olbli; olblj])],2);
   %
   % Sign reverse effect of phase-unwrapping when old-label
   % has the higher region number of the pair.
   %
   signrev1 = [ones(length(olbli),1); -1*ones(length(olblj),1)];
   %
   % Sign reverse contribution from old pair when row-column
   % trades place when old label is traded for new label.
   %
   signrev2 = 2*(i([olbli; olblj])==tmp(:,1)) - 1;
   i([olbli; olblj]) = tmp(:,1);
   j([olbli; olblj]) = tmp(:,2);
   p([olbli; olblj]) = p([olbli; olblj]) + 2*pi*ooff*signrev1.*n([olbli; olblj]);
   p([olbli; olblj]) = signrev2.*p([olbli; olblj]);
   %
   % Merge interfaces to old and new labels
   %
   N = sparse(i,j,n,cn,cn);
   P = sparse(i,j,p,cn,cn);
   [i,j,n] = find(N);
   [i,j,p] = find(P);
   %
   % Calculate new costs for region pairs
   %
   k = -p./(2*pi*n);
   l = round(k);
   c = 8*pi^2*n.*(.5-abs(k-l));
   
   %
   % Make changes in label image and phase-map
   %
   indx = find(rima==olbl); 
   rima(indx) = mlbl;
   upm(indx) = upm(indx) + 2*pi*ooff;

   %
   % Bail out when there are no more connections.
   %
   if length(i) == 0
      break;
   end
end   

varargout{1} = upm;
if nargout > 1
   varargout{2} = rima;
end

return






