function varargout = invert_phasemap(varargin)
% Inverting phasemaps (trickier than it sounds).
% FORMAT ipm = invert_phasemap(pm) or 
% FORMAT ipm = invert_phasemap(P) or
% FORMAT invert_phasemap(P,fname)
%
% This is a gateway function to invert_phasemap_dtj 
% (do the job) which is a mex-file. The job of this
% routine is to handle some of the basic book-keeping
% regarding format and file creation.
%_______________________________________________________________________
% Jesper Andersson 10/1-02

if isfield(varargin{1},'dim') | ischar(varargin{1})
   if isfield(varargin{1},'dim')
      P = varargin{1}; 
   elseif ischar(varargin{1})
      P = spm_vol(varargin{1});
   end
   dim = P.dim(1:3);
   pm = spm_read_vols(P);
else
   pm = varargin{1};
   dim = size(pm); 
end

if length(dim) == 2
   dim(3) = 1;
   sz = size(pm);
   if sz(2) == 1
      % pm = pm'; % Not needed with mex-version
      dim(1:2) = dim(2:-1:1);
   end
end

%
% If dfield is a 2D or 3D matrix inversion is going
% to occurr along the second dimension. If it is a
% vector (row or column) inversion is done along
% (predictably) the only non-singelton dimension.
%

ipm = pm_invert_phasemap_dtj(pm,dim);

%
% The next (dead) section shows the implementation in
% Matlab code for documentation purposes.
%
if 1==0 
   ipm = zeros(dim(1),dim(2),dim(3));
   y = zeros(1,dim(2));
   for sl = 1:dim(3)
      for col = 1:dim(1)
         gy = [1:dim(2)]+pm(col,:,sl);
         for i=1:dim(2)
            indx = find(gy > i);
            if ~isempty(indx) & indx(1) > 1 
               y(i) = (indx(1)-1) + (i-gy(indx(1)-1))*1/(gy(indx(1))-gy(indx(1)-1));
            else
               y(i) = NaN;
            end
         end
         y = y-[1:dim(2)];
         %
         % Let us assume nearest neighbour value for NaNs.
         %
         indx = find(~isnan(y));
         for i=1:indx(1)-1 y(i) = y(indx(1)); end
         for i=indx(end)+1:length(y) y(i) = y(indx(end)); end
         ipm(col,:,sl) = y;
      end
   end
end

%
% The remaining code is alive.
%

if nargout == 1
   varargout{1} = ipm;
end
if length(varargin) == 2 & exist(P) == 1
   oP = struct('fname',     varargin{2},...
               'dim',       [dim spm_type('int16')],...
               'mat',       P.mat,...
               'pinfo',     [1 0 0]',...
               'descrip',   'Inverted displacements from phase map');        
   spm_write_vol(oP,ipm);
end

return



