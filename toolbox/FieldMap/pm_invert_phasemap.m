function varargout = pm_invert_phasemap(varargin)
% Inverting phasemaps (trickier than it sounds).
% FORMAT ipm = invert_phasemap(pm) 
% or 
% FORMAT ipm = invert_phasemap(pm,idim) 
% or
% FORMAT ipm = invert_phasemap(P) 
% or
% FORMAT ipm = invert_phasemap(P,idim) 
% or
% FORMAT invert_phasemap(P,fname)
% or
% FORMAT invert_phasemap(P,fname,idim)
%
% Input:
% pm      1, 2 or 3D array representing a displacement field that
%         is to be inverted along one direction.
% idim    The dimension along which field is to be inverted.
% P       File-struct or -name containing displacement field.
% fname   Name of output file.
%
% Output:
% ipm     Displacement-field inverted along requested direction.
%
% This is a gateway function to invert_phasemap_dtj (do the job) 
% which is a mex-file. The job of this routine is to handle some of 
% the basic book-keeping regarding format and file creation.
%_______________________________________________________________________
% Jesper Andersson 10/1-02
% 
% Added the possibility to specify along which direction
% the field should be inverted.
%_______________________________________________________________________
% Jesper Andersson 17/3-05
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson 
% $Id: pm_invert_phasemap.m 4842 2012-08-15 18:02:30Z guillaume $

%
% Decode first input parameter.
%
if isfield(varargin{1},'dim') || ischar(varargin{1})
   if isfield(varargin{1},'dim')
      P = varargin{1}; 
   elseif ischar(varargin{1})
      P = spm_vol(varargin{1});
   end
   pm = spm_read_vols(P);
else
   pm = varargin{1};
end

%
% Decode second input parameter.
%
if nargin > 1
   if ~ischar(varargin{2})
      idim  = varargin{2};
   else
      if nargin > 2
     idim = varargin{3};
      else
         idim = 2;
      end
   end
else
   idim = 2;
end

%
% Take care of the 1D case.
%
if nargin == 1 && length(size(pm)) == 2 && any(size(pm) == 1)
   sz = size(pm);
   if sz(1) > sz(2)
      idim = 1;
   else
      idim = 2;
   end
end


ipm = pm_invert_phasemap_dtj(pm,idim);

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
            if ~isempty(indx) && indx(1) > 1 
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
if length(varargin) == 2 && exist('P','var') == 1
   oP = struct('fname',     varargin{2},...
               'dim',       [dim spm_type('int16')],...
               'mat',       P.mat,...
               'pinfo',     [1 0 0]',...
               'descrip',   'Inverted displacements from phase map');        
   spm_write_vol(oP,ipm);
end

return
