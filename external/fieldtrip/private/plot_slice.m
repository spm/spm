function plot_slice(data,varargin)
%
% PLOT_SLICE visualizes the slices of a 3D volume
%
% Use as
%   hs = plot_slice(data, varargin)
%
% Some defaults for the additional arguments:
%
%   'nslices'       = number of slices, (default = 20)
%   'slicerange'    = range of slices in data, (default = 'auto')
%                       'auto', full range of data
%                       [min max], coordinates of first and last slice in voxels
%   'slicedim'      = dimension to slice 1 (x-axis) 2(y-axis) 3(z-axis) (default = 3)
%   'title'         = string, title of the figure window
%   'colorbar'      = 'yes' or 'no' (default = 'yes')
%
% Example
%   figure, plot_slice(data,'title','3D volume','colorbar','no')
%
% Copyright (C) 2009, Cristiano Micheli 
%
% $Log: plot_slice.m,v $
% Revision 1.4  2009/04/20 11:22:07  crimic
% integrated help
%
% Revision 1.3  2009/04/20 11:16:49  crimic
% minor changes and corrrections
%
% Revision 1.2  2009/04/20 09:51:50  crimic
% implemented first version
%

% get the optional input arguments
slicerange  = keyval('slicerange',  varargin); if isempty(slicerange),slicerange='auto';end
nslices     = keyval('nslices',  varargin);  if isempty(nslices),nslices=20;end
slicedim    = keyval('slicedim',  varargin);  if isempty(slicedim),slicedim=3;end
title       = keyval('title',  varargin);  if isempty(title),title='';end
colorbar1   = keyval('colorbar',  varargin);  if isempty(colorbar1),colorbar1='yes';end

if isfield(data,'inside')
  % white BG => mskana
  % % TODO: HERE THE FUNCTION THAT MAKES TO SLICE DIMENSION ALWAYS THE THIRD
  % % DIMENSION, AND ALSO KEEP TRANSFORMATION MATRIX UP TO DATE
  % zoiets
  %if hasana; ana = shiftdim(ana,slicedim-1); end;
  %if hasfun; fun = shiftdim(fun,slicedim-1); end;
  %if hasmsk; msk = shiftdim(msk,slicedim-1); end;
  %%%%% select slices
  if ~isstr(slicerange)
    ind_fslice = slicerange(1);
    ind_lslice = slicerange(2);
  elseif isequal(slicerange, 'auto')
    if hasfun %default
      if isfield(data,'inside')
        ind_fslice = min(find(max(max(data.inside,[],1),[],2)));
        ind_lslice = max(find(max(max(data.inside,[],1),[],2)));
      else
        ind_fslice = min(find(~isnan(max(max(fun,[],1),[],2))));
        ind_lslice = max(find(~isnan(max(max(fun,[],1),[],2))));
      end
    elseif hasana %if only ana, no fun
      ind_fslice = min(find(max(max(ana,[],1),[],2)));
      ind_lslice = max(find(max(max(ana,[],1),[],2)));
    else
      error('no functional parameter and no anatomical parameter, can not plot');
    end
  else
    error('do not understand slicerange');
  end
  ind_allslice = linspace(ind_fslice,ind_lslice,nslices);
  ind_allslice = round(ind_allslice);
  % make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
  if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end;
  %if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

  % update the dimensions of the volume
  if hasana; dim=size(ana); else dim=size(fun); end;

  %%%%% make "quilts", that contain all slices on 2D patched sheet
  % Number of patches along sides of Quilt (M and N)
  % Size (in voxels) of side of patches of Quilt (m and n)
  m = dim(1);
  n = dim(2);
  M = ceil(sqrt(dim(3)));
  N = ceil(sqrt(dim(3)));
  num_patch = N*M;
  if slicedim~=3
    error('only supported for slicedim=3');
  end
  num_slice = (dim(slicedim));
  num_empt = num_patch-num_slice;
  % put empty slides on ana, fun, msk, mskana to fill Quilt up
  if hasana; ana(:,:,end+1:num_patch)=0; end;
  if hasfun; fun(:,:,end+1:num_patch)=0; end;
  if hasmsk; msk(:,:,end+1:num_patch)=0; end;
  %if hasmskana; mskana(:,:,end:num_patch)=0; end;
  % put the slices in the quilt
  for iSlice = 1:num_slice
    xbeg = floor((iSlice-1)./M);
    ybeg = mod(iSlice-1, M);
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(:,:,iSlice));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(:,:,iSlice));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(:,:,iSlice));
    end
    %     if hasmskana
    %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(mskana(:,:,iSlice));
    %     end
  end
  % make vols and scales, containes volumes to be plotted (fun, ana, msk) %added ingnie
  if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
  if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
  if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;

  plot2D(vols2D, scales);

  axis off

  if strcmp(colorbar1,  'yes'),
    if hasfun
      % use a normal Matlab coorbar
      hc = colorbar;
      set(hc, 'YLim', [fcolmin fcolmax]);
    else
      warning('no colorbar possible without functional data')
    end
  end

  if ~isempty(title), title(title); end

  % get the output cfg
  cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

end


function plot2D(vols2D, scales);
cla;
% put 2D volumes in fun, ana and msk
hasana = length(vols2D)>0 && ~isempty(vols2D{1});
hasfun = length(vols2D)>1 && ~isempty(vols2D{2});
hasmsk = length(vols2D)>2 && ~isempty(vols2D{3});

% the transpose is needed for displaying the matrix using the Matlab image() function
if hasana; ana = vols2D{1}'; end;
if hasfun; fun = vols2D{2}'; end;
if hasmsk; msk = vols2D{3}'; end;


if hasana
  % scale anatomy between 0 and 1
  fprintf('scaling anatomy\n');
  amin = min(ana(:));
  amax = max(ana(:));
  ana = (ana-amin)./(amax-amin);
  clear amin amax;
  % convert anatomy into RGB values
  ana = cat(3, ana, ana, ana);
  ha = imagesc(ana);
end
hold on

if hasfun
  hf = imagesc(fun);
  caxis(scales{2});
  % apply the opacity mask to the functional data
  if hasmsk
    % set the opacity
    set(hf, 'AlphaData', msk)
    set(hf, 'AlphaDataMapping', 'scaled')
    alim(scales{3});
  elseif hasana
    set(hf, 'AlphaData', 0.5)
  end
end

axis equal
axis tight
axis xy

function [vols2D] = handle_ortho(vols, indx, slicedir, dim);

% put 2Dvolumes in fun, ana and msk
if length(vols)>=1 && isempty(vols{1}); hasana=0; else ana=vols{1}; hasana=1; end;
if length(vols)>=2
  if isempty(vols{2}); hasfun=0; else fun=vols{2}; hasfun=1; end;
else hasfun=0; end
if length(vols)>=3
  if isempty(vols{3}); hasmsk=0; else msk=vols{3}; hasmsk=1; end;
else hasmsk=0; end

% select the indices of the intersection
xi = indx(1);
yi = indx(2);
zi = indx(3);

% select the slice to plot
if slicedir==1
  yi = 1:dim(2);
  zi = 1:dim(3);
elseif slicedir==2
  xi = 1:dim(1);
  zi = 1:dim(3);
elseif slicedir==3
  xi = 1:dim(1);
  yi = 1:dim(2);
end

% cut out the slice of interest
if hasana; ana = squeeze(ana(xi,yi,zi)); end;
if hasfun; fun = squeeze(fun(xi,yi,zi)); end;
if hasmsk; msk = squeeze(msk(xi,yi,zi)); end;

%put fun, ana and msk in vols2D
if hasana; vols2D{1} = ana; end;
if hasfun; vols2D{2} = fun; end;
if hasmsk; vols2D{3} = msk; end;
