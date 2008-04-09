function [vol] = prepare_bemmodel(cfg, mri);

% PREPARE_BEMMODEL constructs triangulations of the boundaries between
% multiple segmented tissue types in an anatomical MRI and subsequently
% computes the BEM system matrix.
%
% Use as
%  [vol] = prepare_bemmodel(cfg, mri), or
%  [vol] = prepare_bemmodel(cfg, vol)
%
% The configuration can contain
%   cfg.tissue         = [1 2 3], segmentation value of each tissue type
%   cfg.numvertices    = [Nskin Nskull Nbrain]
%   cfg.conductivity   = [Cskin Cskull Cbrain]
%   cfg.hdmfile        = string, file containing the volume conduction model (can be empty)
%   cfg.isolatedsource = compartment number, or 0
%   cfg.method         = 'dipoli' or 'brainstorm'
%
% Although the example configuration uses 3 compartments, you can use
% an arbitrary number of compartments.
%
% This function implements 
%   Oostendorp TF, van Oosterom A.
%   Source parameter estimation in inhomogeneous volume conductors of arbitrary shape
%   IEEE Trans Biomed Eng. 1989 Mar;36(3):382-91. 

% Undocumented local options:
% cfg.dipoli

% Copyright (C) 2005, Robert Oostenveld
%
% $Log: prepare_bemmodel.m,v $
% Revision 1.9  2008/03/25 10:56:39  roboos
% use standalone function ama2vol
%
% Revision 1.8  2007/02/13 14:06:19  roboos
% minor change
%
% Revision 1.7  2006/07/24 07:59:16  roboos
% updated documentation
%
% Revision 1.6  2006/07/17 11:04:27  roboos
% show error after dipoli (if any)
%
% Revision 1.5  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.4  2005/12/07 14:41:04  roboos
% added initial support for brainstorms bem_xfer function, it does not work yet
%
% Revision 1.3  2005/11/24 16:27:25  roboos
% added additional default location to search for dipoli
% moved some code around at the end (reading result, cleaning up files)
%
% Revision 1.2  2005/11/07 12:10:32  roboos
% apply the homogenous coordinate transformation on the vertex points
%
% Revision 1.1  2005/11/03 11:15:45  roboos
% new implementation
%

if ~isfield(cfg, 'tissue'),         cfg.tissue = [8 12 14];                  end
if ~isfield(cfg, 'numvertices'),    cfg.numvertices = [1 2 3] * 500;         end
if ~isfield(cfg, 'conductivity'),   cfg.conductivity = [1 1/80 1] * 0.33;    end
if ~isfield(cfg, 'hdmfile'),        cfg.hdmfile = [];                        end
if ~isfield(cfg, 'isolatedsource'), cfg.isolatedsource = [];                 end
if ~isfield(cfg, 'method'),         cfg.method = 'dipoli';                   end

% determine the command-line DIPOLI executable that should be called 
hasdipoli = 0;
if strcmp(cfg.method, 'dipoli') && ~isfield(cfg, 'dipoli')
  target1 = '/Users/roberto/src/thom/bin/dipoli';
  target2 = '/home/coherence/roboos/src/thom/bin/dipoli';
  target3 = '/home/coherence/roboos/bin.linux/dipoli';
  if exist(target1)
    cfg.dipoli = target1;
    hasdipoli = 1;
  elseif exist(target2)
    cfg.dipoli = target2;
    hasdipoli = 1;
  elseif exist(target3)
    cfg.dipoli = target3;
    hasdipoli = 1;
  else
    warning('could not locate dipoli executable');
    hasdipoli = 0;
  end
end

% there are two types of input possible
hasmri = isfield(mri, 'transform');
hasvol = isfield(mri, 'bnd');

if hasvol
  vol = mri;
  clear mri;
else
  vol = [];
end

if ~isfield(vol, 'cond')
  % assign the conductivity of each compartment
  vol.cond = cfg.conductivity;
end

% determine the number of compartments
Ncompartment = length(vol.cond);

if hasmri
  fprintf('using the segmented MRI\n');
  [mrix, mriy, mriz] = ndgrid(1:size(mri.seg,1), 1:size(mri.seg,2), 1:size(mri.seg,3));
  % construct the triangulations of the boundaries from the segmented MRI
  for i=1:Ncompartment
    fprintf('triangulating the boundary of compartment %d\n', i);
    seg = imfill((mri.seg==cfg.tissue(i)), 'holes');
    ori(1) = mean(mrix(seg(:)));
    ori(2) = mean(mriy(seg(:)));
    ori(3) = mean(mriz(seg(:)));
    [pnt, tri] = triangulate_seg(seg, cfg.numvertices(i), ori);
    % apply the coordinate transformation from voxel to head coordinates
    pnt(:,4) = 1;
    pnt = (mri.transform * (pnt'))';
    pnt = pnt(:,1:3);
    vol.bnd(i).pnt = pnt;
    vol.bnd(i).tri = tri;
  end
else
  fprintf('using the pre-specified triangulated boundaries\n');
end

vol.brain = find_innermost_boundary(vol.bnd);
vol.skin  = find_outermost_boundary(vol.bnd);

if isempty(cfg.isolatedsource) && Ncompartment>1
  % the isolated source compartment is by default the most innerone
  cfg.isolatedsource = vol.brain;
elseif isempty(cfg.isolatedsource) && Ncompartment==1
  % the isolated source interface should be contained within at least one other interface
  cfg.isolatedsource = 0;
end

if cfg.isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', cfg.isolatedsource);
else
  fprintf('not using the isolated source approach\n');
end

if strcmp(cfg.method, 'dipoli')
  for i=1:Ncompartment
    bndfile{i} = tempname;
    write_tri(bndfile{i}, vol.bnd(i).pnt, vol.bnd(i).tri);
  end

  if ~isempty(cfg.hdmfile)
    amafile = cfg.hdmfile;
  else
    amafile = tempname;
  end

  exefile = tempname;
  fid = fopen(exefile, 'w');
  fprintf(fid, '#!/bin/sh\n');
  fprintf(fid, '\n');
  fprintf(fid, '%s -i %s << EOF\n', cfg.dipoli, amafile);
  for i=1:Ncompartment
    if i==cfg.isolatedsource
      % the isolated potential approach should be applied using this compartment
      fprintf(fid, '!%s\n', bndfile{i});
    else
      fprintf(fid, '%s\n', bndfile{i});
    end
    fprintf(fid, '%g\n', vol.cond(i));
  end
  fprintf(fid, '\n');
  fprintf(fid, '\n');
  fprintf(fid, 'EOF\n');
  fclose(fid);
  dos(sprintf('chmod +x %s', exefile));

  try
    % execute dipoli and read the resulting file
    system(exefile);
    ama = loadama(amafile);
    vol = ama2vol(ama);
  catch
    warning('an error ocurred while running dipoli');
    disp(lasterr);
  end

  % delete the temporary files
  for i=1:Ncompartment
    delete(bndfile{i})
  end
  delete(exefile);

  % delete the model file, since the user did not make explicit that it should be kept
  if isempty(cfg.hdmfile)
    delete(amafile);
  end

elseif strcmp(cfg.method, 'brainstorm')

  if vol.brain~=1 || vol.skin~=Ncompartment
    % FIXME, it should be possible to rearrange the order of the boundaries here
    error('surfaces should be arranged from innermost to outermost');
  end

  R_eeg = vol.bnd(vol.skin).pnt;       % all vertices of the skin
  R_meg = [];
  O_meg = [];
  for i=1:Ncompartment
    vertices{i} = vol.bnd(i).pnt;
    faces{i}    = vol.bnd(i).tri;
  end
  sigma     = cfg.conductivity;
  mode      = 1;                        % compute EEG only
  basis_opt = 1;                        % linear
  test_opt  = 0;                        % collocation
  ISA       = (cfg.isolatedsource~=0);  % apply ISA
  fn_eeg    = [tempname '.mat'];
  fn_meg    = [tempname '.mat'];
  NVertMax  = [];
  Verbose   = 1;
  checksurf = 1;

  % compute the transfer matrix
  % [Te_ISA,Te,Tm_ISA,Tm,nfv] = bem_xfer(R_eeg,R_meg,O_meg,vertices,faces,sigma,mode,basis_opt,test_opt,ISA,fn_eeg,fn_meg,NVertMax,Verbose,checksurf)
  [Te_ISA] = bem_xfer(R_eeg,R_meg,O_meg,vertices,faces,sigma,mode,basis_opt,test_opt,ISA,fn_eeg,fn_meg,NVertMax,Verbose,checksurf)
  tmp = load(fn_eeg);
  keyboard

end


