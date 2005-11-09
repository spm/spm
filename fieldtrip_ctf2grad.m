function [grad] = fieldtrip_ctf2grad(hdr, dewar);

% CTF2GRAD converts a CTF header to a gradiometer structure that can be
% understood by FieldTrip and Robert Oostenveld's low-level forward and
% inverse routines. This function only works for headers that have been
% read using the READ_CTF_RES4 function.
%
% Use as:
%   [hdr]  = read_ctf_res4(filename);
%   [grad] = ctf2grad(hdr);
%
% This function only computes the hardware (0th order) gradiometer
% definition for the CTF system. If you want to use higher order synthetic
% gradients, you should use the data that is stored in CTF_COEF, but that
% only applies to the F.C. Donders MEG system)
%
% See also NIMH2GRAD, FIF2GRAD

% undocumented option: it will return the gradiometer information in dewar
% coordinates if second argument is present and non-zero

% Copyright (C) 2004, Robert Oostenveld
%
% $Log: ctf2grad.m,v $
% Revision 1.2  2005/06/01 07:59:37  roboos
% added second argument which optionally causes the gradient information to be returend in dewar coordinates
%
% Revision 1.1  2005/05/26 09:55:17  roboos
% renamed the fileio/ctf_grad function to ctf2grad and moved it into fieldtrip/private for consistency with other gradiometer construction functions (nimh2grad and fif2grad)
%
% Revision 1.1  2004/07/02 11:33:47  roboos
% new function that creates a more complete gradiometer definition from the res4 header (compared to read_ctf_res4)
%

% My preferred ordering in the grad structure is:
%   1st 151 coils are bottom coils of MEG channels
%   2nd 151 are the top coils of MEG channels
%   following coils belong to reference channels

if nargin<2 || isempty(dewar)
  dewar = 0;
end

% start with empty gradiometer
grad = [];
grad.pnt = [];
grad.ori = [];
grad.tra = [];

% combine the bottom and top coil of each MEG channel
numMEG = length(hdr.rowMEG);
for i=1:numMEG
  n = hdr.rowMEG(i);
  % get coil positions and orientations of this channel (max. 8)
  if dewar
    pos = cell2mat({hdr.Chan(n).coil.pos}');
    ori = cell2mat({hdr.Chan(n).coil.ori}');
  else
    pos = cell2mat({hdr.Chan(n).coilHC.pos}');
    ori = cell2mat({hdr.Chan(n).coilHC.ori}');
  end
  % determine the number of coils for this channel
  numcoils = sum(sum(pos.^2, 2)~=0);
  if numcoils~=2
    error('unexpected number of coils in MEG channel');
  end
  % add the coils of this channel to the gradiometer array
  grad.pnt(i       ,:) = pos(1,:);
  grad.pnt(i+numMEG,:) = pos(2,:);
  grad.ori(i       ,:) = ori(1,:) .* -sign(hdr.gainV(n));
  grad.ori(i+numMEG,:) = ori(2,:) .* -sign(hdr.gainV(n));
  grad.tra(i,i       ) = 1;
  grad.tra(i,i+numMEG) = 1;
end

% combine the coils of each reference channel
numREF = length(hdr.rowREF);
for i=1:numREF
  n = hdr.rowREF(i);
  % get coil positions and orientations of this channel (max. 8)
  if dewar
    pos = cell2mat({hdr.Chan(n).coil.pos}');
    ori = cell2mat({hdr.Chan(n).coil.ori}');
  else
    pos = cell2mat({hdr.Chan(n).coilHC.pos}');
    ori = cell2mat({hdr.Chan(n).coilHC.ori}');
  end
  % determine the number of coils for this channel
  numcoils = sum(sum(pos.^2, 2)~=0);
  % add the coils of this channel to the gradiometer array
  for j=1:numcoils
    grad.pnt(end+1, :)     = pos(j,:);
    grad.ori(end+1, :)     = ori(j,:) .* -sign(hdr.gainV(n));
    grad.tra(i+numMEG, end+1) = 1;
  end
end

grad.label = hdr.label([hdr.rowMEG hdr.rowREF]);
grad.unit  = 'cm';
