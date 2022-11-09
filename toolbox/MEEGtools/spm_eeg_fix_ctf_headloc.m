function D = spm_eeg_fix_ctf_headloc(S)
% Fix head localization data in a continuous CTF dataset with continuous
% head localization. The tracking has to be valid at least some of the time
%
% The functionality requires the original CTF header (read with CTF toolbox)
% to be present (set S.saveorigheader = 1 at conversion).
%
% FORMAT D = spm_eeg_fix_ctf_headloc(S)
%
% S         - struct (optional)
% (optional) fields of S:
% S.D - meeg object or filename
%
%
% Output:
% D         - MEEG data struct or cell array of MEEG objects with the
%             rejected trials set to bad and sensors corrected (if
%             requested).
%
% Disclaimer: this code is provided as an example and is not guaranteed to work
% with data on which it was not tested. If it does not work for you, feel
% free to improve it and contribute your improvements to the MEEGtools toolbox
% in SPM (http://www.fil.ion.ucl.ac.uk/spm)
%
%__________________________________________________________________________

% Vladimir Litvak, Robert Oostenveld
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Fix CTF head locations',0);

if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select M/EEG mat file');
    S.D = D;
end

D = spm_eeg_load(S.D);

if ~isfield(S, 'mode'), S.mode = 'interpolate'; end

%read HLC-channels
%HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
%HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
%HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
hlc_chan_label = {'HLC0011' 'HLC0012' 'HLC0013'...
    'HLC0021' 'HLC0022' 'HLC0023'...
    'HLC0031' 'HLC0032' 'HLC0033'};

if ~all(ismember(hlc_chan_label, D.chanlabels))
    error('Head localization is not supported for this data')
end

if ~isfield(D, 'origheader')
    error('Original CTF header needs to be present')
end

if ~isfield(S, 'quickfix')
    
    %read HLC-channels
    %HLC0011 HLC0012 HLC0013 x, y, z coordinates of nasion-coil in m.
    %HLC0021 HLC0022 HLC0023 x, y, z coordinates of lpa-coil in m.
    %HLC0031 HLC0032 HLC0033 x, y, z coordinates of rpa-coil in m.
    hlc_chan_label = {'HLC0011' 'HLC0012' 'HLC0013'...
        'HLC0021' 'HLC0022' 'HLC0023'...
        'HLC0031' 'HLC0032' 'HLC0033'};
    
    
    tmpdat  = D(D.indchannel(hlc_chan_label), :);
    
    tmpind = find(~any(tmpdat==0));
    
    cont_fid  = permute(reshape(tmpdat(:, tmpind)', [], 3, 3), [1 3 2]);
    
    if isfield(S, 'valid_fid')
        if isequal(S.valid_fid, 1)
            valid_fid = 0.01*D.origheader.hc.dewar';
        else
            valid_fid = S.valid_fid;             
        end
        
        dewar_fid = 0.01*D.origheader.hc.dewar';
        dewar_dist = [norm(dewar_fid(1,:) - dewar_fid(2,:));
            norm(dewar_fid(2,:) - dewar_fid(3,:));
            norm(dewar_fid(3,:) - dewar_fid(1,:))];
        
        if numel(valid_fid) == 9
            valid_dist = [norm(valid_fid(1,:) - valid_fid(2,:));
                norm(valid_fid(2,:) - valid_fid(3,:));
                norm(valid_fid(3,:) - valid_fid(1,:))];
        elseif numel(valid_fid) == 3
            valid_dist = valid_fid;
        else
            error('Unexpected input');
        end
        
        if max(abs(dewar_dist-valid_dist))<0.01
            dewar_valid = 1;
        else
            dewar_valid = 0;
        end
                    
        dist_dev = [
            (squeeze(sqrt(sum((cont_fid(:, 1, :) - cont_fid(:, 2, :)).^2, 3))) - valid_dist(1))';...
            (squeeze(sqrt(sum((cont_fid(:, 2, :) - cont_fid(:, 3, :)).^2, 3))) - valid_dist(2))';...
            (squeeze(sqrt(sum((cont_fid(:, 3, :) - cont_fid(:, 1, :)).^2, 3))) - valid_dist(3))'];
        
        W = abs(dist_dev) < 0.01;
    else
        %%
        dist = [
            squeeze(sqrt(sum((cont_fid(:, 1, :) - cont_fid(:, 2, :)).^2, 3)))';...
            squeeze(sqrt(sum((cont_fid(:, 2, :) - cont_fid(:, 3, :)).^2, 3)))';...
            squeeze(sqrt(sum((cont_fid(:, 3, :) - cont_fid(:, 1, :)).^2, 3)))' ];
        %%
        rdist = (round(100*dist));
        
        M = mode(rdist');
        
        W = (abs(dist - repmat(0.01*M(:), 1, size(rdist, 2))) < 0.01);
        
        dewar_valid = 1;
    end
    
    if sum(sum(W) == 3) > 2
        W(:, sum(W)<3) = 0;
    end
    
    % Theoretically one should use 'or' here and not 'and' but I found it
    % leaves too much rubbish in
    nasOK = W(1, :) & W(3, :);
    leOK =  W(1, :) & W(2, :);
    reOK =  W(2, :) & W(3, :);
    
    OK = [nasOK;nasOK;nasOK; leOK; leOK;leOK; reOK;reOK;reOK];
    
    if min(sum(OK, 2))>2               
        fixed  = tmpdat;
        marked = tmpdat;

        for i = 1:size(tmpdat, 1)
            if any(~OK(i, :)) || (length(tmpind)<size(tmpdat, 2))                
                fixed(i, :) = interp1(D.time(tmpind(OK(i, :))), tmpdat(i, tmpind(OK(i, :))),  D.time, 'linear', 'extrap');
                marked(i, tmpind(~OK(i, :))) = NaN;
            end
        end
    else
        if dewar_valid
            fixed = repmat(0.01*D.origheader.hc.dewar(:), 1, size(tmpdat, 2));
            warning('Not enough valid head localization data');
        elseif exist('valid_fid', 'var') && numel(valid_fid) == 9
            fixed = repmat(reshape(valid_fid', [], 1), 1, size(tmpdat, 2));
            warning('No valid information in the dataset. Using the valid fiducials provided.');
        else
            error('There is no way to fix the head location data')
        end       
    end
    
    switch S.mode
        case 'interpolate'
            D(D.indchannel(hlc_chan_label), :) = fixed;
        case 'mark'
            D(D.indchannel(hlc_chan_label), :) = marked;
    end
    
    dewarfid = 100*reshape(median(fixed(:, tmpind), 2), 3, 3);
    %%
    D.origheader.hc.dewar = dewarfid;
    
    spm_figure('GetWin','Graphics');clf;
    subplot(2, 1, 1);
    plot(D.time, tmpdat, 'Color', 0.5*[1 1 1], 'LineWidth', 5);
    hold on
    
    switch S.mode
        case 'interpolate'
            plot(D.time, fixed, 'r');
        case 'mark'
            plot(D.time, marked, 'r');
    end
    
    ylim([min(fixed(:)) max(fixed(:))]);
    subplot(2, 1, 2);
else
    spm_figure('GetWin','Graphics');clf;
    dewarfid = D.origheader.hc.dewar;
end

dewargrad = ctf2grad(D.origheader, 1);

M = spm_eeg_inv_headcoordinates(dewarfid(:, 1), dewarfid(:, 2), dewarfid(:, 3));

grad = ft_convert_units(ft_transform_geometry(M, dewargrad), 'mm');

grad.coordsys = 'ctf';

D = sensors(D, 'MEG', grad);

fid = D.fiducials;
fid.pnt = [];
fid.fid.pnt = dewarfid';
fid.unit = 'cm';

fid = ft_convert_units(ft_transform_geometry(M, fid), 'mm');

D = fiducials(D, fid);

save(D);

plot3(grad.coilpos(:, 1), grad.coilpos(:, 2), grad.coilpos(:, 3), 'k.');
axis equal off
hold on
plot3(fid.fid.pnt(:, 1), fid.fid.pnt(:, 2), fid.fid.pnt(:, 3), 'r.', 'MarkerSize', 15);


function [grad, elec] = ctf2grad(hdr, dewar, coilaccuracy)

% CTF2GRAD converts a CTF header to a gradiometer structure that can be understood by
% the FieldTrip low-level forward and inverse routines. The fieldtrip/fileio
% read_header function can use three different implementations of the low-level code
% for CTF data. Each of these implementations is dealt with here.
%
% Use as
%   [grad, elec] = ctf2grad(hdr, dewar, coilaccuracy)
% where
%   dewar        = boolean, whether to return it in dewar or head coordinates (default is head coordinates)
%   coilaccuracy = empty or a number (default is empty)
%
% See also BTI2GRAD, FIF2GRAD, MNE2GRAD, ITAB2GRAD, YOKOGAWA2GRAD,
% FT_READ_SENS, FT_READ_HEADER

% Copyright (C) 2004-2017, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.

if nargin<2 || isempty(dewar)
  dewar = false;
end

if nargin<3 || isempty(coilaccuracy)
  % if empty it will use the original code
  % otherwise it will use the specified accuracy coil definition from the MNE coil_def.dat
  coilaccuracy = [];
end

if isfield(hdr, 'orig')
  hdr = hdr.orig; % use the original CTF header, not the FieldTrip header
end

% start with empty gradiometer structure
grad = [];
grad.coilpos  = [];
grad.coilori  = [];
grad.tra      = [];

% start with empty electrode structure
elec = [];

% meg channels are 5, refmag 0, refgrad 1, ADCs 18
% UPPT001 is 11
% UTRG001 is 11
% SCLK01  is 17
% STIM    is 11
% SCLK01  is 17
% EEG057  is 9
% ADC06   is 18
% ADC07   is 18
% ADC16   is 18
% V0      is 15

if ~isempty(coilaccuracy)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % use the coil definitions from the MNE coil_def.dat file
  % these allow for varying accuracy which is specified by
  % coilaccuracy = 0, 1 or 2
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  ft_hastoolbox('mne', 1);
  [ftver, ftpath] = ft_version;
  def = mne_load_coil_def(fullfile(ftpath, 'external', 'mne', 'coil_def.dat'));
  
  k = 1;
  for i=1:length(hdr.res4.senres)
    switch hdr.res4.senres(i).sensorTypeIndex
      case 5 % 5001
        thisdef = def([def.id]==5001 & [def.accuracy]==coilaccuracy);
      case 0 % 5002
        thisdef = def([def.id]==5002 & [def.accuracy]==coilaccuracy);
      case 1 % 5003
        thisdef = def([def.id]==5003 & [def.accuracy]==coilaccuracy);
      otherwise
        % do not add this as sensor to the gradiometer definition
        thisdef = [];
    end % case
    
    if ~isempty(thisdef)
      % the sensors (i.e. the combination of coils that comprises a channel) is
      % defined at [0 0 0] and with the direction [0 0 1]
      pos0 = [0 0 0];
      ori0 = [0 0 1];
      
      if dewar
        % convert from cm to m
        pos2 = hdr.res4.senres(i).pos0(:,1)' / 100;
        if hdr.res4.senres(i).numCoils==2
          % determine the direction from the position of the two coils
          ori2 = (hdr.res4.senres(i).pos0(:,2) - hdr.res4.senres(i).pos0(:,1))';
        else
          % determine the direction from the orientation of the coil
          ori2 = hdr.res4.senres(i).ori0(:,1);
        end
      else
        % convert from cm to m
        pos2 = hdr.res4.senres(i).pos(:,1)' / 100;
        if hdr.res4.senres(i).numCoils==2
          % determine the direction from the position of the two coils
          ori2 = (hdr.res4.senres(i).pos(:,2) - hdr.res4.senres(i).pos(:,1))';
        else
          % determine the direction from the orientation of the coil
          ori2 = hdr.res4.senres(i).ori(:,1);
        end
      end
      ori2 = ori2/norm(ori2);
      
      for j=1:thisdef.num_points
        weight = thisdef.coildefs(j,1);
        pos1 = thisdef.coildefs(j,2:4);
        ori1 = thisdef.coildefs(j,5:7);
        
        [az0,el0,r0] = cart2sph(ori0(1), ori0(2), ori0(3));
        [az2,el2,r2] = cart2sph(ori2(1), ori2(2), ori2(3));
        % determine the homogenous transformation that rotates [1 0 0] towards ori0
        R0 = rotate([0 0 az0*180/pi]) * rotate([0 -el0*180/pi 0]);
        % determine the homogenous transformation that rotates [1 0 0] towards ori2
        R2 = rotate([0 0 az2*180/pi]) * rotate([0 -el2*180/pi 0]);
        % determine the homogenous transformation that rotates ori0 to ori2
        R = R2/R0;
        T = translate(pos2 - pos0);
        H = T*R;
        grad.tra(i,k)     = weight;
        grad.coilpos(k,:) = ft_warp_apply(T*R, pos1); % first the rotation, then the translation
        grad.coilori(k,:) = ft_warp_apply(  R, ori1); % only the rotation
        k = k+1;
      end % for num_points
      
      grad.chanpos(i,:) = pos2;
      grad.chanori(i,:) = ori2;
      % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
      grad.label{i} = strtok(hdr.res4.chanNames(i,:), '-');
    end % if MEG or MEGREF
  end % for each channel
  
  grad.label = grad.label(:);
  grad.unit  = 'm'; % the coil_def.dat file is in meter
  
  remove = cellfun(@isempty, grad.label);
  grad.label   = grad.label(~remove);
  grad.tra     = grad.tra(~remove,:);
  grad.chanpos = grad.chanpos(~remove,:);
  grad.chanori = grad.chanori(~remove,:);
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
elseif isfield(hdr, 'res4') && isfield(hdr.res4, 'senres')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF p-files, i.e. readCTFds
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  sensType  = [hdr.res4.senres.sensorTypeIndex];
  selMEG    = find(sensType==5);
  selREF    = find(sensType==0 | sensType==1);
  selEEG    = find(sensType==9);
  selMEG    = selMEG(:)';
  selREF    = selREF(:)';
  selEEG    = selEEG(:)';
  numMEG    = length(selMEG);
  numREF    = length(selREF);
  numEEG    = length(selEEG);
  
  % determine the number of channels and coils
  coilcount = 0;
  coilcount = coilcount + sum([hdr.res4.senres(selREF).numCoils]);
  coilcount = coilcount + sum([hdr.res4.senres(selMEG).numCoils]);
  chancount = numMEG + numREF;
  % preallocate the memory
  grad.coilpos = zeros(coilcount, 3);         % this will hold the position of each coil
  grad.coilori = zeros(coilcount, 3);         % this will hold the orientation of each coil
  grad.tra     = zeros(chancount, coilcount); % this describes how each coil contributes to each channel
  
  if numEEG>0
    for i=1:numEEG
      n = selEEG(i);
      if dewar
        pos = hdr.res4.senres(n).pos0';
      else
        pos = hdr.res4.senres(n).pos';
      end
      if hdr.res4.senres(n).numCoils~=1
        ft_error('unexpected number of electrode positions in EEG channel');
      end
      % add this position
      elec.elecpos(i       ,:) = pos(1,:);
    end
    % add the electrode names
    elec.label = cellstr(hdr.res4.chanNames(selEEG,:));
  else
    elec = [];
  end
  
  % combine the bottom and top coil of each MEG channel
  for i=1:numMEG
    n = selMEG(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = hdr.res4.senres(n).pos0';
      ori = hdr.res4.senres(n).ori0';
    else
      pos = hdr.res4.senres(n).pos';
      ori = hdr.res4.senres(n).ori';
    end
    if hdr.res4.senres(n).numCoils~=2
      ft_error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.coilpos(i       ,:) = pos(1,:);
    grad.coilpos(i+numMEG,:) = pos(2,:);
    grad.coilori(i       ,:) = ori(1,:) .* -sign(hdr.res4.senres(n).properGain);
    grad.coilori(i+numMEG,:) = ori(2,:) .* -sign(hdr.res4.senres(n).properGain);
    grad.tra(i,i       ) = 1;
    grad.tra(i,i+numMEG) = 1;
  end
  
  % the MEG channels always have 2 coils, the reference channels vary in the number of coils
  chancount = 1*numMEG;
  coilcount = 2*numMEG;
  
  % combine the coils of each reference channel
  for i=1:numREF
    n = selREF(i);
    % get coil positions and orientations of this channel (max. 8)
    if dewar
      pos = hdr.res4.senres(n).pos0';
      ori = hdr.res4.senres(n).ori0';
    else
      pos = hdr.res4.senres(n).pos';
      ori = hdr.res4.senres(n).ori';
    end
    % determine the number of coils for this channel
    numcoils = hdr.res4.senres(n).numCoils;
    % add the coils of this channel to the gradiometer array
    chancount = chancount+1;
    for j=1:numcoils
      coilcount = coilcount+1;
      grad.coilpos(coilcount, :)         = pos(j,:);
      grad.coilori(coilcount, :)         = ori(j,:) .* -sign(hdr.res4.senres(n).properGain);
      grad.tra(chancount, coilcount) = 1;
    end
  end
  
  label = cellstr(hdr.res4.chanNames);
  for i=1:numel(label)
    % remove the site-specific numbers from each channel name, e.g. 'MZC01-1706' becomes 'MZC01'
    label{i} = strtok(label{i}, '-');
  end
  
  grad.label = label([selMEG selREF]);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
  % convert the balancing coefficients into a montage that can be used with the ft_apply_montage function
  if isfield(hdr.BalanceCoefs, 'G1BR')
    meglabel          = label(hdr.BalanceCoefs.G1BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G1BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G1BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G1BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G2BR')
    meglabel          = label(hdr.BalanceCoefs.G2BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G2BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G2BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G2BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G3BR')
    meglabel          = label(hdr.BalanceCoefs.G3BR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G3BR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G3BR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G3BR = montage;
  end
  
  if isfield(hdr.BalanceCoefs, 'G3AR')
    meglabel          = label(hdr.BalanceCoefs.G3AR.MEGlist);
    reflabel          = label(hdr.BalanceCoefs.G3AR.Refindex);
    nmeg              = length(meglabel);
    nref              = length(reflabel);
    montage.labelold  = cat(1, meglabel, reflabel);
    montage.labelnew  = cat(1, meglabel, reflabel);
    montage.tra       = [eye(nmeg, nmeg), -hdr.BalanceCoefs.G3AR.alphaMEG'; zeros(nref, nmeg), eye(nref, nref)];
    grad.balance.G3AR = montage;
  end
  
  if     all([hdr.res4.senres(selMEG).grad_order_no]==0)
    grad.balance.current = 'none';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==1)
    grad.balance.current = 'G1BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==2)
    grad.balance.current = 'G2BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==3)
    grad.balance.current = 'G3BR';
  elseif all([hdr.res4.senres(selMEG).grad_order_no]==13)
    grad.balance.current = 'G3AR';
  else
    ft_warning('cannot determine balancing of CTF gradiometers');
    grad = rmfield(grad, 'balance');
  end
  
  % sofar the gradiometer definition was the ideal, non-balenced one
  if isfield(grad, 'balance') && ~strcmp(grad.balance.current, 'none')
    % apply the current balancing parameters to the gradiometer definition
    %grad = ft_apply_montage(grad, getfield(grad.balance, grad.balance.current));
    grad = ft_apply_montage(grad, getfield(grad.balance, grad.balance.current), 'keepunused', 'yes');
  end
  
  
elseif isfield(hdr, 'sensType') && isfield(hdr, 'Chan')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the open-source MATLAB code that originates from CTF and that was modified by the FCDC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  selMEG = find(hdr.sensType==5);
  selREF = find(hdr.sensType==0 | hdr.sensType==1);
  selMEG = selMEG(:)';
  selREF = selREF(:)';
  numMEG = length(selMEG);
  numREF = length(selREF);
  
  % combine the bottom and top coil of each MEG channel
  for i=1:numMEG
    n = selMEG(i);
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
      ft_error('unexpected number of coils in MEG channel');
    end
    % add the coils of this channel to the gradiometer array
    grad.coilpos(i       ,:) = pos(1,:);
    grad.coilpos(i+numMEG,:) = pos(2,:);
    grad.coilori(i       ,:) = ori(1,:) .* -sign(hdr.gainV(n));
    grad.coilori(i+numMEG,:) = ori(2,:) .* -sign(hdr.gainV(n));
    grad.tra(i,i       ) = 1;
    grad.tra(i,i+numMEG) = 1;
  end
  numMEGcoils = size(grad.coilpos, 1);
  
  % combine the coils of each reference channel
  for i=1:numREF
    n = selREF(i);
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
      grad.coilpos(numMEGcoils+i, :)     = pos(j,:);
      grad.coilori(numMEGcoils+i, :)     = ori(j,:) .* -sign(hdr.gainV(n));
      grad.tra(numMEG+i, numMEGcoils+i) = 1;
    end
  end
  
  grad.label = hdr.label([selMEG selREF]);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  if dewar
    grad.coordsys = 'dewar';
  else
    grad.coordsys = 'ctf';
  end
  
elseif isfield(hdr, 'sensor') && isfield(hdr.sensor, 'info')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % the header was read using the CTF importer from the NIH and Daren Weber
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if dewar
    % this does not work for Daren Webers implementation
    ft_error('cannot return the gradiometer definition in dewar coordinates');
  end
  
  % only work on the MEG channels
  if isfield(hdr.sensor.index, 'meg')
    sel = hdr.sensor.index.meg;
  else
    sel = hdr.sensor.index.meg_sens;
  end
  
  for i=1:length(sel)
    pnt = hdr.sensor.info(sel(i)).location';
    ori = hdr.sensor.info(sel(i)).orientation';
    numcoils(i) = size(pnt,1);
    if size(ori,1)==1 && size(pnt,1)==1
      % one coil position with one orientation: magnetometer
      ori = ori;
    elseif size(ori,1)==1 && size(pnt,1)==2
      % two coil positions with one orientation: first order gradiometer
      % assume that the orientation of the upper coil is opposite to the lower coil
      ori = [ori; -ori];
    else
      ft_error('do not know how to deal with higher order gradiometer hardware')
    end
    
    % add this channels coil positions and orientations
    grad.coilpos = [grad.coilpos; pnt];
    grad.coilori = [grad.coilori; ori];
    grad.label{i} = hdr.sensor.info(sel(i)).label;
    
    % determine the contribution of each coil to each channel's output signal
    if size(pnt,1)==1
      % one coil, assume that the orientation is correct, i.e. the weight is +1
      grad.tra(i,end+1) = 1;
    elseif size(pnt,1)==2
      % two coils, assume that the orientation for each coil is correct,
      % i.e. the weights are +1 and +1
      grad.tra(i,end+1) = 1;
      grad.tra(i,end+1) = 1;
    else
      ft_error('do not know how to deal with higher order gradiometer hardware')
    end
  end
  
  % prefer to have the labels in a column vector
  grad.label = grad.label(:);
  grad.unit  = 'cm'; % the res4 file represents it in centimeter
  
  % reorder the coils, such that the bottom coils are at the first N
  % locations and the top coils at the last N positions. This makes it
  % easier to use a selection of the coils for topographic plotting
  if all(numcoils==2)
    bot = 1:2:sum(numcoils);
    top = 2:2:sum(numcoils);
    grad.coilpos = grad.coilpos([bot top], :);
    grad.coilori = grad.coilori([bot top], :);
    grad.tra = grad.tra(:, [bot top]);
  end
  
else
  ft_error('unknown header to contruct gradiometer definition');
end

% chantype and chanunit are required as of 2016, see FT_DATATYPE_SENS
if ~isempty(grad)
  grad.chantype = ft_chantype(grad);
  grad.chanunit = ft_chanunit(grad);
end
if ~isempty(elec)
  elec.chantype = ft_chantype(elec);
  elec.chanunit = ft_chanunit(elec);
end
