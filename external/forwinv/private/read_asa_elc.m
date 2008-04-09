function elec = read_asa_elc(fn);

% READ_ASA_ELC reads electrodes from an ASA electrode file
% converting the units to mm

% Copyright (C) 2002, Robert Oostenveld
% 
% $Log: read_asa_elc.m,v $
% Revision 1.3  2003/12/16 10:23:55  roberto
% attemt to make it slightly more robust for the electrode labels
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

Npnt = read_asa(fn, 'NumberPositions=', '%d');
Ndhk = read_asa(fn, 'NumberPolygons=', '%d');
Unit = read_asa(fn, 'UnitPosition', '%s', 1);
pnt  = read_asa(fn, 'Positions', '%f', Npnt);
dhk  = read_asa(fn, 'Polygons', '%d', Ndhk);
lab  = read_asa(fn, 'Labels', '%s', Npnt);

if strcmp(lower(Unit),'mm')
  pnt = 1*pnt;
elseif strcmp(lower(Unit),'cm')
  pnt = 100*pnt;
elseif strcmp(lower(Unit),'m')
  pnt = 1000*pnt;
elseif ~isempty(Unit)
  error(sprintf('Unknown unit of distance for electrodes (%s)', Unit));
end

if length(lab)==1 & iscell(lab)
  % the electrode labels were on a single line
  % reformat the electrode labels into an appropriately sized cell array
  remainder = lab{1};
  lab = {};
  for i=1:size(pnt,1)
    [lab{i}, remainder] = strtok(remainder);
  end
end

elec.pnt = pnt;
elec.dhk = dhk+1;
elec.label = lab;
