function filled=fillstruct(def, prov)
% Check fields from provided struct against default struct
% FORMAT filled=fillstruct(def, prov)
% ======
% Input Arguments
%   def  - Struct containing all required fields and their default values
%   prov - Struct array containing fields that will override defaults
% Output Argument
%   filled - Struct array containing field values from prov, if given,
%            otherwise from def struct.
% Only fields which are in def struct are checked and returned. Thus, prov
% can not add new fields to an existing default struct.
% If prov contains an struct array, fields are checked for each
% individual array member and a filled struct array is returned.
%_______________________________________________________________________
%
% @(#) $Id: fillstruct.m 1436 2008-04-16 15:37:03Z guillaume $

rev = '$Revision: 1436 $';
if isempty(prov)
  filled = def;
else
  fnames = fieldnames(def);
  filled(1:numel(prov)) = deal(def);

  for k = 1:numel(fnames)
    if isfield(prov(1),fnames{k})
      % if 1st element has this field, then all in prov array have it
      for l=1:numel(prov)
        filled(l).(fnames{k}) = prov(l).(fnames{k});
      end;
    end;
  end;
end;
return;
