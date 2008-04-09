function [type] = senstype(sens, desired)

% SENSTYPE determines the type of sensors by looking at the channel names
% and comparing them with predefined lists.
%
% Use as
%   [type] = senstype(sens)
% to get a string describing the type, or
%   [flag] = senstype(sens, desired)
% to get a boolean value.
%
% The output type can be any of the following
%   'electrode'
%   'ctf151'
%   'ctf151_planar'
%   'ctf275'
%   'ctf275_planar'
%   'neuromag122'
%   'neuromag306'
%   'bti148'
%   'bti148_planar'
%   'bti248'
%   'bti248_planar'
%   'yokogawa160'
%   'yokogawa160_planar'
%   'magnetometer'
%
% The optional input argument for the desired type can be any of the above, plus any of the following
%   'eeg'
%   'meg'
%   'meg_planar'
%   'meg_axial'
%   'ctf'
%   'bti'
%   'neuromag'
%   'yokogawa'
%
% Besides specifiying a grad or elec structure as input, also allowed is
% giving a data structure containing a grad or elec field, or giving a list
% of channel names (as cell-arrray). I.e. assuming a FieldTrip data
% structure, all of the following calls would be correct.
%   senstype(data)
%   senstype(data.grad)
%   senstype(data.grad.label)
%
% See also READ_SENS, COMPUTE_LEADFIELD

% Copyright (C) 2007-2008, Robert Oostenveld
%
% $Log: senstype.m,v $
% Revision 1.8  2008/03/18 13:40:27  roboos
% added quick fix for ctf275_planar, the problem is that one channel is missing from the ctf275 list (i.e. it is only 274 long)
%
% Revision 1.7  2008/03/18 12:25:06  roboos
% use sens.type if available, this requires that the content of sens.type is consistent with the strings returned by this function
% preallocate cell-arrays
%
% Revision 1.6  2008/02/29 15:25:31  roboos
% fixed bug for eeg (thanks to Doug), updated documentation
%
% Revision 1.5  2008/02/29 14:04:42  roboos
% added gradiometers to btiref
% added general types to the checks at the end
%
% Revision 1.4  2008/02/29 13:52:12  roboos
% added bti248 and planar
% changed order of arguments for ismember, needed when a lot of EEG channels are present in ctf151
% changed order in which checks are performed, first ctf275 and only then ctf151
%
% Revision 1.3  2008/02/28 09:16:35  roboos
% fixed bug due to many trigger and headloc channels in ctf275 labels
%
% Revision 1.2  2008/02/27 17:01:54  roboos
% added a whole list of channel names, fixed some bugs
%
% Revision 1.1  2007/07/25 08:31:12  roboos
% implemented new helper function
%

% the input may be a data structure which then contains a sens/elec structure
if isa(sens, 'struct')
  if isfield(sens, 'sens')
    sens = sens.sens;
  elseif isfield(sens, 'elec')
    sens = sens.elec;
  elseif isfield(sens, 'grad')
    sens = sens.grad;
  end
elseif isa(sens, 'cell')
  dum.label = sens; sens = dum; clear dum
end

if isfield(sens, 'type')
  % preferably the structure specifies its own type
  type = sens.type;

else
  % first define some sets with channel names, this will be followed by the comparisons

  btiref = {
    'MRxA'
    'MRyA'
    'MRzA'
    'MLxA'
    'MLyA'
    'MLzA'
    'MCxA'
    'MCyA'
    'MCzA'
    'MRxaA'
    'MRyaA'
    'MRzaA'
    'MLxaA'
    'MLyaA'
    'MLzaA'
    'MCxaA'
    'MCyaA'
    'MCzaA'
    'GxxA'
    'GyxA'
    'GzxA'
    'GyyA'
    'GzyA'
    };

  ctfref = {
    'BG1'
    'BG2'
    'BG3'
    'BP1'
    'BP2'
    'BP3'
    'BR1'
    'BR2'
    'BR3'
    'G11'
    'G12'
    'G13'
    'G22'
    'G23'
    'P11'
    'P12'
    'P13'
    'P22'
    'P23'
    'Q11'
    'Q12'
    'Q13'
    'Q22'
    'Q23'
    'R11'
    'R12'
    'R13'
    'R22'
    'R23'
    };

  ctfheadloc = {
    'HLC0011'
    'HLC0012'
    'HLC0013'
    'HLC0021'
    'HLC0022'
    'HLC0023'
    'HLC0031'
    'HLC0032'
    'HLC0033'
    'HLC0018'
    'HLC0028'
    'HLC0038'
    'HLC0014'
    'HLC0015'
    'HLC0016'
    'HLC0017'
    'HLC0024'
    'HLC0025'
    'HLC0026'
    'HLC0027'
    'HLC0034'
    'HLC0035'
    'HLC0036'
    'HLC0037'
    };

  ctf151 = {
    'MLC11'
    'MLC12'
    'MLC13'
    'MLC14'
    'MLC15'
    'MLC21'
    'MLC22'
    'MLC23'
    'MLC24'
    'MLC31'
    'MLC32'
    'MLC33'
    'MLC41'
    'MLC42'
    'MLC43'
    'MLF11'
    'MLF12'
    'MLF21'
    'MLF22'
    'MLF23'
    'MLF31'
    'MLF32'
    'MLF33'
    'MLF34'
    'MLF41'
    'MLF42'
    'MLF43'
    'MLF44'
    'MLF45'
    'MLF51'
    'MLF52'
    'MLO11'
    'MLO12'
    'MLO21'
    'MLO22'
    'MLO31'
    'MLO32'
    'MLO33'
    'MLO41'
    'MLO42'
    'MLO43'
    'MLP11'
    'MLP12'
    'MLP13'
    'MLP21'
    'MLP22'
    'MLP31'
    'MLP32'
    'MLP33'
    'MLP34'
    'MLT11'
    'MLT12'
    'MLT13'
    'MLT14'
    'MLT15'
    'MLT16'
    'MLT21'
    'MLT22'
    'MLT23'
    'MLT24'
    'MLT25'
    'MLT26'
    'MLT31'
    'MLT32'
    'MLT33'
    'MLT34'
    'MLT35'
    'MLT41'
    'MLT42'
    'MLT43'
    'MLT44'
    'MRC11'
    'MRC12'
    'MRC13'
    'MRC14'
    'MRC15'
    'MRC21'
    'MRC22'
    'MRC23'
    'MRC24'
    'MRC31'
    'MRC32'
    'MRC33'
    'MRC41'
    'MRC42'
    'MRC43'
    'MRF11'
    'MRF12'
    'MRF21'
    'MRF22'
    'MRF23'
    'MRF31'
    'MRF32'
    'MRF33'
    'MRF34'
    'MRF41'
    'MRF42'
    'MRF43'
    'MRF44'
    'MRF45'
    'MRF51'
    'MRF52'
    'MRO11'
    'MRO12'
    'MRO21'
    'MRO22'
    'MRO31'
    'MRO32'
    'MRO33'
    'MRO41'
    'MRO42'
    'MRO43'
    'MRP11'
    'MRP12'
    'MRP13'
    'MRP21'
    'MRP22'
    'MRP31'
    'MRP32'
    'MRP33'
    'MRP34'
    'MRT11'
    'MRT12'
    'MRT13'
    'MRT14'
    'MRT15'
    'MRT16'
    'MRT21'
    'MRT22'
    'MRT23'
    'MRT24'
    'MRT25'
    'MRT26'
    'MRT31'
    'MRT32'
    'MRT33'
    'MRT34'
    'MRT35'
    'MRT41'
    'MRT42'
    'MRT43'
    'MRT44'
    'MZC01'
    'MZC02'
    'MZF01'
    'MZF02'
    'MZF03'
    'MZO01'
    'MZO02'
    'MZP01'
    'MZP02'
    };

  ctf151_planar = cell(151, 2);
  for i=1:151
    ctf151_planar{i,1} = sprintf('%s_dH', ctf151{i});
    ctf151_planar{i,2} = sprintf('%s_dV', ctf151{i});
  end

  ctf275 = {
    'MLC11'
    'MLC12'
    'MLC13'
    'MLC14'
    'MLC15'
    'MLC16'
    'MLC17'
    'MLC21'
    'MLC22'
    'MLC23'
    'MLC24'
    'MLC25'
    'MLC31'
    'MLC32'
    'MLC41'
    'MLC42'
    'MLC51'
    'MLC52'
    'MLC53'
    'MLC54'
    'MLC55'
    'MLC61'
    'MLC62'
    'MLC63'
    'MLF11'
    'MLF12'
    'MLF13'
    'MLF14'
    'MLF21'
    'MLF22'
    'MLF23'
    'MLF24'
    'MLF25'
    'MLF31'
    'MLF32'
    'MLF33'
    'MLF34'
    'MLF35'
    'MLF41'
    'MLF42'
    'MLF43'
    'MLF44'
    'MLF45'
    'MLF46'
    'MLF51'
    'MLF52'
    'MLF53'
    'MLF54'
    'MLF55'
    'MLF56'
    'MLF61'
    'MLF62'
    'MLF63'
    'MLF64'
    'MLF65'
    'MLF66'
    'MLF67'
    'MLO11'
    'MLO12'
    'MLO13'
    'MLO14'
    'MLO21'
    'MLO22'
    'MLO23'
    'MLO24'
    'MLO31'
    'MLO32'
    'MLO33'
    'MLO34'
    'MLO41'
    'MLO42'
    'MLO43'
    'MLO44'
    'MLO51'
    'MLO52'
    'MLO53'
    'MLP11'
    'MLP12'
    'MLP21'
    'MLP22'
    'MLP23'
    'MLP31'
    'MLP32'
    'MLP33'
    'MLP34'
    'MLP35'
    'MLP41'
    'MLP42'
    'MLP43'
    'MLP44'
    'MLP45'
    'MLP51'
    'MLP52'
    'MLP53'
    'MLP54'
    'MLP55'
    'MLP56'
    'MLP57'
    'MLT11'
    'MLT12'
    'MLT13'
    'MLT14'
    'MLT15'
    'MLT16'
    'MLT21'
    'MLT22'
    'MLT23'
    'MLT24'
    'MLT25'
    'MLT26'
    'MLT27'
    'MLT31'
    'MLT32'
    'MLT33'
    'MLT34'
    'MLT35'
    'MLT36'
    'MLT37'
    'MLT41'
    'MLT42'
    'MLT43'
    'MLT44'
    'MLT45'
    'MLT46'
    'MLT47'
    'MLT51'
    'MLT52'
    'MLT53'
    'MLT54'
    'MLT55'
    'MLT56'
    'MLT57'
    'MRC11'
    'MRC12'
    'MRC13'
    'MRC14'
    'MRC15'
    'MRC16'
    'MRC17'
    'MRC21'
    'MRC22'
    'MRC23'
    'MRC24'
    'MRC25'
    'MRC31'
    'MRC32'
    'MRC41'
    'MRC42'
    'MRC51'
    'MRC52'
    'MRC53'
    'MRC54'
    'MRC55'
    'MRC61'
    'MRC62'
    'MRC63'
    'MRF11'
    'MRF12'
    'MRF13'
    'MRF14'
    'MRF21'
    'MRF22'
    'MRF23'
    'MRF24'
    'MRF25'
    'MRF31'
    'MRF32'
    'MRF33'
    'MRF34'
    'MRF35'
    'MRF41'
    'MRF42'
    'MRF43'
    'MRF44'
    'MRF45'
    'MRF46'
    'MRF51'
    'MRF52'
    'MRF53'
    'MRF54'
    'MRF55'
    'MRF56'
    'MRF61'
    'MRF62'
    'MRF63'
    'MRF64'
    'MRF65'
    'MRF66'
    'MRF67'
    'MRO11'
    'MRO12'
    'MRO13'
    'MRO14'
    'MRO21'
    'MRO22'
    'MRO23'
    'MRO24'
    'MRO31'
    'MRO32'
    'MRO33'
    'MRO34'
    'MRO41'
    'MRO42'
    'MRO43'
    'MRO44'
    'MRO51'
    'MRO52'
    'MRO53'
    'MRP11'
    'MRP12'
    'MRP21'
    'MRP22'
    'MRP23'
    'MRP32'
    'MRP33'
    'MRP34'
    'MRP35'
    'MRP41'
    'MRP42'
    'MRP43'
    'MRP44'
    'MRP45'
    'MRP51'
    'MRP52'
    'MRP53'
    'MRP54'
    'MRP55'
    'MRP56'
    'MRP57'
    'MRT11'
    'MRT12'
    'MRT13'
    'MRT14'
    'MRT15'
    'MRT16'
    'MRT21'
    'MRT22'
    'MRT23'
    'MRT24'
    'MRT25'
    'MRT26'
    'MRT27'
    'MRT31'
    'MRT32'
    'MRT33'
    'MRT34'
    'MRT35'
    'MRT36'
    'MRT37'
    'MRT41'
    'MRT42'
    'MRT43'
    'MRT44'
    'MRT45'
    'MRT46'
    'MRT47'
    'MRT51'
    'MRT52'
    'MRT53'
    'MRT54'
    'MRT55'
    'MRT56'
    'MRT57'
    'MZC01'
    'MZC02'
    'MZC03'
    'MZC04'
    'MZF01'
    'MZF02'
    'MZF03'
    'MZO01'
    'MZO02'
    'MZO03'
    'MZP01'
    };

  % f.ck, apparently one channel is missing
  ctf275_planar = cell(274,2);
  for i=1:274
    ctf275_planar{i,1} = sprintf('%s_dH', ctf275{i});
    ctf275_planar{i,2} = sprintf('%s_dV', ctf275{i});
  end

  bti148 = cell(148,1);
  for i=1:148
    bti148{i,1} = sprintf('A%d', i);
  end

  bti148_planar = cell(148,1);
  for i=1:148
    bti148_planar{i,1} = sprintf('A%d_dH', i);
    bti148_planar{i,2} = sprintf('A%d_dV', i);
  end

  bti248 = cell(248,1);
  for i=1:248
    bti248{i,1} = sprintf('A%d', i);
  end

  bti248_planar = cell(248,2);
  for i=1:248
    bti248_planar{i,1} = sprintf('A%d_dH', i);
    bti248_planar{i,2} = sprintf('A%d_dV', i);
  end

  neuromag122 = {
    'MEG 001'    'MEG 002'
    'MEG 003'    'MEG 004'
    'MEG 005'    'MEG 006'
    'MEG 007'    'MEG 008'
    'MEG 009'    'MEG 010'
    'MEG 011'    'MEG 012'
    'MEG 013'    'MEG 014'
    'MEG 015'    'MEG 016'
    'MEG 017'    'MEG 018'
    'MEG 019'    'MEG 020'
    'MEG 021'    'MEG 022'
    'MEG 023'    'MEG 024'
    'MEG 025'    'MEG 026'
    'MEG 027'    'MEG 028'
    'MEG 029'    'MEG 030'
    'MEG 031'    'MEG 032'
    'MEG 033'    'MEG 034'
    'MEG 035'    'MEG 036'
    'MEG 037'    'MEG 038'
    'MEG 039'    'MEG 040'
    'MEG 041'    'MEG 042'
    'MEG 043'    'MEG 044'
    'MEG 045'    'MEG 046'
    'MEG 047'    'MEG 048'
    'MEG 049'    'MEG 050'
    'MEG 051'    'MEG 052'
    'MEG 053'    'MEG 054'
    'MEG 055'    'MEG 056'
    'MEG 057'    'MEG 058'
    'MEG 059'    'MEG 060'
    'MEG 061'    'MEG 062'
    'MEG 063'    'MEG 064'
    'MEG 065'    'MEG 066'
    'MEG 067'    'MEG 068'
    'MEG 069'    'MEG 070'
    'MEG 071'    'MEG 072'
    'MEG 073'    'MEG 074'
    'MEG 075'    'MEG 076'
    'MEG 077'    'MEG 078'
    'MEG 079'    'MEG 080'
    'MEG 081'    'MEG 082'
    'MEG 083'    'MEG 084'
    'MEG 085'    'MEG 086'
    'MEG 087'    'MEG 088'
    'MEG 089'    'MEG 090'
    'MEG 091'    'MEG 092'
    'MEG 093'    'MEG 094'
    'MEG 095'    'MEG 096'
    'MEG 097'    'MEG 098'
    'MEG 099'    'MEG 100'
    'MEG 101'    'MEG 102'
    'MEG 103'    'MEG 104'
    'MEG 105'    'MEG 106'
    'MEG 107'    'MEG 108'
    'MEG 109'    'MEG 110'
    'MEG 111'    'MEG 112'
    'MEG 113'    'MEG 114'
    'MEG 115'    'MEG 116'
    'MEG 117'    'MEG 118'
    'MEG 119'    'MEG 120'
    'MEG 121'    'MEG 122'
    };

  neuromag306 = {
    'MEG 0113'    'MEG 0112'   'MEG 0111'
    'MEG 0122'    'MEG 0123'   'MEG 0121'
    'MEG 0132'    'MEG 0133'   'MEG 0131'
    'MEG 0143'    'MEG 0142'   'MEG 0141'
    'MEG 0213'    'MEG 0212'   'MEG 0211'
    'MEG 0222'    'MEG 0223'   'MEG 0221'
    'MEG 0232'    'MEG 0233'   'MEG 0231'
    'MEG 0243'    'MEG 0242'   'MEG 0241'
    'MEG 0313'    'MEG 0312'   'MEG 0311'
    'MEG 0322'    'MEG 0323'   'MEG 0321'
    'MEG 0333'    'MEG 0332'   'MEG 0331'
    'MEG 0343'    'MEG 0342'   'MEG 0341'
    'MEG 0413'    'MEG 0412'   'MEG 0411'
    'MEG 0422'    'MEG 0423'   'MEG 0421'
    'MEG 0432'    'MEG 0433'   'MEG 0431'
    'MEG 0443'    'MEG 0442'   'MEG 0441'
    'MEG 0513'    'MEG 0512'   'MEG 0511'
    'MEG 0523'    'MEG 0522'   'MEG 0521'
    'MEG 0532'    'MEG 0533'   'MEG 0531'
    'MEG 0542'    'MEG 0543'   'MEG 0541'
    'MEG 0613'    'MEG 0612'   'MEG 0611'
    'MEG 0622'    'MEG 0623'   'MEG 0621'
    'MEG 0633'    'MEG 0632'   'MEG 0631'
    'MEG 0642'    'MEG 0643'   'MEG 0641'
    'MEG 0713'    'MEG 0712'   'MEG 0711'
    'MEG 0723'    'MEG 0722'   'MEG 0721'
    'MEG 0733'    'MEG 0732'   'MEG 0731'
    'MEG 0743'    'MEG 0742'   'MEG 0741'
    'MEG 0813'    'MEG 0812'   'MEG 0811'
    'MEG 0822'    'MEG 0823'   'MEG 0821'
    'MEG 0913'    'MEG 0912'   'MEG 0911'
    'MEG 0923'    'MEG 0922'   'MEG 0921'
    'MEG 0932'    'MEG 0933'   'MEG 0931'
    'MEG 0942'    'MEG 0943'   'MEG 0941'
    'MEG 1013'    'MEG 1012'   'MEG 1011'
    'MEG 1023'    'MEG 1022'   'MEG 1021'
    'MEG 1032'    'MEG 1033'   'MEG 1031'
    'MEG 1043'    'MEG 1042'   'MEG 1041'
    'MEG 1112'    'MEG 1113'   'MEG 1111'
    'MEG 1123'    'MEG 1122'   'MEG 1121'
    'MEG 1133'    'MEG 1132'   'MEG 1131'
    'MEG 1142'    'MEG 1143'   'MEG 1141'
    'MEG 1213'    'MEG 1212'   'MEG 1211'
    'MEG 1223'    'MEG 1222'   'MEG 1221'
    'MEG 1232'    'MEG 1233'   'MEG 1231'
    'MEG 1243'    'MEG 1242'   'MEG 1241'
    'MEG 1312'    'MEG 1313'   'MEG 1311'
    'MEG 1323'    'MEG 1322'   'MEG 1321'
    'MEG 1333'    'MEG 1332'   'MEG 1331'
    'MEG 1342'    'MEG 1343'   'MEG 1341'
    'MEG 1412'    'MEG 1413'   'MEG 1411'
    'MEG 1423'    'MEG 1422'   'MEG 1421'
    'MEG 1433'    'MEG 1432'   'MEG 1431'
    'MEG 1442'    'MEG 1443'   'MEG 1441'
    'MEG 1512'    'MEG 1513'   'MEG 1511'
    'MEG 1522'    'MEG 1523'   'MEG 1521'
    'MEG 1533'    'MEG 1532'   'MEG 1531'
    'MEG 1543'    'MEG 1542'   'MEG 1541'
    'MEG 1613'    'MEG 1612'   'MEG 1611'
    'MEG 1622'    'MEG 1623'   'MEG 1621'
    'MEG 1632'    'MEG 1633'   'MEG 1631'
    'MEG 1643'    'MEG 1642'   'MEG 1641'
    'MEG 1713'    'MEG 1712'   'MEG 1711'
    'MEG 1722'    'MEG 1723'   'MEG 1721'
    'MEG 1732'    'MEG 1733'   'MEG 1731'
    'MEG 1743'    'MEG 1742'   'MEG 1741'
    'MEG 1813'    'MEG 1812'   'MEG 1811'
    'MEG 1822'    'MEG 1823'   'MEG 1821'
    'MEG 1832'    'MEG 1833'   'MEG 1831'
    'MEG 1843'    'MEG 1842'   'MEG 1841'
    'MEG 1912'    'MEG 1913'   'MEG 1911'
    'MEG 1923'    'MEG 1922'   'MEG 1921'
    'MEG 1932'    'MEG 1933'   'MEG 1931'
    'MEG 1943'    'MEG 1942'   'MEG 1941'
    'MEG 2013'    'MEG 2012'   'MEG 2011'
    'MEG 2023'    'MEG 2022'   'MEG 2021'
    'MEG 2032'    'MEG 2033'   'MEG 2031'
    'MEG 2042'    'MEG 2043'   'MEG 2041'
    'MEG 2113'    'MEG 2112'   'MEG 2111'
    'MEG 2122'    'MEG 2123'   'MEG 2121'
    'MEG 2133'    'MEG 2132'   'MEG 2131'
    'MEG 2143'    'MEG 2142'   'MEG 2141'
    'MEG 2212'    'MEG 2213'   'MEG 2211'
    'MEG 2223'    'MEG 2222'   'MEG 2221'
    'MEG 2233'    'MEG 2232'   'MEG 2231'
    'MEG 2242'    'MEG 2243'   'MEG 2241'
    'MEG 2312'    'MEG 2313'   'MEG 2311'
    'MEG 2323'    'MEG 2322'   'MEG 2321'
    'MEG 2332'    'MEG 2333'   'MEG 2331'
    'MEG 2343'    'MEG 2342'   'MEG 2341'
    'MEG 2412'    'MEG 2413'   'MEG 2411'
    'MEG 2423'    'MEG 2422'   'MEG 2421'
    'MEG 2433'    'MEG 2432'   'MEG 2431'
    'MEG 2442'    'MEG 2443'   'MEG 2441'
    'MEG 2512'    'MEG 2513'   'MEG 2511'
    'MEG 2522'    'MEG 2523'   'MEG 2521'
    'MEG 2533'    'MEG 2532'   'MEG 2531'
    'MEG 2543'    'MEG 2542'   'MEG 2541'
    'MEG 2612'    'MEG 2613'   'MEG 2611'
    'MEG 2623'    'MEG 2622'   'MEG 2621'
    'MEG 2633'    'MEG 2632'   'MEG 2631'
    'MEG 2642'    'MEG 2643'   'MEG 2641'
    };

  % start with unknown, then try to determine the proper type by looking at the labels
  type = 'unknown';

  if isfield(sens, 'label') && isfield(sens, 'pnt') && ~isfield(sens, 'ori')
    % looks like EEG, there are no further details needed on the sensor type
    type = 'electrode';

  elseif isfield(sens, 'label')
    % probably this is MEG, determine the type of magnetometer/gradiometer system
    % note that the order here is important: first check whether it matches a 275 channel system, then a 151 channel system, since the 151 channels are a subset of the 275
    if (mean(ismember(ctf275(:),        sens.label)) > 0.8) || (mean(ismember(ctfheadloc(:), sens.label)) > 0.8)
      type = 'ctf275';
    elseif (mean(ismember(ctf151(:),        sens.label)) > 0.8)
      type = 'ctf151';
    elseif (mean(ismember(ctf275_planar(:), sens.label)) > 0.8)
      type = 'ctf275_planar';
    elseif (mean(ismember(ctf151_planar(:), sens.label)) > 0.8)
      type = 'ctf151_planar';
    elseif (mean(ismember(bti248(:),        sens.label)) > 0.8)
      type = 'bti248';
    elseif (mean(ismember(bti148(:),        sens.label)) > 0.8)
      type = 'bti148';
    elseif (mean(ismember(bti248_planar(:), sens.label)) > 0.8)
      type = 'bti248_planar';
    elseif (mean(ismember(bti148_planar(:), sens.label)) > 0.8)
      type = 'bti148_planar';
    elseif (mean(ismember(neuromag306(:),   sens.label)) > 0.8)
      type = 'neuromag306';
    elseif (mean(ismember(neuromag122(:),   sens.label)) > 0.8)
      type = 'neuromag122';
    elseif any(ismember(btiref(:), sens.label))
      type = 'bti'; % it might be 148 or 248 channels
    elseif any(ismember(ctfref(:), sens.label))
      type = 'ctf'; % it might be 151 or 275 channels
    elseif isfield(sens, 'pnt') && isfield(sens, 'ori') && numel(sens.label)==numel(sens.pnt)
      type = 'magnetometer';
    elseif isfield(sens, 'pnt') && isfield(sens, 'ori')
      type = 'meg';
    end
  end

end % if isfield(sens, 'type')

if nargin>1
  % return a boolean flag
  switch desired
    case 'eeg'
      type = any(strcmp(type, {'eeg' 'electrode'}));
    case 'meg'
      type = any(strcmp(type, {'meg' 'magnetometer' 'ctf' 'bti' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar' 'neuromag122' 'neuromag306' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar' 'yokogawa160' 'yokogawa160_planar'}));
    case 'ctf'
      type = any(strcmp(type, {'ctf' 'ctf151' 'ctf275' 'ctf151_planar' 'ctf275_planar'}));
    case 'bti'
      type = any(strcmp(type, {'bti' 'bti148' 'bti148_planar' 'bti248' 'bti248_planar'}));
    case 'neuromag'
      type = any(strcmp(type, {'neuromag122' 'neuromag306'}));
    case 'yokogawa'
      type = any(strcmp(type, {'yokogawa160' 'yokogawa160_planar'}));
    case 'meg_axial'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'magnetometer' 'neuromag306' 'ctf151' 'ctf275' 'bti148' 'bti248' 'yokogawa160'}));
    case 'meg_planar'
      % note that neuromag306 is mixed planar and axial
      type = any(strcmp(type, {'neuromag122' 'neuromag306' 'ctf151_planar' 'ctf275_planar' 'bti148_planar' 'bti248_planar' 'yokogawa160_planar'}));
    otherwise
      type = any(strcmp(type, desired));
  end
end

