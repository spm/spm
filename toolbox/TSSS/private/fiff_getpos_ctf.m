function [R,EX,EY,EZ,CH_NAMES] = fiff_getpos_ctf(file,coord_system)
% Copyright (c) 2016, Elekta Oy
% ---------------------------------------
% 
% Redistribution and use of the Software in source and binary forms, with or without 
% modification, are permitted for non-commercial use.
% 
% The Software is provided "as is" without warranties of any kind, either express or
% implied including, without limitation, warranties that the Software is free of defects,
% merchantable, fit for a particular purpose. Developer/user agrees to bear the entire risk 
% in connection with its use and distribution of any and all parts of the Software under this license.
% 
info = fiff_read_meas_info(file);
meg_chs = [];
if strcmp(coord_system,'head')
  disp('Using head coordinates');
  T = info.dev_head_t.trans 
else
  disp('Using device coordinates');
end
for ch = 1:info.nchan
  if (info.chs(ch).coil_type == 5001) || (info.chs(ch).coil_type == 201609 )
    meg_chs = [meg_chs ch];
  end
end
if strcmp(coord_system,'head')
  r = ones(4,1);
  for ch = 1:length(meg_chs)
    r(1:3) = info.chs(meg_chs(ch)).loc(1:3);
    ex = info.chs(meg_chs(ch)).loc(4:6);
    ey = info.chs(meg_chs(ch)).loc(7:9);
    ez = info.chs(meg_chs(ch)).loc(10:12);
    r_this = T*r;
    R(:,ch) = r_this(1:3);
    EX(:,ch) = T(1:3,1:3)*ex;
    EY(:,ch) = T(1:3,1:3)*ey;
    EZ(:,ch) = T(1:3,1:3)*ez;
  end
else
  for ch = 1:length(meg_chs)
    R(:,ch) = info.chs(meg_chs(ch)).loc(1:3);
    EX(:,ch) = info.chs(meg_chs(ch)).loc(4:6);
    EY(:,ch) = info.chs(meg_chs(ch)).loc(7:9);
    EZ(:,ch) = info.chs(meg_chs(ch)).loc(10:12);
    CH_NAMES{ch} = info.chs(meg_chs(ch)).ch_name;
  end
end

  
