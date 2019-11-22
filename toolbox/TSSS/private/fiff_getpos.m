function [R,EX,EY,EZ] = fiff_getpos(file,coord_system,calfile,z_offset)
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
if exist('calfile')
  if exist(calfile) ~= 2
    disp(sprintf('Cannot open calibration file %s',calfile));
    clear calfile;
  end
end
info = fiff_read_meas_info(file);
meg_chs = [];
if strcmp(coord_system,'head')
  disp('Using head coordinates');
  T = info.dev_head_t.trans
else
  disp('Using device coordinates'); 
end
if ~exist('z_offset')
  z_offset = 0;
end
for ch = 1:info.nchan
  if info.chs(ch).kind == 1
    meg_chs = [meg_chs ch];
  end
end
if exist('calfile') 
  disp('Reading array information from calibration file');
  load(calfile);
  if strcmp(coord_system,'head')
    for ch = 1:length(meg_chs)
      r(1:3) = info.chs(meg_chs(ch)).loc(1:3);
      r_this = T*r;
      EX(:,ch) = T(1:3,1:3)*NX(:,ch);
      EY(:,ch) = T(1:3,1:3)*NY(:,ch);;
      EZ(:,ch) = T(1:3,1:3)*NZ(:,ch);
      R(:,ch) = R(:,ch) + z_offset*EZ(:,ch);
    end
  else
    EX = NX; EY = NY; EZ = NZ;
    for ch = 1:length(meg_chs)
      R(:,ch) = info.chs(meg_chs(ch)).loc(1:3);
      R(:,ch) = R(:,ch) + z_offset*EZ(:,ch);
    end
  end
else
  disp('Reading array information from measurement file');
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
      R(:,ch) = R(:,ch) + z_offset*EZ(:,ch);
    end
  else
    for ch = 1:length(meg_chs)
      R(:,ch) = info.chs(meg_chs(ch)).loc(1:3);
      EX(:,ch) = info.chs(meg_chs(ch)).loc(4:6);
      EY(:,ch) = info.chs(meg_chs(ch)).loc(7:9);
      EZ(:,ch) = info.chs(meg_chs(ch)).loc(10:12);
      R(:,ch) = R(:,ch) + z_offset*EZ(:,ch);
    end
  end
end
  
