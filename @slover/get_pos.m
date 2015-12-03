function mmpos = get_pos(obj)
% Return point location from last click, in mm
%__________________________________________________________________________

% Matthew Brett
% $Id: get_pos.m 6623 2015-12-03 18:38:08Z guillaume $

mmpos=[];
pos = get(gca, 'CurrentPoint');
u = get(gca, 'UserData');
if mars_struct('isthere', u, 'type')
    if strcmp(u.type, 'slice') % is slice panel
        mmpos = (pos(1,1:2)'-1).*obj.slicedef(:,2)+obj.slicedef(:,1);
        mmpos = obj.transform \ [mmpos; u.no; 1];
        mmpos = mmpos(1:3,1);
    end
end
