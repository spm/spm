% This is an example script that illustrates how to generate a channel
% template file for use in SPM5.
% The example uses easycap's webpage which provides a list of channel names
% and 3D-coordinates for their easycap montage1.

% go to easycap's webpage http://www.easycap.de/easycap/e/downloads/M1_XYZ.txt
% and save text file ('Save Page as')

% the coordinates are in German formatting, therefore you have to convert ',' to
% '.' first. This can be easily done using the find/replace function of the
% matlab editor:
% call editor (must have Java enabled) by 
% edit M1_XYZ.txt
% and find and replace all , by .
% save the modified file. Then:

% read channel names and coordinates
[Cnames, x, y, z] = textread('M1_XYZ.txt', '%s %f %f %f', 'delimiter','\t ', 'headerlines', 1);
coord = [x y z];

% distance of each electrode from Cz (0,0)
d = sqrt(sum(coord'.^2));

% projection to x-y-plane by dividing through distance
% this gives a projection that is (nearly) identical to the plot shown in
% http://www.easycap.de/easycap/e/downloads/M1_XYZ.htm
coord = coord./repmat(d.^2', 1,3);
x = coord(:,1);
y = coord(:,2);

% take min and max of coordinates
mx = max(abs(x));
my = max(abs(y));

% make coordinates lie between 0.05 and 0.95 for x and y
x = x./(2*mx)*0.9 + 0.5;
y = y./(2*my)*0.9 + 0.5;

Cpos = [x y]';

% add veog, heog
Cpos = [Cpos [0.6; 0.95]];
Cnames{end+1} = {'VEOG', 'VEOGR'};
Cpos = [Cpos [0.4; 0.95]];
Cnames{end+1} = 'HEOG';

Nchannels = length(Cnames);

% display ratio of x- and y-axis (used by SPM's display of M/EEG data)
Rxy = 1.5;

save easycap_Montage1 Cpos Cnames Rxy Nchannels
