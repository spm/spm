function [time, B] = eb_read_lvm(filename)
if num2str(filename(end-3:end)) == '.lvm'
else
    filename = [filename,'.lvm'];
end
data = dlmread(filename, '\t',23,0);
% data = dlmread(filename, '\t',7,0);
time = data(:,1);
B = data(:,2:end);