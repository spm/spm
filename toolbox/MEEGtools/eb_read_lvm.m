function [time, B] = eb_read_lvm(filename,timeind)
if num2str(filename(end-3:end)) == '.lvm'
else
    filename = [filename,'.lvm'];
end
if nargin<2,
    timeind=[];
end;
if isempty(timeind),
    timeind=1; %% default 
end;
data = dlmread(filename, '\t',23,0);
% data = dlmread(filename, '\t',7,0);
time = data(:,timeind);
B = data(:,2:end);