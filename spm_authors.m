function [current, previous] = spm_authors
% Return list of SPM coauthors
% FORMAT [current, previous] = spm_authors
% current  - cell array of SPM coauthors of the current release
% previous - cell array of SPM coauthors of previous releases
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% SPM
% $Id: spm_authors.m 5039 2012-11-06 20:39:58Z guillaume $


fid = fopen(fullfile(spm('Dir'),'AUTHORS.txt'),'rt');
if fid==-1, current = {}; previous = {}; return; end
l = {};
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline), l{end+1} = strtrim(tline(2:end)); end
end
fclose(fid);

e = find(cellfun(@isempty,l));
if numel(e)~=2, l = {}; e = [0 0]; end
current  = l(e(1)+1:e(2)-1)';
previous = l(e(2)+1:end)';
