function [flist] = spm_eeg_firstlevel_searchdir(Dname, flist)
% searches through directory recursively and adds mat-files to a list.
% the mat-files must fullfill certain criteria.

D = dir(Dname);
p = Dname;

for i = 1:length(D)
    if ~D(i).isdir
       % add if matfile
       [tmp1, tmp2, ext] = fileparts(D(i).name);
       if strcmp(ext, '.img')
           flist = strvcat(flist, fullfile(p, D(i).name));
       end
    else
        if ~strcmp(D(i).name, '.') & ~strcmp(D(i).name, '..')
            flist = spm_eeg_firstlevel_searchdir(fullfile(p, D(i).name), flist);
        end
    end    
end