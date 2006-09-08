function [status, Iimg, Limg] = spm_eeg_firstlevel_checkfile(filenames, Nt)
% function status = spm_eeg_firstlevel_checkfile(filename, Nt)
% Function that checks whether the time points (4th dimension) of image
% filename have length Nt. If Nt is 0, the function checks whether all files
% have the same length. The output status is Nt, if all files have the same
% length, or 0, if they don't.
%
% filename: matrix of filenames
% Nt: number of time points (0 if filenames are checked for consistency
% among themselves)
% status: 0, if length of responses not consistent, otherwise length of ith response
% Iimg: If status 0, index of (first) offending image.
% Limg: If status 0, length of offending response.

Iimg = 0; Limg = 0;

for i = 1:size(filenames, 1)
    
    V = nifti(deblank(filenames(i, :)));
    
    if length(V.dat.dim) < 4
        status = 0;
        return;
    end
    
    if Nt == 0
        Nt = V.dat.dim(4);
    end
    
    if Nt ~= V.dat.dim(4)
        status = 0;
        Iimg = i;
        Limg = length(V);
        return;
    end
end
        
status = Nt;


    