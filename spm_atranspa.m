
% Multiplies the transpose of a matrix by itself.
% FORMAT C = spm_atranspa(A)
% A - real matrix
% C - real symmetric matrix resulting from A'A
%___________________________________________________________________________
% This routine was written to save both memory and CPU time.
% The memory saving is achieved by not having to generate A'.
% CPU saving is by only generating half of C, and filling the
% rest in later.
%___________________________________________________________________________

% %W% John Ashburner MRCCU/FIL %E%
