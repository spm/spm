function [I,P,H,Hnames] = spm_design_factorial(fd)
% Extract factorial matrix, file list and H partition of design matrix
% FORMAT [I,P,H,Hnames] = spm_design_factorial(fd)
%
% fd       - structure defined in spm_cfg_factorial_design
%            with fields fact and icell
%
% I        - Nscan x 4 factor matrix
% P        - List of scans
% H        - Component of design matrix describing conditions
% Hnames   - Condition names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny, Guillaume Flandin
% $Id: spm_design_factorial.m 3860 2010-05-04 15:59:25Z guillaume $

% Get number of factors, levels and cells
%--------------------------------------------------------------------------
Nfactors   = numel(fd.fact);
Nlevels    = [fd.fact.levels];
Ncells     = numel(fd.icell);

% Get file list and scan x factor matrix I
%--------------------------------------------------------------------------
I          = [];
P          = {};
for i=1:Ncells
    nc     = numel(fd.icell(i).scans);
    I      = [I; [[1:nc]', repmat(fd.icell(i).levels(:)',nc,1)]];
    P      = {P{:}, fd.icell(i).scans{:}}';
end
[I, A]     = sortrows(I,2:Nfactors+1);
P          = P(A);
I          = [I ones(size(I,1),4-size(I,2))];

% Create H partition of design matrix
%--------------------------------------------------------------------------
[H,Hnames] = spm_DesMtx(I(:,2:Nfactors+1),'-',{fd.fact(:).name});

% Display some warnings
%--------------------------------------------------------------------------
if Ncells ~= prod(Nlevels)
    disp('Some cell(s) are missing.');
end
if size(I,1) ~= size(unique(I,'rows'),1)
    disp('Some cell(s) are defined several times.');
end
