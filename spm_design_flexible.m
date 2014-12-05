function [H,Hnames,B,Bnames] = spm_design_flexible(fblock,I)
% Create H partition of design matrix
% FORMAT [H,Hnames,B,Bnames] = spm_design_flexible(fblock,I)
%
% fblock   - Part of job structure containing within-subject design info
% I        - Nscan x 4 factor matrix
%
% H        - Component of design matrix describing conditions
% Hnames   - Condition names
% B        - Component of design matrix describing blocks
% Bnames   - Block names
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_flexible.m 6280 2014-12-05 14:52:06Z guillaume $

%-Sort main effects and interactions
%--------------------------------------------------------------------------
fmain = struct('fnum',{});
inter = struct('fnums',{});
for k=1:numel(fblock.maininters)
    if isfield(fblock.maininters{k},'fmain')
        fmain(end+1) = fblock.maininters{k}.fmain;
    elseif isfield(fblock.maininters{k},'inter')
        inter(end+1) = fblock.maininters{k}.inter;
    end
end

H = []; Hnames = {};
B = []; Bnames = {};

%-Create main effects
%--------------------------------------------------------------------------
for f=1:length(fmain)
    fcol  = fmain(f).fnum;
    fname = fblock.fac(fcol).name;
    
    % Augment H partition - explicit factor numbers are 1 lower than in I matrix
    [Hf,Hfnames] = spm_DesMtx(I(:,fcol+1),'-',fname);
    if strcmpi(fname,'subject')
        B        = [B,Hf];
        Bnames   = [Bnames; Hfnames];
    else
        H        = [H,Hf];
        Hnames   = [Hnames; Hfnames];
    end
end

%-Create interactions
%--------------------------------------------------------------------------
for i=1:length(inter)
    % Factors for this interaction
    f            = inter(i).fnums;
    
    % Names
    iname        = {fblock.fac(f).name};
    
    % Augment H partition - explicit factor numbers are 1 lower than in I matrix
    [Hf,Hfnames] = spm_DesMtx(I(:,f+1),'-',iname);
    H            = [H,Hf];
    Hnames       = [Hnames; Hfnames];
end
