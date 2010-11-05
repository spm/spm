function [H,Hnames,B,Bnames] = spm_design_flexible(fblock,I)
% Create H partition of design matrix
% FORMAT [H,Hnames] = spm_design_flexible(fblock,I)
%
% fblock   - Part of job structure containing within-subject design info
% I        - Nscan x 4 factor matrix
%
% H        - Component of design matrix describing conditions
% Hnames   - Condition names
% B        - Component of design matrix describing blocks
% Bnames   - Block names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_flexible.m 4111 2010-11-05 12:59:10Z guillaume $

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

H = []; Hnames = [];
B = []; Bnames = [];

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
    % Get the two factors for this interaction
    f1 = inter(i).fnums(1);
    f2 = inter(i).fnums(2);
    
    % Names
    iname{1}     = fblock.fac(f1).name;
    iname{2}     = fblock.fac(f2).name;
    
    % Augment H partition - explicit factor numbers are 1 lower than in I matrix
    Isub         = [I(:,f1+1),I(:,f2+1)];
    [Hf,Hfnames] = spm_DesMtx(Isub,'-',iname);
    H            = [H,Hf];
    Hnames       = [Hnames; Hfnames];
end

if isempty(H)
    error('You have not specified any main effects or interactions');
end
