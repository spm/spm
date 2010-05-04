function [H,Hnames] = spm_design_flexible(fblock,I)
% Create H partition of design matrix
% FORMAT [H,Hnames] = spm_design_flexible(fblock,I)
%
% fblock   - Part of job structure containing within-subject design info
% I        - Nscan x 4 factor matrix
%
% H        - Component of design matrix describing conditions
% Hnames   - Condition names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_flexible.m 3860 2010-05-04 15:59:25Z guillaume $

% Sort main effects and interactions
%--------------------------------------------------------------------------
fmain = struct('fnum',{});
inter = struct('fnums',{});
for k=1:numel(fblock.maininters)
    if isfield(fblock.maininters{k},'fmain')
        fmain(end+1)=fblock.maininters{k}.fmain;
    elseif isfield(fblock.maininters{k},'inter')
        inter(end+1)=fblock.maininters{k}.inter;
    end
end

% Create main effects
%--------------------------------------------------------------------------
H=[];Hnames=[];
nmain=length(fmain);
for f=1:nmain
    fcol=fmain(f).fnum;
    fname=fblock.fac(fcol).name;
    
    % Augment H partition - explicit factor numbers are 1 lower than in I matrix
    [Hf,Hfnames]=spm_DesMtx(I(:,fcol+1),'-',fname);
    H=[H,Hf];
    Hnames=[Hnames;Hfnames];
end

% Create interactions
%--------------------------------------------------------------------------
ni=length(inter);
for i=1:ni
    % Get the two factors for this interaction
    fnums=inter(i).fnums;
    f1=fnums(1);f2=fnums(2);
    
    % Names
    iname{1}     = fblock.fac(f1).name;
    iname{2}     = fblock.fac(f2).name;
    
    % Augment H partition - explicit factor numbers are 1 lower than in I matrix
    Isub         = [I(:,f1+1),I(:,f2+1)];
    [Hf,Hfnames] = spm_DesMtx(Isub,'-',iname);
    H            = [H,Hf];
    Hnames       = [Hnames;Hfnames];
    
end

if nmain==0 && ni==0
    error('You have not specified any main effects or interactions');
end
