function [con] = spm_design_contrasts (SPM)
% Make contrasts for one/two-way ANOVAs
% FORMAT [con] = spm_design_contrasts (SPM)
%
% SPM   SPM data structure
%
% This function generates contrasts on the basis of the current SPM 
% design. This is specified in SPM.factor (how the factors relate to the
% conditions) and SPM.xBF.order (how many basis functions per condition). 
%
% This function generates (transposed) contrast matrices to test:
% (1) average effect of condition, (2) main effect of factor 1
% (3) main effect of factor 2, (4) interaction
%
% con(c).c      Contrast matrix
%       .name   Name
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$

if length(SPM.Sess) > 1
    disp('Warning: spm_design_contrasts only works for single session data');
end

nf=length(SPM.factor);
k1=SPM.factor(1).levels;
if nf==2
    k2=SPM.factor(2).levels;
else
    k2=1;
end

icon=spm_make_contrasts(k1,k2);

if nf==1
    icon(3)=[];
    icon(4)=[];
end

% Get number of basis functions per condition
nbases=SPM.xBF.order;
% Number of regressors in design
k=size(SPM.xX.X,2);

for c=1:length(icon),
    z=kron(icon(c).c,eye(nbases));
    [zr,zc]=size(z);
    con(c).c=zeros(zr,k);
    con(c).c(1:zr,1:zc)=z;
end

con(1).name='Average effect of condition';
con(2).name=['Main effect of ',SPM.factor(1).name];
if nf==2
    con(3).name=['Main effect of ',SPM.factor(2).name];
    con(4).name=['Interaction between ',SPM.factor(1).name,' and ',SPM.factor(2).name];
end