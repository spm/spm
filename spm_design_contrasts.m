function [con] = spm_design_contrasts (SPM)
% Make contrasts for one, two or three-way ANOVAs
% FORMAT [con] = spm_design_contrasts (SPM)
%
% SPM   SPM data structure
%
% This function generates contrasts on the basis of the current SPM 
% design. This is specified in SPM.factor (how the factors relate to the
% conditions) and SPM.xBF.order (how many basis functions per condition). 
%
% This function generates (transposed) contrast matrices to test
% for the average effect of condition, main effects of factors and
% interactions
%
% con(c).c      Contrast matrix
%       .name   Name
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id: spm_design_contrasts.m 184 2005-05-31 13:23:32Z john $

if isempty(SPM.factor)
    % Can't create contrasts if factorial design has not been specified
    con=[];
    return;
end

nf=length(SPM.factor);
kf=[];
for f=1:nf,
    kf=[kf SPM.factor(f).levels];
end

icon=spm_make_contrasts(kf);

% Get number of basis functions per condition
nbases=SPM.xBF.order;

% Number of regressors in first session
%k=length(SPM.Sess(1).col);

for c=1:length(icon),
    con(c).c=kron(icon(c).c,eye(nbases));
end

con(1).name='Average effect of condition';
con(2).name=['Main effect of ',SPM.factor(1).name];
if nf>1
    con(3).name=['Main effect of ',SPM.factor(2).name];
    con(4).name=['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(2).name];
end
if nf>2
    con(5).name=['Main effect of ',SPM.factor(3).name];
    con(6).name=['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(3).name];
    con(7).name=['Interaction: ',SPM.factor(2).name,' x ',SPM.factor(3).name];
    con(8).name=['Interaction: ',SPM.factor(1).name,' x ',SPM.factor(2).name,' x ',SPM.factor(3).name];
end

% If there are multiple sessions, replicate each contrast matrix over
% sessions - this assumes the conditions are identical in each session
nsess=length(SPM.Sess);
ncon=length(con);
if nsess>1
    conds=length(SPM.Sess(1).U);
    for s=1:nsess,
        nconds=length(SPM.Sess(s).U);
        if ~(nconds==conds)
            disp('Error in spm_design_contrasts: number of conditions must be same in all sessions');
            return
        end
    end
    for c=1:ncon,
        con(c).c=kron(ones(1,nsess),con(c).c);
    end
end

% Pad contrasts out - add columns for session effects
k=size(SPM.xX.X,2);
for c=1:ncon,
    [zr,zc]=size(con(c).c);
    cmat=zeros(zr,k);
    cmat(1:zr,1:zc)=con(c).c;
    con(c).c=cmat;
end
