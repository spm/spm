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
% $Id: spm_design_contrasts.m 901 2007-08-30 07:02:57Z volkmar $

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
if isfield(SPM,'xBF')
    nbases=SPM.xBF.order;
else
    % for 2nd level designs .xBF not defined
    nbases=1;
end
for c=1:length(icon),
    con(c).c=kron(icon(c).c,eye(nbases));
end

if isfield(SPM,'Sess')
    % For 2nd level designs .Sess won't be specified
    % and we don't need to worry about parametric mod/multi-sessions
    
    % Pad out parametric modulation columns with zeros
    for c=1:length(icon),
        conds=length(SPM.Sess(1).U);
        nr=size(con(c).c,1);
        col=1;
        block=[];
        for cc=1:conds,
            block=[block,con(c).c(:,col:col+nbases-1)];
            block=[block,zeros(nr,SPM.Sess(1).U(cc).P.h)];
            col=col+nbases;
        end
        con(c).c=block;
    end
    
    % If there are multiple sessions, replicate each contrast matrix over
    % sessions - this assumes the conditions are identical in each session
    nsess=length(SPM.Sess);
    ncon=length(con);
    if nsess>1
        conds=length(SPM.Sess(1).U);
	covs =zeros(1,nsess);
        for s=1:nsess,
            nconds=length(SPM.Sess(s).U);
            if ~(nconds==conds)
                disp('Error in spm_design_contrasts: number of conditions must be same in all sessions');
                return
            end
	    covs(s) = size(SPM.Sess(s).C.C,2);
        end
        for c=1:ncon,
	    c1 = [con(c).c zeros(1,covs(1))];
	    for s=2:nsess
		c1 = [c1 con(c).c zeros(1,covs(s))];
	    end;
	    con(c).c = c1;
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