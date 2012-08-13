function [I,P,cov] = spm_design_within_subject(fblock,cov)
% Set up within-subject design when specified subject by subject
% FORMAT [I,P,cov] = spm_design_within_subject(fblock,cov)
%
% fblock   - Part of job structure containing within-subject design info
% cov      - Part of job structure containing covariate info
%
% I        - Nscan x 4 factor matrix
% P        - List of scans
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_within_subject.m 4837 2012-08-13 18:10:04Z guillaume $

%-Specify design subject-by-subject
%--------------------------------------------------------------------------
P      = {};
[I, subj] = deal([]);
for s  = 1:numel(fblock.fsuball.fsubject)
    P  = [P; fblock.fsuball.fsubject(s).scans];
    ns = length(fblock.fsuball.fsubject(s).scans);
    cc = fblock.fsuball.fsubject(s).conds;
    
    [ccr,ccc] = size(cc);
    if ccr~=ns && ccc~=ns
        error('Error for subject %d: conditions not specified for each scan',s);
    elseif ccr~=ccc && ccc==ns
        %warning('spm:transposingConditions',['Condition matrix ',...
        %    'appears to be transposed. Transposing back to fix.\n',...
        %    'Alert developers if it is not actually transposed.'])
        cc=cc';
    end
    subj=[subj;s*ones(ns,1)];
    % get real replications within each subject cell
    [unused cci  ccj] = unique(cc,'rows');
    repl = zeros(ns, 1);
    for k=1:max(ccj)
        repl(ccj==k) = 1:sum(ccj==k);
    end
    I = [I; [repl cc]];
end

%-Handle keyword factors 'repl' and 'subject'
%--------------------------------------------------------------------------
nf = length(fblock.fac);
subject_factor = false;
for i=1:nf
    if strcmpi(fblock.fac(i).name,'repl')
        % Copy 'replications' column to create explicit 'replications' factor
        I = [I(:,1:i) I(:,1) I(:,i+1:end)];
    end
    if strcmpi(fblock.fac(i).name,'subject')
        % Create explicit 'subject' factor
        I = [I(:,1:i) subj I(:,i+1:end)];
        subject_factor = true;
    end  
end

%-Re-order scans conditions and covariates into standard format
% This is to ensure compatibility with how variance components are created
%--------------------------------------------------------------------------
if subject_factor
    U=unique(I(:,2:nf+1),'rows');
    Un=length(U);
    Uc=zeros(Un,1);
    r=1;rj=[];
    for k=1:Un
        for j=1:size(I,1)
            match=sum(I(j,2:nf+1)==U(k,:))==nf;
            if match
                Uc(k)=Uc(k)+1;
                Ir(r,:)=[Uc(k),I(j,2:end)];
                r=r+1;
                rj=[rj;j];
            end
        end
    end
    P=P(rj); % -scans
    I=Ir;    % -conditions
    for k=1:numel(cov) % -covariates
        cov(k).c = cov(k).c(rj);
    end
end

%-Pad out factorial matrix to cover the four canonical factors
%--------------------------------------------------------------------------
I = [I ones(size(I,1),4-size(I,2))];
