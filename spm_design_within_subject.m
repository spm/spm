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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_design_within_subject.m 3860 2010-05-04 15:59:25Z guillaume $

nsub=length(fblock.fsuball.fsubject);
% Specify design subject-by-subject
P=[];I=[];
subj=[];
for s=1:nsub,
    P  = [P; fblock.fsuball.fsubject(s).scans];
    ns = length(fblock.fsuball.fsubject(s).scans);
    cc = fblock.fsuball.fsubject(s).conds;
    
    [ccr,ccc] = size(cc);
    if ~(ccr==ns) && ~(ccc==ns)
        disp(sprintf('Error for subject %d: conditions not specified for each scan',s));
        return
    elseif ~(ccr==ccc) && (ccc==ns)
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
    end;
    I = [I; [repl cc]];
end

nf=length(fblock.fac);
subject_factor=0;
for i=1:nf,
    if strcmpi(fblock.fac(i).name,'repl')
        % Copy `replications' column to create explicit `replications' factor
        nI=I(:,1:i);
        nI=[nI,I(:,1)];
        nI=[nI,I(:,i+1:end)];
        I=nI;
    end
    if strcmpi(fblock.fac(i).name,'subject')
        % Create explicit `subject' factor
        nI=I(:,1:i);
        nI=[nI,subj];
        nI=[nI,I(:,i+1:end)];
        I=nI;
        subject_factor=1;
    end
    
end

% Re-order scans conditions and covariates into standard format
% This is to ensure compatibility with how variance components are created
if subject_factor
    U=unique(I(:,2:nf+1),'rows');
    Un=length(U);
    Uc=zeros(Un,1);
    r=1;rj=[];
    for k=1:Un,
        for j=1:size(I,1),
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
    end;
end

% Pad out factorial matrix to cover the four canonical factors
[ns,nI]=size(I);
if nI < 4
    I = [I, ones(ns,4-nI)];
end
