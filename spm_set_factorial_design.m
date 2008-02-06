function [I,P,H,Hnames] = spm_set_factorial_design (job)
% Extract factorial matrix, file list and H partition of design matrix
% FORMAT [I,P,H,Hnames] = spm_set_factorial_design (job)
%
% job       job structure defined in spm_config_factorial_design
%
% I         Nscan x 4 factor matrix
% P         List of scans
% H         Component of design matrix describing conditions
% Hnames    Condition names
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_set_factorial_design.m 1131 2008-02-06 11:17:09Z spm $

% Get number of factors, names and levels
sF{1}='Repl'; % - first `factor' (for I) is always replication
Nfactors=length(job.des.fd.fact);
for i=1:Nfactors,
    % For creating design matrix and nonsphericity matrix
    sF{i+1}=job.des.fd.fact(i).name;
end

% Get number of scans in each cell, get names and create file list
P=[];
Ncells=length(job.des.fd.icell);  % Assume all cells defined for now !
switch Nfactors,
    case 1,
        % One-way designs
        cc=1;
        for i=1:job.des.fd.fact(1).levels,
            % Find corresponding cell and augment file list
            cell_num=0;
            for c=1:Ncells,
                levs=job.des.fd.icell(c).levels;
                if i==levs(1) 
                    cell_num=c;
                end
            end
            if cell_num==0
                disp(sprintf('Error: cell %d not defined',i));
            else
                P=[P;job.des.fd.icell(cell_num).scans];
                ns(cc)=length(job.des.fd.icell(cell_num).scans);
                cc=cc+1;
            end
        end
        
    case 2,
        % Two-way designs
        cc=1;
        for i=1:job.des.fd.fact(1).levels,
            for j=1:job.des.fd.fact(2).levels,
                % Find corresponding cell and augment file list
                cell_num=0;
                for c=1:Ncells,
                    levs=job.des.fd.icell(c).levels;
                    if i==levs(1) & j==levs(2)
                        cell_num=c;
                    end
                end
                if cell_num==0
                    disp(sprintf('Error: cell %d %d not defined',i,j));
                else
                    P=[P;job.des.fd.icell(cell_num).scans];
                    ns(cc)=length(job.des.fd.icell(cell_num).scans);
                    cc=cc+1;
                end
            end
        end
    case 3,
        % Three-way designs
        cc=1;
        for i=1:job.des.fd.fact(1).levels,
            for j=1:job.des.fd.fact(2).levels,
                for k=1:job.des.fd.fact(3).levels,
                    % Find corresponding cell and augment file list
                    cell_num=0;
                    for c=1:Ncells,
                        levs=job.des.fd.icell(c).levels;
                        if i==levs(1) & j==levs(2) & k==levs(3)
                            cell_num=c;
                        end
                    end
                    if cell_num==0
                        disp(sprintf('Error: cell %d %d %d not defined',i,j,k));
                    else
                        P=[P;job.des.fd.icell(cell_num).scans];
                        ns(cc)=length(job.des.fd.icell(cell_num).scans);
                        cc=cc+1;
                    end
                end
            end
        end
    otherwise,
        disp(sprintf('%d-way full-factoral designs are not supported',Nfactors));
end

% Create scan x factor matrix I
I=[];
switch Nfactors,
    case 1,
        % One-way designs
        c=1;
        for i=1:job.des.fd.fact(1).levels,
            I=[I;[1:ns(c)]',repmat(i,[ns(c) 1]),ones(ns(c),2)];
            c=c+1;
        end
        % Create H partition of design matrix
        [H,Hnames]=spm_DesMtx(I(:,2),'-',sF{2});
    case 2,
        % Two-way designs
        c=1;
        for i=1:job.des.fd.fact(1).levels,
            for j=1:job.des.fd.fact(2).levels,
                I=[I;[1:ns(c)]',repmat([i j],[ns(c) 1]),ones(ns(c),1)];
                c=c+1;
            end
        end
        % Create H partition of design matrix
        F2names={sF{2},sF{3}};
        [H,Hnames]=spm_DesMtx(I(:,2:3),'-',F2names);
    case 3,
        % Three-way designs
        c=1;
        for i=1:job.des.fd.fact(1).levels,
            for j=1:job.des.fd.fact(2).levels,
                for k=1:job.des.fd.fact(3).levels,
                    I=[I;[1:ns(c)]',repmat([i j k],[ns(c) 1])];
                    c=c+1;
                end
            end
        end
        % Create H partition of design matrix
        F2names={sF{2},sF{3},sF{4}};
        [H,Hnames]=spm_DesMtx(I(:,2:4),'-',F2names);
     otherwise,
        disp(sprintf('%d-way full-factoral designs are not supported',Nfactors));
end
