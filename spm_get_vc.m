function SPM=spm_get_vc(SPM)
% generate error variance components for SPM design
% FORMAT SPM = spm_get_vc(SPM);
%
% in: SPM   - SPM struct
% required fields:
%    factor  - factor description
%    .levels   - # levels for given factor
%    .variance - 0: all levels have identical error variance
%                1: error variance different across levels
%    .dept     - 0: error variances are independent between levels
%                1: error variances are dependent between levels
%    xVi.I       - nscan-by-4 index matrix that associates each scan
%                  with levels for up to 3 user defined factors (column 1
%                  reserved for implicit replication factor)
%
% out: SPM.xVi.Vi   - cell vector of covariance components
%_______________________________________________________________________
% 
% spm_get_vc generates variance components for a given design. For each
% factor, the user specifies whether its levels have identical variances
% and are uncorrelated.
%_______________________________________________________________________
% Copyright (C) 2006 Freiburg Brain Imaging 

% Volkmar Glauche
% $Id: spm_get_vc.m 622 2006-09-13 06:04:27Z volkmar $

nscan = size(SPM.xVi.I,1);
depVi={};
VVi={};
Vi={};

% identical/non-identical variances
% for each factor, create a single variance component if variances are
% identical across levels, and level specific variance components if
% variances are non-identical
for f=1:numel(SPM.factor)
    if ~SPM.factor(f).variance
        varVi{f}{1} = speye(nscan,nscan);
    else
        for l=1:max(SPM.factor(f).levels)
            varVi{f}{l}=double(spdiags(SPM.xVi.I(:,f+1)==l,0,nscan,nscan));
        end;
    end;
end;

% split factor specific variance components if they contain non-identical
% components of another factor - this may produce duplicate or all-zero
% variance components
switch numel(SPM.factor),
case 1,
    VVi = varVi{1};
case 2,
    for l1=1:numel(varVi{1})
        for l2=1:numel(varVi{2})
            VVi{end+1}=varVi{1}{l1}*varVi{2}{l2};
        end;
    end;
case 3,
    for l1=1:numel(varVi{1})
        for l2=1:numel(varVi{2})
            for l3=1:numel(varVi{3})
                VVi{end+1}=varVi{1}{l1}*varVi{2}{l2}*varVi{3}{l3};
            end;
        end;
    end;
case 4,
    for l1=1:numel(varVi{1})
        for l2=1:numel(varVi{2})
            for l3=1:numel(varVi{3})
                for l4=numel(varVi{4})
                    VVi{end+1}=varVi{1}{l1}*varVi{2}{l2}* ...
                              varVi{3}{l3}*varVi{4}{l4};
                end;
            end;
        end;
    end;
end;

%independent/dependent variances
for f=1:numel(SPM.factor)
    if SPM.factor(f).dept
        tmpI = SPM.xVi.I;
        tmpI(:,f+1) = [];
        % find levels, where non-independent factor varies, but the
        % others do not. For each combination of the other factors,
        % create a temporary factor
        [unused unused tmplvls] = unique(tmpI,'rows');
        for l1=1:max(SPM.factor(f).levels)
            ind1 = find(SPM.xVi.I(:,f+1)==l1);
            for l2=1:(l1-1)
                ind2 = find(SPM.xVi.I(:,f+1)==l2);
                tmpdepVi = sparse([ind1; ind2],...
                                      [ind2; ind1],...
                                      1,nscan,nscan);
                for tmpl=1:max(tmplvls)
                    depVi{end+1}=sparse(nscan,nscan);
                    rows=find(tmplvls==tmpl);
                    depVi{end}(rows,rows)=tmpdepVi(rows,rows);
                end;
            end;
        end
    end;
end;

Vi=[VVi depVi];

% remove dupl & all-zero variance components
dupl = zeros(size(Vi));
for cVi = 1:numel(Vi)
    dupl(cVi) = dupl(cVi)||full(sum(Vi{cVi}(:)) == 0);
    if ~dupl(cVi)
        for cVi1 = (cVi+1):numel(Vi)
            dupl(cVi1) = dupl(cVi1)||full(all(Vi{cVi}(:) == Vi{cVi1}(:)));
        end
    end;
end;
Vi = Vi(~dupl);
if numel(Vi) == 1
    SPM.xVi.V = Vi{1};
else
    SPM.xVi.Vi = Vi;
end;