function [P,g,prior] = spm_bms_anova_img (P,g,prior)
% Log Bayes Factor against null for ANOVA; functional imaging data
% FORMAT [P,g,prior] = spm_bms_anova_img (P,g,prior)
%
% P         Cell array of filenames eg from SPM.xY.P with N cells
% g         [N x 1] vector with entries 1,2,3 etc denoting group membership
% prior     Specification of a single group is equivalent to a one sample t-test.
%           For this case you can specify 'unit' or 'jzs' (default) priors
%           See spm_bms_ttest.m and spm_bms_anova.m for more details
%__________________________________________________________________________
% Copyright (C) 2014-2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_bms_anova_img.m 6654 2015-12-22 12:55:36Z spm $

% Select files and groups if not provided
if nargin < 2 || isempty(P) || isempty(g)
    Ng=input('Enter number of groups ');
    P=[];g=[];
    for j=1:Ng,
        str=sprintf('Select images for group %d',j);
        Pj=cellstr(spm_select(Inf,'image',str));
        nj=length(Pj);
        P=[P;Pj];
        g=[g;j*ones(nj,1)];
    end
else
    Ng=length(unique(g));
end

if nargin < 3
    prior = 'jzs';
end

VY      = spm_data_hdr_read(P);
M       = VY{1}.mat;
DIM     = VY{1}.dim(1:3);
YNaNrep = spm_type(VY{1}.dt(1),'nanrep');
mask    = true(DIM);

fn = ['logBF_alt_',prior];
V = struct(...
    'fname',   [fn spm_file_ext],...
    'dim',     DIM,...
    'dt',      [spm_type('float32') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_bms_anova:LogBF against null');
V = spm_data_hdr_write(V);

nScan     = length(VY);
chunksize = floor(spm_get_defaults('stats.maxmem') / 8 / nScan);
nbchunks  = ceil(prod(DIM) / chunksize);
chunks    = min(cumsum([1 repmat(chunksize,1,nbchunks)]),prod(DIM)+1);

spm_progress_bar('Init',nbchunks,'Parameters estimation','Chunks');

for i=1:nbchunks
    chunk = chunks(i):chunks(i+1)-1;
    
    %-Report progress
    %======================================================================
    if i > 1, fprintf(repmat(sprintf('\b'),1,72)); end                  %-# 
    fprintf('%-40s: %30s', sprintf('Chunk %3d/%-3d',i,nbchunks),...
                           '...processing');                            %-#
                       
    %-Get data & construct analysis mask
    %======================================================================
    Y     = zeros(nScan,numel(chunk));
    cmask = mask(chunk);
    for j=1:nScan
        if ~any(cmask), break, end                 %-Break if empty mask
        
        Y(j,cmask) = spm_data_read(VY{j},chunk(cmask));%-Read chunk of data
        
        cmask(cmask) = Y(j,cmask) > -Inf;      %-Threshold (& NaN) mask
        if ~YNaNrep        %-Use implicit mask
            cmask(cmask) = abs(Y(j,cmask)) > eps;
        end
    end
    cmask(cmask) = any(diff(Y(:,cmask),1));        %-Mask constant data
    Y            = Y(:,cmask);                     %-Data within mask
    
    Nvoxels=size(Y,2);
    for n=1:Nvoxels
        if Ng==1
            logBF(n)=spm_bms_ttest(Y(:,n),prior);
        else
            logBF(n)=spm_bms_anova(Y(:,n),g);
        end
        if rem(n,1000)==0
            fprintf('Voxel %d out of %d\n',n,Nvoxels);                  %-#
        end
    end
    
    %-Write output file
    %======================================================================
    c = NaN(numel(chunk),1);
    c(cmask) = logBF;
    V = spm_data_write(V, c, chunk); 
    
    %-Report progress
    %======================================================================
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');             %-#
    spm_progress_bar('Set',i);
end

fprintf('\n');                                                          %-#
spm_progress_bar('Clear');
