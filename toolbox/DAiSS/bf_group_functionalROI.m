function res = bf_group_functionalROI(BF, S)
% Computes Minimum Norm projectors
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, using the code from Matti Stenroos and Olaf Hauk
% http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools
%
% Please cite:
% Hauk O, Stenroos M.
% A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.
% Human Brain Mapping 2013
% $Id: bf_group_functionalROI.m 7706 2019-11-22 16:30:29Z spm $

%--------------------------------------------------------------------------
if nargin == 0
    measure = cfg_menu;
    measure.tag = 'measure';
    measure.name = 'Measure';
    measure.labels = {'masked covariance', 'variance'};
    measure.val = {'lJcov'};
    measure.values = {'lJcov', 'var'};
    measure.help = {'How to estimate measure'};
    
    spread = cfg_entry;
    spread.tag = 'spread';
    spread.name = 'Spread';
    spread.strtype = 'r';
    spread.num = [1 1];
    spread.val = {4};
    spread.help = {'bla bla bla'};
    
    threshold = cfg_entry;
    threshold.tag = 'threshold';
    threshold.name = 'Threshold';
    threshold.strtype = 'r';
    threshold.num = [1 1];
    threshold.val = {0.01};
    threshold.help = {'bla bla bla'};
    
    mincorr = cfg_entry;
    mincorr.tag = 'mincorr';
    mincorr.name = 'Minimal correlation inside basin';
    mincorr.strtype = 'r';
    mincorr.num = [1 1];
    mincorr.val = {0.7};
    mincorr.help = {'bla bla bla'};
    
    maxsize = cfg_entry;
    maxsize.tag = 'maxsize';
    maxsize.name = 'Maximal size of basin';
    maxsize.strtype = 'r';
    maxsize.num = [1 1];
    maxsize.val = {50};
    maxsize.help = {'bla bla bla'};
 
    distratio1 = cfg_entry;
    distratio1.tag = 'distratio1';
    distratio1.name = 'First ratio between correlation and geodesic based ditances';
    distratio1.strtype = 'r';
    distratio1.num = [1 1];
    distratio1.val = {4};
    distratio1.help = {'bla bla bla'};
    
    distratio2 = cfg_entry;
    distratio2.tag = 'distratio2';
    distratio2.name = 'Second ratio between correlation and geodesic based ditances';
    distratio2.strtype = 'r';
    distratio2.num = [1 1];
    distratio2.val = {2};
    distratio2.help = {'bla bla bla'};
    
    maxclustsize = cfg_entry;
    maxclustsize.tag = 'maxclustsize';
    maxclustsize.name = 'Maxclust size';
    maxclustsize.strtype = 'r';
    maxclustsize.num = [1 1];
    maxclustsize.val = {30};
    maxclustsize.help = {'bla bla bla'};
    
    maxclust = cfg_branch;
    maxclust.tag = 'maxclust';
    maxclust.name = 'Maxclust';
    maxclust.val  = {maxclustsize};

    cutoffthresh = cfg_entry;
    cutoffthresh.tag = 'cutoffthresh';
    cutoffthresh.name = 'Cutoff threshold';
    cutoffthresh.strtype = 'r';
    cutoffthresh.num = [1 1];
    cutoffthresh.val = {7};
    cutoffthresh.help = {'bla bla bla'};
   
    cutoff = cfg_branch;
    cutoff.tag = 'cutoff';
    cutoff.name = 'Cutoff';
    cutoff.val  = {cutoffthresh};
    
    cluster = cfg_choice;
    cluster.tag = 'cluster';
    cluster.name = 'Cluster method';
    cluster.values = {maxclust, cutoff};
    cluster.val = {cutoff};
    
    linkmeth = cfg_menu;
    linkmeth.tag = 'linkmeth';
    linkmeth.name = 'Linkage method';
    linkmeth.labels = {'complete', 'average'};
    linkmeth.val = {'complete'};
    linkmeth.values = {'complete', 'average'};
    linkmeth.help = {'bla bla bla'};

    similarity = cfg_entry;
    similarity.tag = 'similarity';
    similarity.name = 'Similarity threshold';
    similarity.strtype = 'r';
    similarity.num = [1 1];
    similarity.val = {0};
    similarity.help = {'bla bla bla'};

    functionalROI      = cfg_branch;
    functionalROI.tag  = 'functionalROI';
    functionalROI.name = 'Functional ROI';
    functionalROI.val  = {measure, spread, threshold, mincorr, maxsize, distratio1, distratio2, cluster, linkmeth, similarity};
    functionalROI.help = {'bla bla bla'};
    res = functionalROI;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

S.modality = 'MEGPLANAR'; %temp

Nl = length(BF);
S.Nl = Nl;

% make full data
J = [];
BFp = [];
dist_ex = 0;
for p=1:Nl
    
    BFp{p} = load(BF{p});
    
    if isfield(BFp{p}.inverse.(S.modality),'J')
        Nr = size(BFp{p}.inverse.(S.modality).S,2);
        Nt = size(BFp{p}.inverse.(S.modality).J,2)/Nr;
        Nd = size(BFp{p}.inverse.(S.modality).J,1);
        Jt = [];
        for i=1:Nt
            Jt(:,:,i) =  BFp{p}.inverse.(S.modality).J(:,1+(i-1)*Nr:Nr+(i-1)*Nr);
        end
        J = cat(1,J,Jt);
        
%     elseif isfield(BFp.inverse.(modality),'W') for MNE and others

    else
        error('There is no inversion');
    end
    
    % check distance field at least for one of the subjects
    if isfield(BFp{p}.inverse.(S.modality),'distance')
        dist_ex = 1;
        dist_in = p;
    end
        
end

% copy or calculate distance if needed; 
if dist_ex==0
   mesh = [];
   Dt = spm_eeg_load(BFp{1}.data.D);
   mesh.Vertices  = Dt.inv{Dt.val}.mesh.tess_mni.vert;
   mesh.Faces  = Dt.inv{Dt.val}.mesh.tess_mni.face;
   BFp{1}.inverse.(S.modality).distance = GALA_calculate_distance(mesh);
   outdir = spm_file(BF{1}, 'fpath');
   bf_save_path(BFp{1},fullfile(outdir, 'BF.mat'));
   for p=2:Nl
       if ~isfield(BFp{p}.inverse.(S.modality),'distance')
            BFp{p}.inverse.(S.modality).distance = BFp{1}.inverse.(S.modality).distance ;
            outdir = spm_file(BF{p}, 'fpath');
            bf_save_path(BFp{p},fullfile(outdir, 'BF.mat'));
       end
   end
else
   for p=1:Nl
       if ~isfield(BFp{p}.inverse.(S.modality),'distance')
            BFp{p}.inverse.(S.modality).distance = BFp{dist_in}.inverse.(S.modality).distance ;
            outdir = spm_file(BF{p}, 'fpath');
            bf_save_path(BFp{p},fullfile(outdir, 'BF.mat'));
       end
   end
end

% make measure vector
measure = [];
J1 = squeeze(mean(J,3));

switch S.measure
    case 'lJcov'
        
        Dt = spm_eeg_load(BFp{1}.data.D);
        vert  = Dt.inv{Dt.val}.mesh.tess_mni.vert;
        face  = Dt.inv{Dt.val}.mesh.tess_mni.face;
        A = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
        A = A^S.spread;
        ssQ1 = kron(ones(Nl,Nl),spones(A));

        lJcov=zeros(1,Nd*Nl);
        for i=1:Nd*Nl
            if any(J1(i,:))
                Jcovi = J1*J1(i,:)';
                sJcovi = Jcovi.*ssQ1(:,i);
                lJcov(i) = squeeze(sum(sJcovi));
            end
        end
        
    case 'var'
        
        lJcov = var(J1,2);
end

distance = BFp{1}.inverse.(S.modality).distance;

res = GALA_clustering(lJcov,J1, S, distance, A);
% load('G:\4VK\Processed\Sub02\MEEG\BF\res.mat');


Ncl = size(res.pclvi,2);
label = [];
tcourse = [];
Nr = size(BFp{1}.inverse.(S.modality).S,2); % Number of temporal modes
Nt = size(BFp{1}.inverse.(S.modality).J,2)/Nr; % Number of trials

for i=1:Ncl
    tvert = zeros(1,3);
    for p=1:Nl
        tvert = tvert+vert(res.pmaxs{p}(i),:);
        for j=1:Nt
            tcourse{p}(i,:,j) = mean(BFp{p}.inverse.(S.modality).J(res.pclvi{p,i},1+(j-1)*Nr:Nr+(j-1)*Nr),1)...
            *BFp{p}.inverse.(S.modality).S';
        end
    end
    label{i} = sprintf('%.2f_%.2f_%.2f', tvert/Nl);
end

res = cell(1, numel(BF));
for p=1:Nl
    ftdata.label = label(:);
    for j=1:Nt
        ftdata.trial{j} = tcourse{p}(:,:,j);
        D = BFp{p}.data.D;
        ftdata.time{j}  = D.time(1:D.nsamples);
    end
        
    BFp{p}.output.sourcedata.(S.modality).ftdata = ftdata;
    
    outdir = spm_file(BF{p}, 'fpath');
    bf_save_path(BFp{p},fullfile(outdir, 'BF.mat'));
    res(p) = {fullfile(outdir, 'BF.mat')};
end


end


