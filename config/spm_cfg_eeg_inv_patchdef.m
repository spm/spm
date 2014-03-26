function inv_patchdef = spm_cfg_eeg_inv_patchdef
% Configuration file for taking a number of previous inversion results (maybe based on different data), smoothing and creating an approximate posterior
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_patchdef.m 5927 2014-03-26 16:02:56Z gareth $

D = cfg_files;
D.tag = 'D';
D.name = 'Prior patch density file';
D.filter = 'gii';
D.num = [1 Inf];
D.help = {'Select the prior patch distribution file'};



npatch = cfg_entry;
npatch.tag = 'npatch';
npatch.name = 'Number of patches';
npatch.strtype = 'n';
npatch.help = {'Total number of patches'};
npatch.val = {512};

niter = cfg_entry;
niter.tag = 'niter';
niter.name = 'Number of iterations';
niter.strtype = 'n';
niter.help = {'Total number of random patch sets'};
niter.val = {32};



inv_patchdef = cfg_exbranch;
inv_patchdef.tag = 'inv_patchdef';
inv_patchdef.name = 'Set prior sampling';
inv_patchdef.val = {D,npatch,niter};
inv_patchdef.help = {'Use to make up a patch definition file'};
inv_patchdef.prog = @run_inv_patchdef;
inv_patchdef.vout = @vout_inv_patchdef;
inv_patchdef.modality = {'MEG'};

function  out = run_inv_patchdef(job)


inverse = [];
if numel(job.D)<1,
    error('Need to add a number of files to combine');
end;


prior_mesh=gifti(job.D);

pprob=prior_mesh.cdata(:);
[sprior,sortind]=sort(pprob);
cumprior=cumsum(sprior);



figure;
subplot(2,1,1);
plot(1:length(pprob),pprob,'o');


rng('shuffle');
Ip=zeros(job.niter,job.npatch);
dum=zeros(size(sprior));
for j=1:job.niter,
    for k=1:job.npatch,
         pval=(k-1)./job.npatch;
         sind=find(cumprior>=pval);
         randind=randperm(length(sind));
         Ip(j,k)=sortind(sind(randind(1)));
         dum(Ip(j,k))=dum(Ip(j,k))+1;
    end;
    subplot(2,1,2); hold on;
    plot(1:length(pprob),dum,'.');
    drawnow;
    
end;

[a1,b1,c1]=fileparts(job.D{1});
patchfilename=[a1 filesep b1 sprintf('N%d_Np%d_%d.mat',job.niter,job.npatch,randi(1e6))]
save(patchfilename,'Ip');


out.D = {patchfilename};

function dep = vout_inv_patchdef(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Smooth posterior from a number of imaging source reconstructions';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

