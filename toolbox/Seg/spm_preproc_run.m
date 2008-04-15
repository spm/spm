function varargout = spm_preproc_run(job,arg)
% Segment a bunch of images
% FORMAT spm_preproc(job)
% job.channel(n).vols{m}
% job.channel(n).biasreg
% job.channel(n).biasfwhm
% job.channel(n).write
% job.tissue(k).tpm
% job.tissue(k).ngaus
% job.tissue(k).native
% job.tissue(k).warped
% job.warp.affreg
% job.warp.reg
% job.warp.samp
% job.warp.write
% job.warp.bb
% job.warp.vox
%
% See the user interface for a description of the fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_preproc_run.m 1424 2008-04-15 20:49:14Z john $

if nargin==1,
    run_job(job);
elseif strcmpi(arg,'check'),
    varargout{:} = check_job(job);
elseif strcmpi(arg,'vfiles'),
    varargout{:} = vfiles_job(job);
else
    error('Unknown argument ("%s").', arg);
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function run_job(job)

tpm    = strvcat(cat(1,job.tissue(:).tpm));
tpm    = spm_load_priors8(tpm);

nit = 1;

for iter=1:nit,
    if nit>1,
        % Sufficient statistics for possible generation of group-specific
        % template data.
        SS = zeros([size(tpm.dat{1}),numel(tpm.dat)],'single');
    end
    for subj=1:numel(job.channel(1).vols),
        images = '';
        for n=1:numel(job.channel),
            images = strvcat(images,job.channel(n).vols{subj});
        end
        obj.image    = spm_vol(images);
        spm_check_orientations(obj.image);

        obj.fudge    = 5;
        obj.biasreg  = cat(1,job.channel(:).biasreg);
        obj.biasfwhm = cat(1,job.channel(:).biasfwhm);
        obj.tpm      = tpm;
        obj.lkp      = [];
        for k=1:numel(job.tissue),
            obj.lkp = [obj.lkp ones(1,job.tissue(k).ngaus)*k];
        end;
        obj.reg      = job.warp.reg;
        obj.samp     = job.warp.samp;

        if iter==1,
            % Initial affine registration.
            obj.Affine  = eye(4);
            if ~isempty(job.warp.affreg),
                obj.Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge*4,tpm,obj.Affine,job.warp.affreg);
                obj.Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge,  tpm,obj.Affine,job.warp.affreg);
            end;
        else
            % Load results from previous iteration for use with next round of
            % iterations, with the new group-specific tissue probability map.
            [pth,nam] = fileparts(job.channel(1).vols{subj});
            res       = load(fullfile(pth,[nam '_seg8.mat']));
            obj.Affine = res.Affine;
            obj.Twarp  = res.Twarp;
            obj.Tbias  = res.Tbias;
            obj.mg     = res.mg;
            obj.mn     = res.mn;
            obj.vr     = res.vr;
        end

        res = spm_preproc8(obj);

        try,
            [pth,nam] = fileparts(job.channel(1).vols{subj});
            savefields(fullfile(pth,[nam '_seg8.mat']),res);
        catch
        end

        if iter==nit,
            % Final iteration, so write out the required data.
            tmp2 =  cat(1,job.channel(:).write);
            tmp1 = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
            tmp2 =  cat(1,job.channel(:).write);
            tmp3 = job.warp.write;
            spm_preproc_write8(res,tmp1,tmp2,tmp3);
        else
            % Not the final iteration, so compute sufficient statistics for
            % re-estimating the template data.
            N    = numel(job.channel);
            K    = numel(job.tissue);
            cls  = spm_preproc_write8(res,zeros(K,4),zeros(N,2),[0 0]);
            for k=1:K,
                SS(:,:,:,k) = SS(:,:,:,k) + cls{k};
            end
        end

    end
    if iter<nit && nit>1,
         % Treat the tissue probability maps as Dirichlet priors, and compute the 
         % MAP estimate of group tissue probability map using the sufficient
         % statistics.
         alpha = 1.0;
         for k=1:K,
             SS(:,:,:,k) = SS(:,:,:,k) + spm_bsplinc(tpm.V(k),[0 0 0  0 0 0])*alpha + eps;
         end
         s = sum(SS,4);
         for k=1:K,
             tmp        = SS(:,:,:,k)./s;
             tpm.bg(k)  = mean(mean(tmp(:,:,1)));
             tpm.dat{k} = spm_bsplinc(log(tmp+tpm.tiny),[ones(1,3)*(tpm.deg-1)  0 0 0]);
         end
    end
end
return
%_______________________________________________________________________

%_______________________________________________________________________
function msg = check_job(job)
msg = {};
return
%_______________________________________________________________________

%_______________________________________________________________________
function vf = vfiles_job(job)
vf = {};
return
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if spm_matlab_version_chk('7') >= 0
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;

return;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function vf = vfiles(job)
opts  = job.output;
sopts = [opts.GM;opts.WM;opts.CSF];
vf    = cell(numel(job.data),2);
for i=1:numel(job.data),
    [pth,nam,ext,num] = spm_fileparts(job.data{i});
    vf{i,1} = fullfile(pth,[nam '_seg_sn.mat']);
    vf{i,2} = fullfile(pth,[nam '_seg_inv_sn.mat']);
    j       = 3;
    if opts.biascor,
        vf{i,j} = fullfile(pth,['m' nam ext ',1']);
        j       = j + 1;
    end;
    for k1=1:3,
        if sopts(k1,3),
            vf{i,j} = fullfile(pth,[  'c', num2str(k1), nam, ext, ',1']);
            j       = j + 1;
        end;
        if sopts(k1,2),
            vf{i,j} = fullfile(pth,[ 'wc', num2str(k1), nam, ext, ',1']);
            j       = j + 1;
        end;
        if sopts(k1,1),
            vf{i,j} = fullfile(pth,['mwc', num2str(k1), nam, ext, ',1']);
            j       = j + 1;
        end;
    end;
end;
vf = vf(:);

