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

% John Ashburner
% $Id: spm_preproc_run.m 1151 2008-02-14 17:36:47Z john $

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

function run_job(job)

tpm    = strvcat(cat(1,job.tissue(:).tpm));
tpm    = spm_load_priors8(tpm);

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

    % Initial affine registration.
    obj.Affine  = eye(4);
    if ~isempty(job.warp.affreg),
        obj.Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge*4,tpm,obj.Affine,job.warp.affreg);
        obj.Affine  = spm_maff8(obj.image(1),job.warp.samp,obj.fudge,  tpm,obj.Affine,job.warp.affreg);
    end;
    res = spm_preproc8(obj);
    savefields('results.mat',res);

    tmp1 =  cat(1,job.channel(:).write);
    tmp2 = [cat(1,job.tissue(:).native) cat(1,job.tissue(:).warped)];
    tmp3 = [tmp2(:,1) any(tmp2(:,2:end),2)];
    cls  = spm_preproc_write8(res,tmp1,tmp3);

end

return

function msg = check_job(job)
msg = [];
return

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

