function opts = spm_config_slice_timing
% configuration file for slice timing
%____________________________________________________________________
% D Gitelman 03/07/01
% DRG CS-RCS: $Id: spm_config_slice_timing.m,v 1.2 2005-02-08 21:03:45-06 drg Exp drg $

w = spm_jobman('HelpWidth');

% ---------------------------------------------------------------------
scans.type = 'files';
scans.name = 'Sessions';
scans.tag  = 'scans';
scans.filter = 'image';
scans.num  = [2 Inf];
scans.help = {'Select images to acquisition correct.'};
% ---------------------------------------------------------------------

data.type = 'repeat';
data.name = 'Data';
data.values = {scans};
data.help = spm_justify(w,...
    'Subjects or sessions. The same parameters specified below will',...
    'be applied to all sessions');

% ---------------------------------------------------------------------

nslices.type = 'entry';
nslices.name = 'Number of Slices';
nslices.tag  = 'nslices';
nslices.strtype = 'n';
nslices.num  = [1 1];
nslices.help = {'Enter the number of slices'};
% ---------------------------------------------------------------------

refslice.type = 'entry';
refslice.name = 'Reference Slice';
refslice.tag  = 'refslice';
refslice.strtype = 'n';
refslice.num  = [1 1];
refslice.help = {'Enter the reference slice'};
% ---------------------------------------------------------------------

TR.type = 'entry';
TR.name = 'TR';
TR.tag  = 'tr';
TR.strtype = 'r';
TR.num  = [1 1];
TR.help = {'Enter the TR in seconds'};
% ---------------------------------------------------------------------

TA.type = 'entry';
TA.name = 'TA';
TA.tag  = 'ta';
TA.strtype = 'e';
TA.num  = [1 1];
TA.help = spm_justify(w,'The TA (in secs) must be entered by the user. ',...
    'It is usually calculated as TR-(TR/nslices). You can simply enter',...
    'this equation with the variables replaced by appropriate numbers');

% ---------------------------------------------------------------------

sliceorder.type = 'entry';
sliceorder.name = 'Slice order';
sliceorder.tag = 'so';
sliceorder.strtype = 'e';
sliceorder.num = [1 Inf];
p1 = spm_justify(w,'Enter the slice order. Bottom slice = 1. Sequence types ',...
        'and examples of code to enter are given below.');
p2 = spm_justify(w,'ascending (first slice=bottom): [1:1:nslices]');
p3 = spm_justify(w,'descending (first slice=top): [nslices:-1:1]');
p4 = spm_justify(w,'interleaved (middle-top):  for k = 1:nslices, ',...
        'round((nslices-k)/2 + (rem((nslices-k),2) * (nslices - 1)/2)) + 1, end ');
p5 = spm_justify(w, 'interleaved (bottom -> up): [1:2:nslices 2:2:nslices]');
p6 = spm_justify(w, 'interleaved (top -> down): [nslices:-2:1, nslices-1:-2:1]');
sliceorder.help = {p1{:},'',p2{:},'',p3{:},'',p4{:},'',p5{:},'',p6{:}};
    
% ---------------------------------------------------------------------

opts.type = 'branch';
opts.name = 'Slice Timing';
opts.tag  = 'st';
opts.val  = {data,nslices,TR,TA,sliceorder,refslice};
opts.prog = @slicetiming;
opts.vfiles = @vfiles;
opts.modality = {'FMRI'};
opts.help = {'Correct differences in image acquisition time between slices. '};
% ---------------------------------------------------------------------

function slicetiming(varargin)
job = varargin{1};
Seq = job.so;
TR  = job.tr;
TA  = job.ta;
nslices   = job.nslices;
refslice  = job.refslice;
timing(2) = TR - TA;
timing(1) = TA / (nslices -1);

for i = 1:length(job.scans)
    P   = strvcat(job.scans{i});
    spm_slice_timing(P,Seq,refslice,timing)
end
return;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
function vf = vfiles(varargin)
job = varargin{1};
vf  = cell(numel([job.scans{:}]),1);
n = 1;
for i=1:numel(job.scans),
    for j = 1:numel(job.scans{i})
    [pth,nam,ext,num] = spm_fileparts(job.scans{i}{j});
    vf{n} = fullfile(pth,['a' nam '.img' num]);
    n = n+1;
    end
end;
return;
% ---------------------------------------------------------------------

% ---------------------------------------------------------------------
