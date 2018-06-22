function V4 = spm_file_merge(V,fname,dt,RT)
% Concatenate 3D volumes into a single 4D volume
% FUNCTION V4 = spm_file_merge(V,fname,dt)
% V      - images to concatenate (char array or spm_vol struct)
% fname  - filename for output 4D volume [defaults: '4D.nii']
%          Unless explicit, output folder is the one containing first image
% dt     - datatype (see spm_type) [defaults: 0]
%          0 means same datatype than first input volume
% RT     - Interscan interval {seconds} [defaults: NaN]
%
% V4     - spm_vol struct of the 4D volume
%__________________________________________________________________________
%
% For integer datatypes, the file scale factor is chosen as to maximise
% the range of admissible values. This may lead to quantization error
% differences between the input and output images values.
%__________________________________________________________________________
% Copyright (C) 2009-2018 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_file_merge.m 7354 2018-06-22 10:44:22Z guillaume $

%-Input: V
%--------------------------------------------------------------------------
if ~nargin
    [V,sts] = spm_select([1 Inf],'image','Select images to concatenate');
    if ~sts, return; end
end
if ischar(V)
    V = spm_vol(V);
elseif iscellstr(V)
    V = spm_vol(char(V));
end
ind  = cat(1,V.n);
N    = cat(1,V.private);

%-Input: fname
%--------------------------------------------------------------------------
cwd = fileparts(V(1).fname);
if isempty(cwd), cwd = pwd; end

if nargin < 2
    fname = fullfile(cwd,'4D.nii');
else
    [p,n,e] = fileparts(fname);
    if isempty(e), fname = [fname '.nii']; end
    if isempty(p), fname = fullfile(cwd, fname); end
end

%-Input: dtype
%--------------------------------------------------------------------------
if nargin < 3
    dt = 0;
end
if dt == 0
    dt = V(1).dt(1);
end

%-Input: RT
%--------------------------------------------------------------------------
if nargin < 4
    RT = NaN;
end
if isnan(RT) && ...
   isfield(V(1).private,'timing') && isfield(V(1).private.timing,'tspace')
     RT = V(1).private.timing.tspace;
end

%-Set scalefactors and offsets
%==========================================================================
d = cat(1,V.dt); d = d(:,1);
s = cat(2,V.pinfo); s = s(1,:);
o = cat(2,V.pinfo); o = o(2,:);

%-Reuse parameters of input images if same scalefactor, offset and datatype
%--------------------------------------------------------------------------
if length(unique(s)) == 1 && length(unique(o)) == 1 ...
        && length(unique(d)) == 1 && d(1) == dt
    sf  = V(1).pinfo(1);
    off = V(1).pinfo(2);
else

    dmx  = spm_type(dt,'maxval');
    dmn  = spm_type(dt,'minval');

    %-Integer datatypes: scale to min/max of overall data
    %----------------------------------------------------------------------
    if isfinite(dmx)
        spm_progress_bar('Init',numel(V),'Computing scale factor','Volumes Complete');
        mx      = -Inf;
        mn      = Inf;
        for i=1:numel(V)
            dat = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
            dat = dat(isfinite(dat));
            mx  = max(mx,max(dat(:)));
            mn  = min(mn,min(dat(:)));
            spm_progress_bar('Set',i);
        end
        spm_progress_bar('Clear');
        if isempty(mx), mx = 0; end
        if isempty(mn), mn = 0; end

        if mx~=mn
            if dmn < 0
                sf = max(mx/dmx,-mn/dmn);
            else
                sf = mx/dmx;
            end
            off    = 0;
        else
            sf     = mx/dmx;
            off    = 0;
        end
    
    %-floating precison: no scaling
    %----------------------------------------------------------------------
    else
        sf         = 1;
        off        = 0;
    end
end

%-Create and write 4D volume image
%==========================================================================
spm_unlink(fname);

%-Create NifTI header
%--------------------------------------------------------------------------
ni         = nifti;
ni.dat     = file_array(fname,...
                        [V(1).dim numel(V)],...
                        [dt spm_platform('bigend')],...
                        0,...
                        sf,...
                        off);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = '4D image';
if ~isnan(RT)
    ni.timing = struct('toffset',0, 'tspace',RT);
end
create(ni);

%-Write 4D data
%--------------------------------------------------------------------------
spm_progress_bar('Init',size(ni.dat,4),'Saving 4D image','Volumes Complete');
for i=1:size(ni.dat,4)
    ni.dat(:,:,:,i) = N(i).dat(:,:,:,ind(i,1),ind(i,2));
    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

%-Fix ?form_code in header (mat_intent is changed by spm_get_space above)
%--------------------------------------------------------------------------
ni = nifti(fname);
ni.mat_intent  = N(1).mat_intent;
ni.mat0_intent = N(1).mat0_intent;
create(ni);

%-Remove .mat file if present and not necessary
%--------------------------------------------------------------------------
matfname = spm_file(fname,'ext','mat');
if spm_existfile(matfname)
    M = load(matfname);
    if isequal(fieldnames(M),{'mat'}) % contains only 'mat'
        if sum(sum(M.mat(:,:,1).^2))==0
            M.mat(:,:,1) = N(1).mat;
        end
        if sum(sum(diff(M.mat,1,3).^2))<1e-8
            spm_unlink(matfname);
        end
    end
end

%-Return spm_vol structure
%--------------------------------------------------------------------------
if nargout
    V4 = spm_vol(fname);
end
