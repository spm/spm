function h = nifti(varargin)
% Create a NIFTI-1 object
% _______________________________________________________________________
% $Id$

switch nargin
case 0,
    org = niftistruc;
    hdr = [];
    for i=1:length(org),
        hdr.(org(i).label) = feval(org(i).dtype.conv,org(i).def);
    end;
    h = struct('hdr',hdr,'dat',[],'extras',[]);
    h = class(h,'nifti');
case 1
    if ischar(varargin{1})
        if size(varargin{1},1)>1,
            h = nifti(cellstr(varargin{1}));
            return;
        end;
        fname = varargin{1};
        vol   = read_hdr(fname);

        % Over-ride sform if a .mat file exists
        extras = read_extras(fname);
        %if isfield(extras,'mat') && size(extras.mat,3)>=1,
        %    mat            = extras.mat(:,:,1)*[eye(4,3) [1 1 1 1]'];
        %    vol.hdr.srow_x = mat(1,:);
        %    vol.hdr.srow_y = mat(2,:);
        %    vol.hdr.srow_z = mat(3,:);
        %    if vol.hdr.sform_code == 0, vol.hdr.sform_code = 2; end;
        %    if vol.hdr.sform_code == vol.hdr.qform_code,
        %        vol.hdr = encode_qform(extras.mat(:,:,1),vol.hdr);
        %    end;
        %end;

        dim   = double(vol.hdr.dim);
        dim   = dim(2:(dim(1)+1));
        dt    = double(vol.hdr.datatype);
        offs  = double(vol.hdr.vox_offset);

        if isfield(vol.hdr,'magic'),
            if ~vol.hdr.scl_slope && ~vol.hdr.scl_inter,
                vol.hdr.scl_slope = 1;
            end;
            slope = double(vol.hdr.scl_slope);
            inter = double(vol.hdr.scl_inter);
        else,
            if ~vol.hdr.roi_scale && ~vol.hdr.funused1,
                vol.hdr.roi_scale = 1;
            end;
            slope = double(vol.hdr.roi_scale);
            inter = double(vol.hdr.funused1);
        end;

        dat   = file_array(vol.iname,dim,[dt,vol.be],offs,slope,inter);
        h     = struct('hdr',vol.hdr,'dat',dat,'extras',extras);
        h     = class(h,'nifti');
    elseif isstruct(varargin{1})
        h     = class(varargin{1},'nifti');
    elseif iscell(varargin{1})
        fnames = varargin{1};
        h(numel(fnames)) = struct('hdr',[],'dat',[],'extras',[]);
        h     = class(h,'nifti');
        for i=1:numel(fnames),
            h(i) = nifti(fnames{i});
        end;
    else
        error('Dont know what to do yet.');
    end;
otherwise
    error('Dont know what to do yet');
end;
return;
