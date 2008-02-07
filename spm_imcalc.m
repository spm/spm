function Vo = spm_imcalc(Vi,Vo,f,flags,varargin)
% Perform algebraic functions on images
% FORMAT Vo = spm_imcalc(Vi,Vo,f,flags,Xtra_vars...)
% Vi            - vector of mapped image volumes to work on (from spm_vol)
% Vo (input)    - handle structure containing information on output image
%                 ( pinfo field is computed for the resultant image data, )
%                 ( and can be omitted from Vo on input.  See spm_vol     )
% f             - expression to be evaluated
% flags         - cell vector of flags: {dmtx,mask,hold}
% dmtx          - Read images into data matrix?
%                 [defaults (missing or empty) to 0 - no]
% mask          - implicit zero mask?
%                 [defaults (missing or empty) to 0]
%                  ( negative value implies NaNs should be zeroed )
% hold          - interpolation hold (see spm_slice_vol)
%                 [defaults (missing or empty) to 0 - nearest neighbour]
% Xtra_vars...  - additional variables which can be used in expression
% Vo (output)   - handle structure of output image volume after modifications
%                 for writing
%
% With no arguments, spm_imcalc_ui is called.
%_______________________________________________________________________
%
% The images specified in Vi, are referred to as i1, i2, i3,...  in the
% expression to be evaluated, unless the dmtx flag is setm in which
% case the images are read into a data matrix X, with images in rows.
%
% Computation is plane by plane, so in data-matrix mode, X is a NxK
% matrix, where N is the number of input images [prod(size(Vi))], and K
% is the number of voxels per plane [prod(Vi(1).dim(1:2))].
%
% For data types without a representation of NaN, implicit zero masking
% assummes that all zero voxels are to be treated as missing, and
% treats them as NaN. NaN's are written as zero (by spm_write_plane),
% for data types without a representation of NaN.
%
% See spm_imcalc_ui for example usage...
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Andrew Holmes
% $Id: spm_imcalc.m 1143 2008-02-07 19:33:33Z spm $



%-Parameters & arguments
%=======================================================================
if nargin<3, error('insufficient arguments'), end
if nargin<4, flags={}; end

if length(flags)<3, hold=[]; else hold=flags{3}; end
if isempty(hold), hold=0; end
if length(flags)<2, mask=[]; else mask=flags{2}; end
if isempty(mask), mask=0; end
if length(flags)<1, dmtx=[]; else dmtx=flags{1}; end
if isempty(dmtx), dmtx=0; end


%-Process any additional variables
%-----------------------------------------------------------------------
if nargin>4
    reserved = {'Vi','Vo','f','flags','hold','mask','dmtx','varargin',...
            'reserved','n','Y','p','B','X','i','M','d','sf'};
    for i=5:nargin
        if any(strcmp(inputname(i),reserved))
            error(['additional parameter (',inputname(i),...
                ') clashes with internal variable'])
        end
        eval([inputname(i),' = varargin{i-4};'])
    end
end


%=======================================================================
%-Computation
%=======================================================================
n   = numel(Vi);                %-#images
if n==0, error('no input images specified'), end
Y   = zeros(Vo.dim(1:3));       %-result of calculations


%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',Vo.dim(3),f,'planes completed');


%-Loop over planes computing result Y
%-----------------------------------------------------------------------
for p = 1:Vo.dim(3),
    B = spm_matrix([0 0 -p 0 0 0 1 1 1]);

    if dmtx, X=zeros(n,prod(Vo.dim(1:2))); end
    for i = 1:n
        M = inv(B*inv(Vo.mat)*Vi(i).mat);
        d = spm_slice_vol(Vi(i),M,Vo.dim(1:2),[hold,NaN]);
        if (mask<0), d(isnan(d))=0; end;
        if (mask>0) && ~spm_type(Vi(i).dt(1),'nanrep'), d(d==0)=NaN; end
        if dmtx, X(i,:) = d(:)'; else eval(['i',num2str(i),'=d;']); end
    end

    eval(['Yp = ' f ';'],['error([''Can''''t evaluate "'',f,''".'']);']);
    if prod(Vo.dim(1:2)) ~= numel(Yp),
        error(['"',f,'" produced incompatible image.']); end
    if (mask<0), Yp(isnan(Yp))=0; end
    Y(:,:,p) = reshape(Yp,Vo.dim(1:2));

    spm_progress_bar('Set',p);
end


%-Write output image (uses spm_write_vol - which calls spm_write_plane)
%-----------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y);

%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear')
