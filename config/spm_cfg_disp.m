function disp = spm_cfg_disp
% SPM Configuration file for Image Display
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_disp.m 6952 2016-11-25 16:03:13Z guillaume $


%-------------------------------------------------------------------------
% data Image to Display
%-------------------------------------------------------------------------
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Image to Display';
data.help    = {'Image to display.'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 1];

%--------------------------------------------------------------------------
% disp Display Image
%--------------------------------------------------------------------------
disp         = cfg_exbranch;
disp.tag     = 'disp';
disp.name    = 'Display Image';
disp.val     = {data};
disp.help    = {
    'Interactive display of the orthogonal sections from an image volume.'
    'Clicking the cursor on either of the three images moves the point around which the orthogonal sections are viewed.  The co-ordinates of the cursor are shown both in voxel co-ordinates and millimetres within some fixed framework. The intensity at that point in the image (sampled using the current interpolation scheme) is also given. The position of the cross-hairs can also be moved by specifying the co-ordinates in millimetres to which they should be moved.  Clicking on the horizontal bar above these boxes will move the cursor back to the origin  (analogous to setting the cross-hair position (in mm) to [0 0 0]).'
    ''
    'The images can be re-oriented by entering appropriate translations, rotations and zooms into the panel on the left.  The transformations can then be saved by hitting the "Reorient images..." button.  The transformations that were applied to the image are saved to the header information of the selected images.  The transformations are considered to be relative to any existing transformations that may be stored.  Note that the order that the transformations are applied in is the same as in spm_matrix.m.'
    ''
    'The "Reset..." button next to it is for setting the orientation of images back to transverse.  It retains the current voxel sizes, but sets the origin of the images to be the centre of the volumes and all rotations back to zero.'
    ''
    'The right panel shows miscellaneous information about the image. This includes:'
    '   Dimensions - the x, y and z dimensions of the image.'
    '   Datatype   - the computer representation of each voxel.'
    '   Intensity  - scale-factors and possibly a DC offset.'
    '   Miscellaneous other information about the image.'
    '   Vox size   - the distance (in mm) between the centres of neighbouring voxels.'
    '   Origin     - the voxel at the origin of the co-ordinate system'
    '   DIr Cos    - Direction cosines.  This is a widely used representation of the orientation of an image.'
    ''
    'There are also a few options for different resampling modes, zooms etc.  You can also flip between voxel space (as would be displayed by Analyze) or world space (the orientation that SPM considers the image to be in).  If you are re-orienting the images, make sure that world space is specified.  Blobs (from activation studies) can be superimposed on the images and the intensity windowing can also be changed.'
    ''
    'If you have put your images in the correct file format, then (possibly after specifying some rigid-body rotations):'
    '    The top-left image is coronal with the top (superior) of the head displayed at the top and the left shown on the left. This is as if the subject is viewed from behind.'
    '    The bottom-left image is axial with the front (anterior) of the head at the top and the left shown on the left. This is as if the subject is viewed from above.'
    '    The top-right image is sagittal with the front (anterior) of the head at the left and the top of the head shown at the top. This is as if the subject is viewed from the left.'
    '/*\begin{figure} \begin{center} \includegraphics[width=150mm]{images/disp1} \end{center} \caption{The Display routine. \label{disp1}}\end{figure} */'
    }';
disp.prog = @disp_image;


%==========================================================================
function disp_image(job)
spm_image('Init', job.data{1});
