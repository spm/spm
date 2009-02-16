function [] = spm_run_bms_vis(varargin)
% Show results for Bayesian Model Selection Maps
% SPM job execution function
% takes a harvested job data structure (or no input) and calls SPM 
% functions to show results from Bayesian Model Selection of 
% Log. Evidence Maps  
%
% Input:
% Varargin - can be harvested job data structure (see matlabbatch help).
% Output:
% No output.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Maria Joao Rosa
% $Id: spm_run_bms_vis.m 2751 2009-02-16 15:50:26Z maria $

% Input
% -------------------------------------------------------------------------
if isempty(varargin)
   job.file{1}=''; job.img{1}=''; job.thres = []; job.scale = [];
else
   job = varargin{1};
end
file_str = length(job.file{1});
image    = length(job.img{1});
thres_b  = length(job.thres);
odds_b   = length(job.scale);

% Initialise SPM Interactive window
% -------------------------------------------------------------------------
[Finter,F,CmdLine] = spm('FnUIsetup','SPM BMC Toolbox: Visualise options',0);
F = spm_figure('FindWin','Graphics');
spm('Clear',Finter);
WS = spm('WinScale');
FS = spm('FontSizes');

% Select results (.img) from BMS Maps
% -------------------------------------------------------------------------
if file_str
   file_fname = job.file{1};
else
   file_fname = spm_select(1,'mat','Load BMS file (.mat)');
end

% Load BMS
load(file_fname);

% Select results (.img) from BMS Maps
% -------------------------------------------------------------------------
if image
   post = job.img{1};
else
   post = spm_select(1,'image','Select BMS Maps results (ex: ppm, alpha or epm image)');
end

% Select threshold to apply to image
% -------------------------------------------------------------------------
if thres_b
   threshold = job.thres;
else
   threshold = spm_input('Probability Threshold:', '+1', 'r', '0.5', [1, inf]);
end

% Select scale
% -------------------------------------------------------------------------
if odds_b
   odds = job.scale;
else
   odds = spm_input('Change scale to Log-Odds?', '+1', 'y/n', [1,0], 2);
end

% Image dimensions
% -------------------------------------------------------------------------
V             = spm_vol(post);
M             = V.mat;
iM            = inv(M);
vox           = sqrt(sum(M(1:3,1:3).^2));
DIM           = V.dim(1:3)'; 
xdim          = DIM(1); ydim  = DIM(2); zdim  = DIM(3);
[xords,yords] = ndgrid(1:xdim,1:ydim);
xords         = xords(:)';  yords = yords(:)';
I             = 1:xdim*ydim;
zords_init    = repmat(1,1,xdim*ydim);

% Legend of results being plotted
% -------------------------------------------------------------------------
b_str        = find(V.fname == '_');
separators   = find(V.fname == '\');
if length(separators) == 0
   separators=find(V.fname == '/');
end
res_name     = V.fname(b_str(end)+1);
res_model    = V.fname(b_str(end)-1);
subset_model = V.fname(separators(end)+1:b_str(end)-2);

% Show results being displayed on graphics window
switch res_name
       case 'a'        
            xSPM.str = sprintf('Alphas: %s %s',subset_model,res_model);
       case 'p'
            xSPM.str = sprintf('PPM: %s %s',subset_model,res_model);
       case 'x'
            xSPM.str = sprintf('PPM: %s %s',subset_model,res_model);
       case 'e'
            xSPM.str = sprintf('EPM: %s %s',subset_model,res_model);
       otherwise
            xSPM.str = '';
end

% Loop over slices (get voxels above threshold)
% -------------------------------------------------------------------------
xyz_above = [];
z_above   = [];

for z = 1:zdim,
    
    zords = z*zords_init;
    xyz   = [xords(I); yords(I); zords(I)];
    zvals = spm_get_data(V,xyz);
    above = find(zvals>threshold);
    if ~isempty(above)
        xyz_above = [xyz_above,xyz(:,above)];
        z_above   = [z_above,zvals(above)];
    end

end

% SPM figure: SPM, MIP and GUI
% -------------------------------------------------------------------------
if ~isempty(z_above)
    
    z_ps = z_above;
    % Log odds transform
    if odds
       z_odds = log(z_above./(ones(1,length(z_above))-z_above));
       if isempty(find(z_odds == Inf))  % Data don't contain Infs
          z_above = z_odds;
        else
          odds = 0;                     % Don't plot odds
        end
    end

    % Save data in xSPM structure
    voxels_mm     = M*[xyz_above;ones(1,size(xyz_above,2))];
    voxels_mm     = voxels_mm(1:3,:);
    xSPM.swd      = file_fname;
    xSPM.STAT     = '';
    xSPM.Z        = z_above;
    xSPM.XYZ      = xyz_above;
    xSPM.XYZmm    = voxels_mm;
    xSPM.xVol.M   = M;
    xSPM.xVol.DIM = DIM;
    xSPM.DIM      = DIM;
    xSPM.M        = M;
    xSPM.iM       = iM;
    xSPM.n        = 1;
    xSPM.k        = 0;
    xSPM.df       = [];
    xSPM.u        = threshold;
    xSPM.STAT     = '';
    xSPM.VOX      = vox;        
    xSPM.thres    = threshold;
    xSPM.vols     = post;
    xSPM.scale    = odds;
    xSPM.units    = {'mm'  'mm'  'mm'};
    xSPM.Ps       = [];
    xSPM.S        = 0;
    xSPM.z_ps     = z_ps;
    BMS.xSPM      = xSPM;
    
    % Display results
    spm_bms_display(BMS,'Init');
    
else
    
    % No voxels for selected threshold
    disp('No voxels above threshold!');
    uicontrol(Finter,'Style','Text','String','No voxels above threshold!',...
            'Position',[020 010 200 020].*WS,...
            'FontAngle','Italic',...
            'FontSize',FS(12),...
            'HorizontalAlignment','Left',...
            'ForegroundColor','w')
    do_results = spm_input('Return to Results?', '+1', 'y/n', [1,0], 2);
    
    % Go back to BMS Maps (Results)
    if do_results  
        
        % Return to display routine 
        spm_run_bms_vis();
        
    else
        
        spm('Clear',Finter); % Clear Finter window
        return               % Finish 
        
    end
    
    
end