function res = bf_group_GALA(BF, S)
% Computes Minimum Norm projectors
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, using the code from Matti Stenroos and Olaf Hauk
% http://imaging.mrc-cbu.cam.ac.uk/meg/AnalyzingData/DeFleCT_SpatialFiltering_Tools
%
% Please cite:
% Hauk O, Stenroos M.
% A framework for the design of flexible cross-talk functions for spatial filtering of EEG/MEG data: DeFleCT.
% Human Brain Mapping 2013
% $Id: bf_group_GALA.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    iter = cfg_entry;
    iter.tag = 'iter';
    iter.name = 'Number of iteration';
    iter.strtype = 'r';
    iter.num = [1 1];
    iter.val = {3};
    iter.help = {'Number of iteration'};
    
   
    GALA      = cfg_branch;
    GALA.tag  = 'GALA';
    GALA.name = 'GALA';
    GALA.val  = {iter};
    GALA.help = {'bla bla bla'};
    res = GALA;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

S.modality = 'MEGPLANAR'; % temp

Nl = length(BF);
S.Nl = Nl;

[J, St] = GALA_invert(BF,S);

Nd = length(J)/Nl;
res = cell(1, numel(BF));
for p=1:Nl
    BFp = load(BF{p});
    % extract sunbject's J
    Jp = J(1+(p-1)*Nd:Nd+(p-1)*Nd,:);
    BFp.inverse.(S.modality).J = Jp;
    BFp.inverse.(S.modality).S = St;
    BFp.inverse.(S.modality).channels = BFp.sources.channels.MEG; % temp
    
    outdir = spm_file(BF{p}, 'fpath');
    bf_save_path(BFp,fullfile(outdir, 'BF.mat'));
    res(p) = {fullfile(outdir, 'BF.mat')};
end

end


