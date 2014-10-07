function [Ds, D] = spm_eeg_inv_extract(D)
% Exports source activity using the MAP projector
% FORMAT [Ds] = spm_eeg_inv_extract(D)
% Requires:
%
%     D.inv{i}.source.XYZ   - (n x 3) matrix of MNI coordinates
%
% Optional:
%
%     D.inv{i}.source.rad   - radius (mm) of VOIs (default 5 mm)
%     D.inv{i}.source.label - label(s) for sources (cell array)
%     D.inv{i}.source.fname - output file name
%     D.inv{i}.source.type  - output type ('evoked'/'trials')
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Vladimir Litvak, Laurence Hunt, Karl Friston
% $Id: spm_eeg_inv_extract.m 6231 2014-10-07 13:42:16Z vladimir $
 
% SPM data structure
%==========================================================================
try
    inv = D.inv{D.val};
catch
    inv = D.inv{end};
    D.val = numel(D.inv);
end
 
Ds = [];
 
% defaults
%--------------------------------------------------------------------------
try, XYZ   = inv.source.XYZ;   catch,  return;                          end
try, rad   = inv.source.rad;   catch,  rad   = 5;                       end
try, label = inv.source.label; catch,  label = {};                      end
try
    fname  = inv.source.fname; 
catch
    fname  = fullfile(D.path, ['i' D.fname]); 
end
try
    type   =  inv.source.type;
catch
    if isequal(D.type, 'evoked')
        type = 'evoked';
    else
        type = 'trials';
    end
end
 
if numel(label)<size(XYZ, 1)
    for i = (numel(label)+1):size(XYZ, 1)
        label{i} = ['Source ' num2str(i)];
    end
end
 
% find relevant vertices
%==========================================================================
vert  = inv.mesh.tess_mni.vert(inv.inverse.Is, :);        % vertices
Ns    = size(XYZ, 1);                        % number of sources
svert = {};
for i = 1:Ns
    dist = sqrt(sum([vert(:,1) - XYZ(i,1), ...
                     vert(:,2) - XYZ(i,2), ...
                     vert(:,3) - XYZ(i,3)].^2, 2));
    if rad > 0
        svert{i} = find(dist < rad);
    else
        [junk,svert{i}] = min(dist);
        XYZ(i, :) = vert(svert{i}, :);
    end
end
 

% report
%--------------------------------------------------------------------------
Iy = find(cellfun('isempty', svert));
if ~isempty(Iy)
    disp(['No proximal vertices found for source(s) ' num2str(Iy)]);
    disp('These sources will be discarded');
    
    svert(Iy) = [];
    label(Iy) = [];
    Ns        = numel(label);
    
    if Ns == 0
        disp('No valid sources - Exiting');
        return;
    end
end
 
js      = spm_vec(svert);
 
% get info from inv
%--------------------------------------------------------------------------
scale   = inv.inverse.scale;
J       = inv.inverse.J;  
U       = inv.inverse.U; % spatial projector (contains multiple cells if multimodal)
T       = inv.inverse.T; % S in Friston 2008 - temporal projector
TT      = T*T';
M       = inv.inverse.M(js, :);
Ic      = inv.inverse.Ic;
It      = inv.inverse.It;
Np      = length(It);
 
try
    trial = inv.inverse.trials;
catch
    trial = D.condlist;
end
 
 
% get source data
%==========================================================================
MYi    = 1;
switch(type)
    case 'evoked'
        Ne     = length(trial);
        MY     = zeros(size(M,1), Ne*Np);   
        
        for i = 1:Ne
            MY(:, MYi:(MYi+Np-1)) = J{i}(js,:)*T';
            MYi                   = MYi + Np;
        end
        
        clabel = trial;
    case 'trials'
        Ne     = length(D.indtrial(trial, 'GOOD'));
        MY     = zeros(size(M,1), Ne*Np);
        clabel = {};
        
        for i = 1:numel(trial)
            
            c      = D.indtrial(trial{i}, 'GOOD');
            clabel = [clabel D.conditions(c)];
            
            % conditional expectation of contrast (J*W) and its energy
            %--------------------------------------------------------------
            Nt    = length(c);
            spm_progress_bar('Init',Nt,sprintf('extracting data condition %d',i),'trials');
            
            for j = 1:Nt
                
                if ~strcmp(D.modality(1,1), 'Multimodal')
                    
                    % unimodal data
                    %------------------------------------------------------
                    Y     = D(Ic{1},It,c(j));
                    Y     = U{1}*Y*scale;
                    
                else
                    
                    % multimodal data
                    %------------------------------------------------------
                    for k = 1:length(U)
                        Y       = D(Ic{k},It,c(j));
                        UY{k,1} = U{k}*Y*scale(k);
                    end
                    Y = spm_cat(UY);
                end
                
                MY(:, MYi:(MYi+Np-1)) = M*Y;
                MYi                   = MYi+Np;
                
                spm_progress_bar('Set',j)
            end
            spm_progress_bar('Clear')
        end
end
 
% compute regional response in terms of first eigenvariate
%==========================================================================
iS    = [0 cumsum(cellfun('length', svert))];
Y     = zeros(Ns, size(MY, 2));
 
spm_progress_bar('Init', Ns, 'extracting eigenvariates', 'sources');
 
for i = 1:Ns
    y     = MY((iS(i)+1):iS(i+1),:)';
    
    [m,n]   = size(y);
    if n>1
        if m > n
            [v,s,v] = svd(y'*y);
            s       = diag(s);
            v       = v(:,1);
            u       = y*v/sqrt(s(1));
        elseif m>1
            [u,s,u] = svd(y*y');
            s       = diag(s);
            u       = u(:,1);
            v       = y'*u/sqrt(s(1));
        end
        d       = sign(sum(v));
        u       = u*d;
        v       = v*d;
        Y(i, :) = u'*sqrt(s(1)/n);
    else
        Y(i, :) = y';
    end
    
    spm_progress_bar('Set',i);
end

spm_progress_bar('Clear')
 
% create source dataset
%-----------------------------------------------------------------------
Ds = clone(D, fname, [Ns Np Ne], 3);
Ds = chanlabels(Ds, 1:Ns, label);
Ds = timeonset(Ds, D.time(It(1)));
Ds = chantype(Ds, 1:Ns, 'LFP');
Ds = conditions(Ds, 1:Ne, clabel);
Ds = sensors(Ds, 'MEG', []);
Ds = sensors(Ds, 'EEG', []);
Ds = fiducials(Ds, []);
Ds = rmfield(Ds, fieldnames(Ds));
Ds = condlist(Ds, trial);
 
spm_progress_bar('Init', Ne, 'writing out data', 'trials');
for i = 1:Ne
    Ds(:,:,i) = Y(:,((i - 1)*Np + 1):i*Np);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');
 
% store results in the original inv
%-----------------------------------------------------------------------
D.inv{D.val}.source.type  = type;
D.inv{D.val}.source.XYZ   = XYZ;
D.inv{D.val}.source.label = label;
D.inv{D.val}.source.rad   = rad;
D.inv{D.val}.source.fname = fname;
 
% save headers
%-----------------------------------------------------------------------
save(Ds);


