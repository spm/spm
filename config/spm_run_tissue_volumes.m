function out = spm_run_tissue_volumes(cmd, job)
% SPM job execution function for Tissue Volumes
%
% See also: spm_cfg_tissue_volumes, spm_summarise
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Ged Ridgway
% $Id: spm_run_tissue_volumes.m 5800 2013-12-10 18:33:15Z guillaume $

switch lower(cmd)
    %----------------------------------------------------------------------    
    case 'exec'
        
        %%
        mat = job.matfiles;
        T   = job.tmax;
        msk = char(job.mask);
        outf = job.outf;
        if isempty(msk)
            msk = 'all';
        end
        
        %%
        N = numel(mat);
        vol = nan(N, T);
        for n = 1:N
            res = load(mat{n});
            
            % check for previously computed volumes (to save time if there)
            if isfield(res, 'volumes') && strcmp(res.volumes.mask, msk) ...
                    && numel(res.volumes.litres) >= T
                vol(n, :) = res.volumes.litres(1:T);
                continue
            end
            
            % determine number of tissue classes Kb
            if isfield(res, 'mg'),
                Kb  = max(res.lkp);
            else
                Kb  = size(res.intensity(1).lik, 2);
            end
            tc  = false(Kb, 4);
            
            % look for existing mwc files
            fnm = res.image.fname;
            if ~exist(fnm, 'file')
                error('Original image no longer found at:\n%s\n', fnm)
            end
            mwc = cell(T, 1);
            for t = 1:T
                mwc{t} = spm_file(fnm, 'prefix', ['mwc' num2str(t)],'ext','nii');
                if ~exist(mwc{t}, 'file')
                    tc(t, 4) = true; % i.e. need to produce this mwc
                end
            end
            
            % produce mwc files if required
            if any(tc(:, 4))
                spm_preproc_write8(res, tc);
            end
            
            % compute tissue volumes
            vol(n, :) = spm_summarise(mwc, msk, 'litres');
            
            % add to mat file for future reuse
            volumes.litres  = vol(n, :);
            volumes.mask    = msk;
            save(mat{n}, 'volumes', '-append')
            
            % if mwc newly created above, delete now
            for t = 1:T
                if tc(t, 4), spm_unlink(mwc{t}), end
            end
        end
        
        %% Put into output structure for use with dependencies
        for t = 1:T
            out.(sprintf('vol%d', t)) = vol(:, t);
        end
        out.vol_sum = sum(vol, 2); % (total intracranial volume if T=1:3)
        
        %% Optionally save in CSV format
        if ~isempty(outf)
            [pth, nam, ext] = spm_fileparts(outf);
            if isempty(ext), ext = '.csv'; end
            fnm = fullfile(pth, [nam ext]);
            fid = fopen(fnm, 'wt');
            if fid < 0, error('Failed to open %s\n', fnm); end
            delim = ',';
            fprintf(fid, 'File');
            fprintf(fid, [delim 'Volume%d'], 1:T);
            for n = 1:N
                fprintf(fid, '\n''%s''', mat{n});
                fprintf(fid, [delim '%d'], vol(n, :));
            end
            fprintf(fid, '\n');
            fclose(fid);
        end
        
        %% Display in command window
        fprintf('\nSegmentation files:\n');
        fprintf('\t%s\n', mat{:});
        fprintf('\nVolumes (litres):\n');
        disp(vol);
        %------------------------------------------------------------------
    case 'vout'
        try
            T = job.tmax;
        catch
            T = 3;
        end
        out(T+1) = cfg_dep;
        for t = 1:T
            out(t).sname        = num2str(t);
            out(t).src_output   = substruct('.', sprintf('vol%d', t));
            out(t).tgt_spec     = cfg_findspec({
                {'strtype','e', 'strtype','r'}
                });
        end
        out(T+1).sname      = 'Sum';
        out(T+1).src_output = substruct('.', 'vol_sum');
        out(T+1).tgt_spec   = cfg_findspec({
            {'strtype','e', 'strtype','r'}
            });
end
