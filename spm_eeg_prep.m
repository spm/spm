function D = spm_eeg_prep(S)
% spm_eeg_prep function performs several tasks
% for preparation of converted MEEG data for further analysis
% FORMAT spm_eeg_prep(S)
%   S - configuration struct (obligatory)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep.m 1336 2008-04-09 14:57:45Z vladimir $

D = S.D;

switch S.task
    case 'settype'
        D = chantype(D, S.ind, S.type);
    case {'loadtemplate', 'setcoor2d'}
        if strcmp(S.task, 'loadtemplate')
            template = load(S.P); % must contain Cpos, Cnames
            xy = template.Cpos;
            label = template.Cnames;
        else
            xy = S.xy;
            label = S.label;
        end

        [sel1, sel2] = spm_match_str(lower(D.chanlabels), lower(label));

        if ~isempty(sel1)

            megind = strmatch('MEG', chantype(D), 'exact');
            eegind = strmatch('EEG', chantype(D), 'exact');

            if ~isempty(intersect(megind, sel1)) && ~isempty(setdiff(megind, sel1))
                error('2D locations not found for all MEG channels');
            end

            if ~isempty(intersect(eegind, sel1)) && ~isempty(setdiff(eegind, sel1))
                warning(['2D locations not found for all EEG channels, changing type of channels', ...
                    num2str(setdiff(eegind, sel1)) ' to ''Other''']);

                D = chantype(D, setdiff(eegind, sel1), 'Other');
            end
            
            if any(any(coor2D(D, sel1) - xy(:, sel2)))
                D = coor2D(D, sel1, num2cell(xy(:, sel2)));
            end
        end
    case 'loadeegsens'
        switch S.source
            case 'mat'
                fiducials = load(S.sensfile{1});
                name    = fieldnames(fiducials);
                fiducials = getfield(fiducials,name{1});

                senspos = load(S.sensfile{2});
                name    = fieldnames(senspos);
                senspos = getfield(senspos,name{1});
                
                label = chanlabels(D, sort(strmatch('EEG', D.chantype, 'exact')));

            case 'filpolhemus'
                [fiducials, senspos] = spm_eeg_inv_ReadPolhemus(S.sensfile); 
                
                label = chanlabels(D, sort(strmatch('EEG', D.chantype, 'exact')));
                
            case 'locfile'
                elec = read_sens(S.sensfile);
                
                [junk fidind] = spm_match_str({'fidnz', 'fidt9', 'fidt10'}, lower(elec.label));
                
                if length(fidind) ~= 3
                    error('Locations file should contain 3 fiducials labeled ''fidnz'', ''fidt9'' and ''fidt10''');
                end                                                               
                
                sensind = setdiff(1:length(elec.label), fidind);
                
                fiducials = elec.pnt(fidind, :);
                label = elec.label(sensind);
                senspos = elec.pnt(sensind, :);
        end        

        if size(senspos, 1) ~= length(label)
            error('To read sensor positions without labels the numbers of sensors and EEG channels should match.');
        end

        elec = [];
        elec.pnt = senspos;
        elec.label = label;
        elec.fid = fiducials;

        D = elec2sens(D, elec);
    case 'headshape'
        switch S.source
            case 'mat'
                headshape = load(S.headshapefile);
                name    = fieldnames(headshape);
                headshape = getfield(headshape,name{1});
            case 'filpolhemus'
                [fiducials, headshape] = spm_eeg_inv_ReadPolhemus(S.headshapefile);
        end

        D.headshape = headshape;
    case 'coregister'
        D = D.sensorcoreg;
end