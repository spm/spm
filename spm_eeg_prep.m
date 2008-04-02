function D = spm_eeg_prep(S)

D = S.D;

switch S.task
    case 'settype'
        D = chantype(D, S.ind, S.type);
    case 'loadtemplate'
        template = load(S.P); % must contain Cpos, Cnames

        [sel1, sel2] = spm_match_str(lower(D.chanlabels), lower(template.Cnames));

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

            D = coor2D(D, sel1, num2cell(template.Cpos(:, sel2)));
        end
    case 'loadeegsens'
        switch S.source
            case 'mat'
                senspos = load(S.sensfile{1});
                name    = fieldnames(senspos);
                senspos = getfield(senspos,name{1});

                fiducials = load(S.sensfile{2});
                name    = fieldnames(fiducials);
                fiducials = getfield(fiducials,name{1});         
                
                label = chanlabels(D, strmatch('EEG', D.chantype, 'exact'));
            case 'filpolhemus'
                [fiducials, senspos] = spm_eeg_inv_ReadPolhemus(S.sensfile);
                label = chanlabels(D, strmatch('EEG', D.chantype, 'exact'));
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