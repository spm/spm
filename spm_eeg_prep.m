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
end