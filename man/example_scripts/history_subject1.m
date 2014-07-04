spm('defaults', 'eeg');

S = [];
S.dataset = 'subject1.bdf';
S.mode = 'continuous';
S.channels = {'EEG'};
S.eventpadding = 0;
S.blocksize = 3276800;
S.checkboundary = 1;
S.saveorigheader = 0;
S.outfile = 'spmeeg_subject1';
S.timewin = [];
S.conditionlabels = {'Undefined'};
S.inputformat = [];
D = spm_eeg_convert(S);


S = [];
S.D = 'spmeeg_subject1.mat';
S.mode = 'write';
S.blocksize = 655360;
S.prefix = 'M';
S.montage = 'avref_montage.mat';
S.keepothers = 0;
S.keepsensors = 1;
S.updatehistory = 1;
D = spm_eeg_montage(S);


S = [];
S.D = 'Mspmeeg_subject1.mat';
S.type = 'butterworth';
S.band = 'high';
S.freq = 0.1;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


S = [];
S.D = 'fMspmeeg_subject1.mat';
S.fsample_new = 200;
S.prefix = 'd';
D = spm_eeg_downsample(S);


S = [];
S.D = 'dfMspmeeg_subject1.mat';
S.type = 'butterworth';
S.band = 'low';
S.freq = 30;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
D = spm_eeg_filter(S);


S = [];
S.D = 'fdfMspmeeg_subject1.mat';
S.timewin = [-100 400];
S.trialdef(1).conditionlabel = 'std';
S.trialdef(1).eventtype = 'STATUS';
S.trialdef(1).eventvalue = 1;
S.trialdef(1).trlshift = 0;
S.trialdef(2).conditionlabel = 'odd';
S.trialdef(2).eventtype = 'STATUS';
S.trialdef(2).eventvalue = 3;
S.trialdef(2).trlshift = 0;
S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
D = spm_eeg_epochs(S);


S = [];
S.D = 'efdfMspmeeg_subject1.mat';
S.mode = 'reject';
S.badchanthresh = 0.2;
S.methods.channels = {'EEG'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 80;
S.methods.settings.excwin = 1000;
S.append = true;
S.prefix = 'a';
D = spm_eeg_artefact(S);


S = [];
S.D = 'aefdfMspmeeg_subject1.mat';
S.robust.ks = 3;
S.robust.bycondition = false;
S.robust.savew = false;
S.robust.removebad = true;
S.circularise = false;
S.prefix = 'm';
D = spm_eeg_average(S);


