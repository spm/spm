function [eeg,t] = read_eep_trial(filename, triggernumber, interval);
%
% READ_EEP_TRIAL reads a data from an EEProbe *.cnt file
%
% [eeg,t] = read_eep_trial(filename, triggernumber, interval);
%
% interval = [ -2, 5] reads a window of -2 to 5 seconds around
% the given trigger number
%
% Script returns eeg data structure: it contains the data, labels etc
% for the given interval
%
% t is the time in milliseconds of the trial (at the trigger, i.e., stimulus-time)
%
% See also READ_EEP_TRG, READ_EEP_REJ, READ_EEP_AVR
%
% ANT Software BV, Netherlands, www.ant-software.nl

% read short piece of data to get sampling rate, channels etc
eeg = read_eep_cnt([filename '.cnt'],100,101);
samplerate = eeg.rate;
eeglabels = eeg.label;

% read trigger information from external trigger file
trg = read_eep_trg([filename '.trg']);
if triggernumber > length(trg)
	error('Invalid trigger number, trigger does not exist!'); return
else
   trg = trg(triggernumber);
end

% compute interval to extract from file in milliseconds
sample1 = trg.time + interval(1)*1000;
sample2 = trg.time + interval(2)*1000;

% convert interval to samples
sample1 = sample1 /1000 * samplerate + 1;
sample2 = sample2 /1000 * samplerate + 1;

% read data from file
eeg = read_eep_cnt([filename '.cnt'],sample1,sample2);
t = trg.time;

