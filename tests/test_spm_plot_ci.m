function tests = test_spm_plot_ci
% Unit Tests for spm_plot_ci
% Ensures that all the different plot types run without error
%__________________________________________________________________________

% Copyright (C) 2016-2022 Wellcome Centre for Human Neuroimaging


tests = functiontests(localfunctions);

function test_barchart(testCase)

E = (1:5)';
C = [1 1 1 1/4 1/4]';

do_plot(E,C);

function test_barchart_exp(testCase)

E = (1:5)';
C = [1 1 1 1/4 1/4]';
do_plot(E,C,'exp');

function test_grouped_barchart(testCase)

E = [1 2 3 4 5
     6 7 8 9 10]';
C = [1   1   1 1/4 1/4;
     1/4 1/4 1 1   1]';
do_plot(E,C);

function test_grouped_barchart_exp(testCase)

E = [1 2 3 4 5
     6 7 8 9 10]';
C = [1   1   1 1/4 1/4;
     1/4 1/4 1 1   1]';
do_plot(E,C,'exp');

function test_confidence_region(testCase)

E = [3 4];
C = [1 1/4];
do_plot(E,C);

function test_linechart(testCase)

rng(1);
E = randn(3,8); % At least 8 columns (i.e., measurements)
C = exp(randn(3,8));
do_plot(E,C);

function do_plot(E,C,str)
% Run spm_plot_ci within a figure
f=figure;
if nargin < 3
    spm_plot_ci(E,C);
else
    spm_plot_ci(E,C,[],[],str);
end
close(f);
