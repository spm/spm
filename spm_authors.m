function [current, previous] = spm_authors
% Return list of SPM coauthors
% FORMAT [current, previous] = spm_authors
% current  - cell array of coauthors of the current SPM release
% previous - cell array of former SPM coauthors
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% SPM
% $Id: spm_authors.m 4025 2010-07-29 11:10:15Z guillaume $


current = {...
'John Ashburner'
'Gareth Barnes'
'Chun-Chuan Chen'
'Justin Chumbley'
'Jean Daunizeau'
'Guillaume Flandin'
'Karl Friston'
'Darren Gitelman'
'Volkmar Glauche'
'Lee Harrison'
'Rik Henson'
'Chloe Hutton'
'Maria Joao Rosa'
'Stefan Kiebel'
'James Kilner'
'Vladimir Litvak'
'Rosalyn Moran'
'Tom Nichols'
'Robert Oostenveld'
'Will Penny'
'Christophe Phillips'
'Ged Ridgway'
'Klaas Enno Stephan'
};

previous = {...
'Jesper Andersson'
'Matthew Brett'
'Christian Buechel'
'Jon Heather'
'Andrew Holmes'
'Jeremie Mattout'
'Jean-Baptiste Poline'
'Keith Worsley'
};
