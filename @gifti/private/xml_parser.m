function tree = xml_parser(xml)
% XML Parser
% FORMAT tree = xml_parser(xml)
% xml         - XML-encoded string or filename of an XML document
% tree        - struct array representation of the XML document
%__________________________________________________________________________
%
% This C-MEX file is a wrapper around yxml:
%   https://dev.yorhel.nl/yxml
% by Yoran Heling:
%   https://yorhel.nl/
%
% A pure MATLAB implementation of a similar XML parser is available at:
%   https://github.com/gllmflndn/xmltree
%__________________________________________________________________________
%
% The tree representation of the XML document is stores as a struct array
% with fields:
%  - type:       'element' or 'chardata'
%  - value:      tag name of an 'element' or content of a 'chardata'
%  - attributes: key/value struct array of element's attributes
%  - children:   array of uids of element's children
%  - uid:        unique identifier (index in the struct array)
%  - parent:     uid of parent ([] if root)
%
% This corresponds to an XML string of the sort:
% <value key="value">value</value>
%
% Processing instructions and comments are not reported.
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging


error('A compiled version of "xml_parser" is not available for your platform.');
