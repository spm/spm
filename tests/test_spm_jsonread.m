function tests = test_spm_jsonread
% Unit Tests for spm_jsonread
%__________________________________________________________________________
% Copyright (C) 2015-2017 Wellcome Trust Centre for Neuroimaging

% $Id: test_spm_jsonread.m 7077 2017-05-23 09:13:10Z guillaume $

tests = functiontests(localfunctions);


function test_jsonread_from_string_1(testCase)
json = '["one", "two", "three"]';

exp = {'one';'two';'three'};
act = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

function test_jsonread_from_string_2(testCase)
json = '{"Width":800,"Height":600,"Title":"View from the 15th Floor","Animated":false,"IDs":[116,943,234,38793]}';

exp = struct('Width',800,'Height',600,'Title','View from the 15th Floor','Animated',false,'IDs',[116;943;234;38793]);
act = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

function test_jsonread_all_types(testCase)

%   JSON Data Type            | MATLAB Data Type

%   null, in numeric arrays   | NaN
json = '[1, 2, null, 4]';
exp  = [1; 2; NaN; 4];
act  = spm_jsonread(json);
%testCase.verifyTrue(isequaln(exp, act));

%   null, in nonnumeric arrays| empty double []
json = '{"null": null}';
exp  = struct('null',[]);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Boolean                   | scalar logical
json = '{"logical": false}';
exp  = struct('logical',false);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '{"logical": true}';
exp  = struct('logical',true);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Number                    | scalar double
json = '{"number": 3.14}';
exp  = struct('number',3.14);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   String                    | character vector
json = '{"string": "string"}';
exp  = struct('string','string');
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Object (In JSON, object   | scalar structure
%    means an unordered set   |  (Names are made
%    of name-value pairs.)    |  valid.)
json = '{"object": {"field1": 1, "field-2": 2, "3field": 3}}';
exp  = struct('object',struct('field1',1,'field_2',2,'x3field',3));
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '{"object": {"field 1": 1, "field two": 2, "  field  Three  ": 3}}';
exp  = struct('object',struct('field1',1,'fieldTwo',2,'fieldThree',3));
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '{"object": {"": 1}}';
exp  = struct('object',struct('x',1));
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array, when elements are  | cell array
%    of different data types  |
json = '{"array": ["a", 1]}';
exp  = struct('array',{{'a';1}});
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of booleans         | logical array
json = '{"logical_array": [true, false]}';
exp  = struct('logical_array',[true;false]);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of numbers          | double array
json = '{"number_array": [1, 2, 3, 5, 8, 13]}';
exp  = struct('number_array',[1; 2; 3; 5; 8; 13]);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of strings          | cellstr 
json = '{"cellstr": ["Statistical","Parametric","Mapping"]}';
exp  = struct('cellstr',{{'Statistical';'Parametric';'Mapping'}});
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of objects, when    | structure array
%    all objects have the     |
%    same set of names        |
json = '{"structarray": [{"a":1,"b":2},{"a":3,"b":4}]}';
exp  = struct('structarray',struct('a',{1;3},'b',{2;4}));
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of objects, when    | cell array of
%    objects have different   | scalar structures
%    names                    |
json = '{"cellarray": [{"a":1,"b":2},{"a":3,"c":4}]}';
exp  = struct('cellarray',{{struct('a',1,'b',2);struct('a',3,'c',4)}});
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

%   Array of objects, when    | cell array of
%    all objects have the     |scalar structures
%     same set of names        |
%    but different order
json = '{"structarray": [{"a":1,"b":2},{"b":3,"a":4}]}';
exp  = struct('structarray',{{struct('a',1,'b',2);struct('b',3,'a',4)}});
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

% empty struct
json = '{"a":"aa","b":{},"c":"cc"}';
exp = struct('a','aa','b',struct(),'c','cc');
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

% empty array
json = '{"a":"aa","b":[],"c":"cc"}';
exp = struct('a','aa','b',[],'c','cc');
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

% NaN and Inf
json = '[1,NaN,Inf]';
exp  = [1;NaN;Inf];
act  = spm_jsonread(json);
testCase.verifyTrue(isequaln(exp, act));

% numbers and booleans
json = '[1,NaN,Inf,true]';
exp  = {1;NaN;Inf;true};
act  = spm_jsonread(json);
testCase.verifyTrue(isequaln(exp, act));

% numbers and arrays of numbers
json = '[1,2,[3,4]]';
exp  = {1;2;[3;4]};
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '[[1,2]]';
exp  = [1 2];
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '[[1,2],[3,4]]';
exp  = [1 2;3 4]; 
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '[[1,2],[3,4,5]]';
exp  = {[1;2];[3;4;5]}; 
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));

json = '[[[1,2],[3,4]],[[5,6],[7,8]]]';
exp  = cat(3,[1,3;5,7],[2,4;6,8]);
act  = spm_jsonread(json);
testCase.verifyTrue(isequal(exp, act));
