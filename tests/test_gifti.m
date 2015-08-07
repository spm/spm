function tests = test_gifti
% Unit Tests for gifti
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% $Id: test_gifti.m 6516 2015-08-07 17:28:33Z guillaume $

tests = functiontests(localfunctions);


function test_gifti_constructor(testCase)
import matlab.unittest.constraints.*
g1 = gifti;
g2 = gifti(fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii'));
g3 = gifti(g2);
s  = struct(g3);
g4 = gifti(s);
g5 = gifti(rand(16,1));

testCase.verifyThat(g1, IsOfClass('gifti'));
testCase.verifyThat(fieldnames(g1), IsEmpty);

testCase.verifyThat(g2, IsOfClass('gifti'));
testCase.verifyThat(s, HasField('faces'));
testCase.verifyThat(s, HasField('vertices'));
testCase.verifyThat(s, HasField('mat'));

testCase.verifyThat(g3, IsOfClass('gifti'));
testCase.verifyThat(g3, IsEqualTo(g2));

testCase.verifyThat(s, IsOfClass('struct'));
testCase.verifyThat(s, HasField('faces'));
testCase.verifyThat(s, HasField('vertices'));
testCase.verifyThat(s, HasField('mat'));

testCase.verifyThat(g4, IsOfClass('gifti'));
testCase.verifyThat(struct(g4), IsEqualTo(s));

testCase.verifyThat(g5, IsOfClass('gifti'));
testCase.verifyThat(g5, HasField('cdata'));


function test_gifti_accessor(testCase)
import matlab.unittest.constraints.*
g1 = gifti(fullfile(spm('Dir'),'canonical','cortex_8196.surf.gii'));
s.faces    = g1.faces;
s.vertices = g1.vertices;
s.mat      = g1.mat;
testCase.verifyEqual(s.faces(1000,:), g1.faces(1000,:));
testCase.verifyEqual(s.vertices(1000,:), g1.vertices(1000,:));
testCase.verifyEqual(s.mat(:,4), g1.mat(:,4));

cdata = single(rand(16,1));
g2 = gifti(cdata);
testCase.verifyEqual(cdata,g2.cdata);
testCase.verifyEqual(cdata(8),g2.cdata(8));
testCase.verifyEqual(cdata(8,1),g2.cdata(8,1));


function test_gifti_mutator(testCase)
import matlab.unittest.constraints.*
g = gifti(fullfile(spm('Dir'),'canonical','cortex_5124.surf.gii'));
g.mat = eye(4);
g.mat(logical(eye(4))) = 2;
g.mat(4,4) = 1;
testCase.verifyEqual(g.mat, diag([2 2 2 1]));

faces = g.faces;
faces = faces(:,[2 1 3]);
g.faces = faces;
g.faces(64,:) = [1 2 3];
g.faces(:,1) = faces(:,1);
testCase.verifyEqual(g.faces(64,:), [faces(64,1) 2 3]);

vertices = g.vertices;
vertices = vertices(:,[2 1 3]);
g.vertices = vertices;
g.vertices(64,:) = zeros(1,3);
g.vertices(:,2:3) = ones(size(vertices,1),2);
testCase.verifyEqual(g.vertices(64,:), single([0 1 1]));

g.cdata = rand(size(g.vertices,1),1);
g.cdata(1) = pi;
g.cdata(2:end) = exp(1);
testCase.verifyEqual(g.cdata(1), single(pi));


function test_gifti_export(testCase)
import matlab.unittest.constraints.*
mri = load('mri');
g = gifti(isosurface(smooth3(squeeze(mri.D)),5));
s = export(g);
testCase.verifyThat(s, HasField('vertices'));
testCase.verifyThat(s, HasField('faces'));
s = export(g,'matlab');
testCase.verifyThat(s, HasField('vertices'));
testCase.verifyThat(s, HasField('faces'));
s = export(g,'patch');
testCase.verifyThat(s, HasField('vertices'));
testCase.verifyThat(s, HasField('faces'));
testCase.verifyThat(s, ~HasField('mat'));
s = export(g,'fieldtrip');
testCase.verifyThat(s, HasField('pnt'));
testCase.verifyThat(s, HasField('tri'));
s = export(g,'spm');
testCase.verifyThat(s, HasField('vert'));
testCase.verifyThat(s, HasField('face'));


function test_gifti_load(testCase)
import matlab.unittest.constraints.*
d = fullfile(spm('Dir'),'canonical');
files = dir(fullfile(d,'*.gii'));
for i=1:numel(files)
    g = gifti(fullfile(d,files(i).name));
    testCase.verifyThat(evalc('g'), ~IsEmpty); % check display()
end


function test_gifti_save(testCase)
mri = load('mri');
g = gifti(isosurface(smooth3(squeeze(mri.D)),5));
g.cdata = rand(size(g.vertices,1),1);
basename = tempname;
file = [basename '.gii'];
save(g,file,'ASCII');
g = gifti(file);
delete(file);
save(g,file,'Base64Binary');
g = gifti(file);
delete(file);
save(g,file,'GZipBase64Binary');
g = gifti(file);
delete(file);
save(g,file,'ExternalFileBinary');
g = gifti(file);
delete(file);
delete([basename '.dat']);
