% Tests the internal consistency of the kernel rotation functions.
%
% ilya.kuprov@oerc.ox.ac.uk

function rotations_test()

% Generate a random symmetric traceless 3x3 matrix
A=rand(3); A=A+A'; A=A-trace(A)/3;

% Generate a random set of Euler angles
eulers=rand(1,3); eulers=eulers.*[2*pi pi 2*pi];

%% Test 1: euler2dcm, euler2wigner, mat2sphten

% DCM rotation followed by a transformation into irreducible components
[~,~,rank2_a]=mat2sphten(euler2dcm(eulers)*A*euler2dcm(eulers)');

% Transformation into irreducible components followed by a Wigner rotation
[~,~,rank2]=mat2sphten(A); rank2_b=(euler2wigner(eulers)*rank2')';

% Check the difference
if norm(rank2_a-rank2_b)<1e-10
    disp('Test 1 passed.');
else
    disp('Test 1 inconsistency detected');
end

%% Test 2: euler2dcm, dcm2euler

% Transforms Euler angles into DCM
R=euler2dcm(eulers);

% Transform the DCM back into Euler angles
new_eulers=dcm2euler(R);

% Check the difference
if norm(eulers-new_eulers)<1e-3
    disp('Test 2 passed.');
else
    disp('Test 2 inconsistency detected');
end

%% Test 3: euler2dcm, euler2wigner, dcm2wigner

% Transforms Euler angles into DCM, then DCM to Wigner matrix
W_a=dcm2wigner(euler2dcm(eulers));

% Transform Euler angles to Wigner matrix
W_b=euler2wigner(eulers);

% Check the difference
if norm(W_a-W_b)<1e-10
    disp('Test 3 passed.');
else
    disp('Test 3 inconsistency detected');
end

%% Test 4: mat2sphten, sphten2mat

% Transform the matrix into spherical tensors
[rank0,rank1,rank2]=mat2sphten(A);

% Transform back
B=sphten2mat(rank0,rank1,rank2);

% Check the difference
if norm(A-B)<1e-10
    disp('Test 4 passed.');
else
    disp('Test 4 inconsistency detected');
end


end