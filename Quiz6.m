clear; clc;

A = [0.4264, 0.5222, -0.7385;
     -0.8165, 0.5774, 0;
     -0.5345, -0.6547, -0.5345];

% Checking columns are unit vectors
norm(A(1,:))
norm(A(2,:))
norm(A(3,:))

% Checking if rows are unit vectors
A = A'
norm(A(1,:))
norm(A(2,:))
norm(A(3,:))

% Checking orthogonality
inv(A) == A'

det(A) == 1

%%
clear; clc;

A = [0.4264, 0.5222, -0.7385; -0.8165, 0.5774, 0 ; -0.5345, -0.6547, -0.5345 ]